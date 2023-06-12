from functools import cached_property
from typing import Dict, List, Optional

from abnumber.chain import Chain, Alignment, Position
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from dataclasses import dataclass
from flask import current_app
import numpy as np
from numpy.random import default_rng
import pandas as pd
import sapiens

from biophi.common.utils.formatting import get_valid_filename


@dataclass
class HumanizationParams:
    method = None
    scheme: str = 'kabat'
    cdr_definition: str = 'kabat'

    def get_export_name(self):
        raise NotImplementedError()


@dataclass
class SapiensHumanizationParams(HumanizationParams):

    method = 'sapiens'
    model_version: str = 'latest'
    humanize_cdrs: bool = False
    backmutate_vernier: bool = False
    iterations: int = 1

    def get_export_name(self):
        name = f'Sapiens_{self.iterations}iter_'
        if self.humanize_cdrs:
            name += f'humanized_{self.cdr_definition.title()}_CDRs_'
        else:
            name += f'parental_{self.cdr_definition.title()}_CDRs_'
        return name


@dataclass
class CDRGraftingHumanizationParams(HumanizationParams):
    method = 'cdr_grafting'
    heavy_v_germline: str = 'auto'
    light_v_germline: str = 'auto'
    backmutate_vernier: bool = True
    sapiens_iterations: int = 0

    def get_export_name(self):
        if self.backmutate_vernier:
            name = f'{self.cdr_definition.title()}_Vernier_graft_'
        else:
            name = f'{self.cdr_definition.title()}_Straight_graft_'

        if self.sapiens_iterations:
            name += f'Sapiens_{self.sapiens_iterations}iter_'

        name += f'{get_valid_filename(self.heavy_v_germline)}_'
        name += f'{get_valid_filename(self.light_v_germline)}_'
        return name


@dataclass
class ManualHumanizationParams(HumanizationParams):
    method = 'manual'

    def get_export_name(self):
        return 'Designer_'


@dataclass
class HumanizedResidueAnnot:
    parental_score: float
    pred_score: float


@dataclass
class ChainHumanization:
    # Chain object containing the parental sequence
    parental_chain: Chain
    # Chain object containing the humanized sequence
    humanized_chain: Chain
    # Scores predicted from last parental sequence, can be used to explain the prediction
    scores: Dict[Position, Dict[str, float]]
    # Scores predicted from last humanized sequence, can be used to propose next mutations
    next_scores: Dict[Position, Dict[str, float]]

    @cached_property
    def alignment(self) -> Alignment:
        return self.parental_chain.align(self.humanized_chain)

    def num_mutations(self):
        return self.alignment.num_mutations()

    def get_alignment_string(self):
        chain_label = 'VH' if self.parental_chain.is_heavy_chain() else 'VL'
        return f'{self.parental_chain.name} {chain_label}\n{self.alignment}'

    def get_top_scores(self, n, next=False):
        scores = self.next_scores if next else self.scores
        top_scores = []
        for i in range(n):
            top_n_scores = []
            for pos, aa_scores in scores.items():
                aa, score = sorted(aa_scores.items(), key=lambda d: -d[1])[i]
                top_n_scores.append((pos, aa, score))
            top_scores.append(top_n_scores)
        return top_scores

    def to_score_dataframe(self, aligned_index=False) -> pd.DataFrame:
        scores = pd.DataFrame(self.scores.values())
        if aligned_index:
            scores.index = [p.format() for p in self.scores]
            scores.index.name = self.humanized_chain.scheme.title()
        else:
            scores.index.name = 'raw_pos'
        return scores


@dataclass
class AntibodyHumanization:
    vh: ChainHumanization
    vl: ChainHumanization

    def get_parental_chains(self) -> List[Chain]:
        chains = []
        if self.vh:
            chains.append(self.vh.parental_chain)
        if self.vl:
            chains.append(self.vl.parental_chain)
        return chains

    def get_humanized_chains(self) -> List[Chain]:
        chains = []
        if self.vh:
            chains.append(self.vh.humanized_chain)
        if self.vl:
            chains.append(self.vl.humanized_chain)
        return chains

    def get_alignment_string(self) -> str:
        strings = []
        if self.vh:
            strings.append(self.vh.get_alignment_string())
        if self.vl:
            strings.append(self.vl.get_alignment_string())
        return '\n'.join(strings)

    def to_score_dataframe(self) -> pd.DataFrame:
        scores = []
        if self.vh:
            scores.append(self.vh.to_score_dataframe())
        if self.vl:
            scores.append(self.vl.to_score_dataframe())
        return pd.concat(scores)


def humanize_antibody(vh: Optional[Chain], vl: Optional[Chain], params: HumanizationParams) -> AntibodyHumanization:
    return AntibodyHumanization(
        vh=humanize_chain(vh, params=params) if vh else None,
        vl=humanize_chain(vl, params=params) if vl else None
    )


def humanize_chain(parental_chain: Chain, params: HumanizationParams) -> ChainHumanization:
    assert parental_chain.cdr_definition == params.cdr_definition, \
        f'Expected chain with {params.cdr_definition} CDR definition, got {parental_chain.cdr_definition}'

    if isinstance(params, SapiensHumanizationParams):
        return sapiens_humanize_chain(parental_chain, params=params)
    if isinstance(params, CDRGraftingHumanizationParams):
        return cdr_grafting_humanize_chain(parental_chain, params=params)
    if isinstance(params, ManualHumanizationParams):
        return sapiens_humanize_chain(parental_chain, params=SapiensHumanizationParams(iterations=0))
    raise NotImplementedError(f'Unknown humanization method params {type(params)}')


def cdr_grafting_humanize_chain(parental_chain: Chain, params: CDRGraftingHumanizationParams) -> ChainHumanization:
    v_gene = params.heavy_v_germline if parental_chain.is_heavy_chain() else params.light_v_germline
    if v_gene == 'auto':
        v_gene = None

    humanized_chain = parental_chain.graft_cdrs_onto_human_germline(
        v_gene=v_gene,
        backmutate_vernier=params.backmutate_vernier
    )

    # Compute Sapiens scores (used for Designer and final Sapiens pass if enabled)
    sapiens_humanization = sapiens_humanize_chain(
        humanized_chain,
        params=SapiensHumanizationParams(
            iterations=params.sapiens_iterations,
            backmutate_vernier=params.backmutate_vernier
        )
    )
    if params.sapiens_iterations:
        humanized_chain = sapiens_humanization.humanized_chain

    return ChainHumanization(
        parental_chain=parental_chain,
        humanized_chain=humanized_chain,
        scores=sapiens_humanization.scores,
        next_scores=sapiens_humanization.next_scores
    )


def sapiens_humanize_chain(parental_chain: Chain, params: SapiensHumanizationParams) -> ChainHumanization:
    # Repeat Sapiens multiple times if requested, we start with the parental chain
    humanized_chain = parental_chain.clone()

    pred = None
    for iteration in range(params.iterations):
        # Get Sapiens scores as a positions (rows) by amino acids (columns) matrix
        pred = sapiens_predict_chain(humanized_chain, model_version=params.model_version)

        # Create humanized sequence by taking the amino acid with highest score at each position
        humanized_seq = ''.join(pred.idxmax(axis=1).values)

        # Reuse numbering from parental chain to make sure that we get the same numbering
        # Note: Chain can be aligned from scratch, but downstream analysis needs to cope with different numbering
        # pred_chain = Chain(pred_seq, scheme=chain.scheme, cdr_definition=chain.cdr_definition, name=chain.name)
        humanized_chain = parental_chain.clone(humanized_seq)

        # Graft parental CDRs into the humanized sequence, unless humanizing CDRs as well
        if not params.humanize_cdrs:
            humanized_chain = parental_chain.graft_cdrs_onto(
                humanized_chain,
                backmutate_vernier=params.backmutate_vernier
            )
        else:
            if params.backmutate_vernier:
                raise ValueError('Cannot backmutate Vernier regions when humanizing CDRs')

    if pred is None:
        # Support case with 0 iterations, still return probabilities
        pred = sapiens_predict_chain(parental_chain, model_version=params.model_version)

    # Predict scores of potential next mutation from the final humanized sequence
    pred_next = sapiens_predict_chain(humanized_chain, model_version=params.model_version)

    return ChainHumanization(
        parental_chain=parental_chain,
        humanized_chain=humanized_chain,
        scores={pos: row.to_dict() for pos, (i, row) in zip(humanized_chain.positions, pred.iterrows())},
        next_scores={pos: row.to_dict() for pos, (i, row) in zip(humanized_chain.positions, pred_next.iterrows())},
    )


def sapiens_predict_chain(chain, model_version='latest', return_all_hiddens=False):
    return sapiens.predict_scores(
        seq=chain.seq,
        chain_type=chain.chain_type,
        model_version=model_version,
        return_all_hiddens=return_all_hiddens
    )


def generate_samples(input_df: pd.DataFrame, initial_n_samples: int) -> List[str]:
    input_ndarray = input_df.to_numpy(dtype="float32")
    n_rows = input_ndarray.shape[0]
    n_possible_aa_residues = input_ndarray.shape[1]
    aa_array=input_df.columns.values
    rng = default_rng(seed=42)
    
    # Ensure that we the probabilities for each residue sum to 1
    fixed_probabilities = input_ndarray / input_ndarray.sum(axis=1,keepdims=1)
    # make an ndarray of shape (n_samples, n_rows), where each row corresponds to a possible sequence,
    # and each value in the row corresponds to the aa to be chosen from the aa_array
    # If sampled[0] = [1,4,5], then that sequence would correspond to CFG
    # this could be simplified with something like
    # https://stackoverflow.com/questions/34187130/fast-random-weighted-selection-across-all-rows-of-a-stochastic-matrix?noredirect=1&lq=1
    # https://stackoverflow.com/questions/47722005/vectorizing-numpy-random-choice-for-given-2d-array-of-probabilities-along-an-a
    # https://stackoverflow.com/questions/40474436/how-to-apply-numpy-random-choice-to-a-matrix-of-probability-values-vectorized-s
    list_of_random_samples = [
        rng.choice(
            n_possible_aa_residues,
            initial_n_samples,
            p=fixed_probabilities[i]
        ) for i in range(n_rows)
    ]
    sampled = np.stack(list_of_random_samples, axis=1)
    # Get unique_sequences, many of the sequences will be duplicates given that each
    # residue normally has one aa with a very high probability
    unique_rows = np.unique(sampled, axis=0)
    # Convert indices of aa_array to strings
    sequences = ["".join(np.take(aa_array, unique_rows[i])) for i in range(unique_rows.shape[0])]
    return sequences

def generate_seq_records_from_sequence_strings(
    sequences: List[str],
    id_prefix: str,
    description: str,
) -> List[SeqRecord]:
        return [
            SeqRecord(
                Seq(seq),
                id=id_prefix,
                description=(
                    f'{description} {i}'
                )
            )
            for i, seq in enumerate(sequences)
        ]
