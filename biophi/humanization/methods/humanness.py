import os
from dataclasses import dataclass
from functools import cached_property, lru_cache
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from abnumber.chain import Chain, Position
import pandas as pd
from abnumber.germlines import get_imgt_v_chains, get_imgt_j_chains
from sqlalchemy import create_engine
from typing import Dict, Union, List, Optional, Tuple
import numpy as np
from sqlalchemy.engine import Engine
from biophi.humanization.methods.stats import get_oasis_percentile, \
    get_germline_family_residue_frequency, get_chain_type_residue_frequency

OASIS_MIN_SUBJECTS_THRESHOLDS = {
    'loose': 0.01,
    'relaxed': 0.1,
    'medium': 0.5,
    'strict': 0.9
}

@dataclass
class OASisParams:
    oasis_db_path: str
    min_fraction_subjects: float

    @property
    def min_fraction_subjects_label(self):
        label = f'{self.min_fraction_subjects:.0%}'
        for label, threshold in OASIS_MIN_SUBJECTS_THRESHOLDS.items():
            if threshold == self.min_fraction_subjects:
                break
        return label


@dataclass
class PeptideHumanness:
    seq: str
    num_oas_subjects: int
    fraction_oas_subjects: float
    num_oas_occurrences: int

    def is_human(self, min_fraction_subjects):
        assert 0 < min_fraction_subjects < 1, 'min_fraction_subjects should be between 0 and 1'
        if self.fraction_oas_subjects is None:
            return np.nan
        return self.fraction_oas_subjects >= min_fraction_subjects

    def to_series(self):
        return pd.Series({
            'Seq': self.seq,
            'Num OAS Subjects': self.num_oas_subjects,
            'Fraction OAS Subjects': self.fraction_oas_subjects,
            'Num OAS Hits': self.num_oas_occurrences
        })


@dataclass
class ChainHumanness:
    chain: Chain
    imgt_chain: Chain
    peptides: Dict[Position, PeptideHumanness]
    num_germline_residues: int
    v_germline_names: List[str]
    j_germline_names: List[str]
    v_germline_family: str
    v_germline_suffix: str
    germline_family_residue_frequency: Dict[Position, Dict[str, float]]
    chain_type_residue_frequency: Dict[Position, Dict[str, float]]

    @cached_property
    def pos_imgt_mapping(self) -> Dict[Position, Position]:
        return {pos: imgt_pos for pos, imgt_pos in zip(self.chain.positions, self.imgt_chain.positions)}

    def get_oasis_curve(self, frequency=True, cumulative=True) -> pd.Series:
        bins = [get_fraction_subjects_bin(peptide.fraction_oas_subjects) for peptide in self.peptides.values()]
        curve = pd.Series(bins, name=self.chain.name).value_counts().reindex(FRACTION_SUBJECTS_BINS).fillna(0)

        if frequency:
            curve = curve / curve.sum()
        else:
            curve = curve.astype(np.int)

        if cumulative:
            curve = curve[::-1].cumsum()[::-1]
        return curve[FRACTION_SUBJECTS_BINS]

    def get_oasis_identity(self, min_fraction_subjects) -> float:
        return self.get_num_human_peptides(min_fraction_subjects) / self.get_num_peptides()

    def get_oasis_percentile(self, min_fraction_subjects) -> float:
        return get_oasis_percentile(
            chain_type=self.chain.chain_type,
            oasis_identity=self.get_oasis_identity(min_fraction_subjects),
            min_fraction_subjects=min_fraction_subjects
        )

    def get_v_germline_chains(self, limit=None) -> List[Chain]:
        germlines = get_imgt_v_chains(self.chain.chain_type)
        return [germlines[name] for name in self.v_germline_names[:limit]]

    def get_j_germline_chains(self, limit=None) -> List[Chain]:
        germlines = get_imgt_j_chains(self.chain.chain_type)
        return [germlines[name] for name in self.j_germline_names[:limit]]

    def get_num_human_peptides(self, min_fraction_subjects) -> int:
        return sum(p.is_human(min_fraction_subjects) for p in self.peptides.values())

    def get_num_nonhuman_peptides(self, max_fraction_subjects) -> int:
        return sum(not p.is_human(max_fraction_subjects) for p in self.peptides.values())

    def get_num_peptides(self) -> int:
        return len(self.peptides)

    def has_position(self, position) -> bool:
        return position in self.peptides

    def get_peptide(self, position, edges=True) -> PeptideHumanness:
        if self.has_position(position):
            return self.peptides[position]
        if not edges:
            raise KeyError(str(position))
        positions = list(self.peptides.keys())
        return self.peptides[positions[0] if position < positions[0] else positions[-1]]

    def get_peptide_length(self) -> int:
        lengths = [len(peptide.seq) for peptide in self.peptides.values()]
        assert len(set(lengths)) == 1, f'Peptides should have same lengths, got: {set(lengths)}'
        return lengths[0]

    def get_positional_humanness(self, min_fraction_subjects) -> List[Tuple[Position, str, List[PeptideHumanness]]]:
        chain_positions = list(self.chain.positions)
        peptide_len = self.get_peptide_length()
        assert peptide_len % 2 == 1, 'Peptide needs to have odd length'
        annots = []
        for raw_pos, (pos, aa) in enumerate(self.chain):
            window = [chain_positions[i] for i in range(max(0, raw_pos-peptide_len+1), raw_pos+1)]
            peptides = [self.peptides[peptide_pos]
                         for peptide_pos in window if peptide_pos in self.peptides]
            non_human_peptides = [p for p in peptides if not p.is_human(min_fraction_subjects)]
            annots.append((pos, aa, non_human_peptides))
        return annots

    def to_peptide_dataframe(self) -> pd.DataFrame:
        scheme = list(self.peptides)[0].scheme
        imgt_positions = self.imgt_chain.positions.keys()
        germline = {
            **self.get_j_germline_chains(1)[0].positions,
            **self.get_v_germline_chains(1)[0].positions
        }
        scheme_col = f'{scheme.title()}'
        df = pd.DataFrame([{
            'Antibody': self.chain.name,
            'Chain type': self.chain.chain_type,
            'Region': pos.get_region(),
            scheme_col: pos,
            'AA': aa,
            'Germline': germline.get(imgt_pos, '-'),
            'Chain Freq': self.chain_type_residue_frequency[pos].get(aa, 0)
                                if self.chain_type_residue_frequency.get(pos) else None,
            'Family Freq': self.germline_family_residue_frequency[pos].get(aa, 0)
                                if self.germline_family_residue_frequency.get(pos) else None
        } for (pos, aa), imgt_pos in zip(self.chain, imgt_positions)])
        df = df.merge(pd.DataFrame(
            [pep.to_series() for pep in self.peptides.values()],
            index=list(self.peptides.keys())
        ).add_prefix('Peptide '), left_on=scheme_col, right_index=True, how='left')
        return df

    @classmethod
    def to_sequence_dataframe(cls, humanness_list: List['ChainHumanness'], species=True) -> pd.DataFrame:
        chains = [humanness.chain for humanness in humanness_list]
        df = Chain.to_dataframe(chains)
        df.insert(0, 'scheme', [chain.scheme for chain in chains])
        if not species and 'species' in df.columns:
            df.drop(columns=['species'], inplace=True)
        df.index.name = 'Antibody'
        return df

    def get_germline_content(self):
        return self.num_germline_residues / len(self.chain)

    def get_top_freqs(self, n):
        top_freqs = []
        for i in range(n):
            top_n_freqs = []
            for pos, aa_freqs in self.germline_family_residue_frequency.items():
                if aa_freqs is not None and len(aa_freqs) > i:
                    aa, freq = sorted(aa_freqs.items(), key=lambda d: -d[1])[i]
                else:
                    aa, freq = 'X', 0
                top_n_freqs.append((pos, aa, freq))
            top_freqs.append(top_n_freqs)
        return top_freqs


@dataclass
class AntibodyHumanness:
    vh: ChainHumanness
    vl: ChainHumanness

    def get_oasis_identity(self, min_fraction_subjects) -> float:
        return self.get_num_human_peptides(min_fraction_subjects) / self.get_num_peptides()

    def get_oasis_percentile(self, min_fraction_subjects) -> float:
        return get_oasis_percentile(
            chain_type='mean',
            oasis_identity=self.get_oasis_identity(min_fraction_subjects),
            min_fraction_subjects=min_fraction_subjects
        )

    def get_oasis_curve(self, frequency=True) -> pd.Series:
        if self.vh and self.vl:
            counts_curve = self.vh.get_oasis_curve(frequency=False, cumulative=False) \
                           + self.vl.get_oasis_curve(frequency=False, cumulative=False)
            if frequency:
                counts_curve = counts_curve / counts_curve.sum()
            return counts_curve[::-1].cumsum()[::-1]
        if self.vh:
            return self.vh.get_oasis_curve(frequency=frequency)
        if self.vl:
            return self.vl.get_oasis_curve(frequency=frequency)

    def get_num_peptides(self) -> int:
        return (self.vh.get_num_peptides() if self.vh else 0) + (self.vl.get_num_peptides() if self.vl else 0)

    def get_num_human_peptides(self, min_fraction_subjects) -> int:
        return (self.vh.get_num_human_peptides(min_fraction_subjects) if self.vh else 0) \
               + (self.vl.get_num_human_peptides(min_fraction_subjects) if self.vl else 0)

    def get_num_nonhuman_peptides(self, min_fraction_subjects) -> int:
        return (self.vh.get_num_nonhuman_peptides(min_fraction_subjects) if self.vh else 0) \
               + (self.vl.get_num_nonhuman_peptides(min_fraction_subjects) if self.vl else 0)

    def to_peptide_dataframe(self) -> pd.DataFrame:
        dfs = []
        if self.vh:
            dfs.append(self.vh.to_peptide_dataframe())
        if self.vl:
            dfs.append(self.vl.to_peptide_dataframe())
        return pd.concat(dfs)

    def get_germline_content(self):
        num_germline_residues = 0
        num_total_residues = 0
        if self.vh:
            num_germline_residues += self.vh.num_germline_residues
            num_total_residues += len(self.vh.chain)
        if self.vl:
            num_germline_residues += self.vl.num_germline_residues
            num_total_residues += len(self.vl.chain)
        if num_total_residues == 0:
            return None
        return num_germline_residues / num_total_residues


def chop_seq_peptides(seq: Union[SeqRecord, Chain], peptide_length):
    if isinstance(seq, SeqRecord):
        seq = str(seq.seq)
        # numbered list from 1 to length of sequence
        positions = list(range(1, len(seq) + 1))
    elif isinstance(seq, Chain):
        # numbers are coming from Chain Position objects
        positions = list(seq.positions.keys())
        seq = ''.join(seq.positions.values())
    else:
        raise ValueError(f'Unsupported sequence type: {type(seq)}')

    return [(pos, seq[i:i + peptide_length]) for i, pos in enumerate(positions[:-peptide_length+1])]


def get_antibody_humanness(vh: Optional[Chain], vl: Optional[Chain], params: OASisParams) -> AntibodyHumanness:
    return AntibodyHumanness(
        vh=get_chain_humanness(vh, params=params) if vh else None,
        vl=get_chain_humanness(vl, params=params) if vl else None
    )


def get_chain_oasis_peptides(chain, params: OASisParams):
    pos_peptides = chop_seq_peptides(chain, peptide_length=9)

    if params.oasis_db_path is None:
        return {pos: PeptideHumanness(
                    seq=peptide,
                    num_oas_subjects=None,
                    fraction_oas_subjects=None,
                    num_oas_occurrences=None
                ) for pos, peptide in pos_peptides}

    if params.oasis_db_path.endswith('.gz'):
        raise ValueError('The OASis DB file needs to be unzipped (use "gunzip DB_PATH.db.gz")')

    if not os.path.exists(params.oasis_db_path):
        raise FileNotFoundError(f'The OASis DB path does not exist: {params.oasis_db_path}')

    oas_engine = create_engine('sqlite:///' + os.path.abspath(params.oasis_db_path), echo=False)

    oas_filter_chain = "Heavy" if chain.is_heavy_chain() else "Light"
    result = oas_engine.execute(f'SELECT COUNT(*) FROM subjects WHERE Complete{oas_filter_chain}Seqs >= 10000')
    num_total_oas_subjects = int(result.fetchall()[0][0])

    oas_hits = get_oas_hits(
        [peptide for pos, peptide in pos_peptides],
        engine=oas_engine,
        filter_chain=oas_filter_chain
    )
    oas_grouped = oas_hits.groupby('peptide')

    return {
        pos: parse_peptide_humanness(
            pep,
            oas_hits=oas_grouped.get_group(pep) if pep in oas_grouped.groups else oas_hits[:0],
            num_total_oas_subjects=num_total_oas_subjects
        ) for pos, pep in pos_peptides
    }


def get_chain_humanness(chain: Chain, params: OASisParams) -> ChainHumanness:

    peptides = get_chain_oasis_peptides(chain, params=params)

    imgt_chain = chain if chain.scheme == 'imgt' and chain.cdr_definition == 'imgt' else chain.renumber('imgt')

    v_germline_chains, j_germline_chains = chain.find_human_germlines(limit=10)
    top_v, top_j = v_germline_chains[0], j_germline_chains[0]
    num_germline_residues = sum(top_v.positions.get(pos) == aa or top_j.positions.get(pos) == aa
                                    for pos, aa in imgt_chain)

    v_germline_family = top_v.name.split('-')[0].split('/')[0]
    v_germline_suffix = top_v.name.replace(v_germline_family, '')
    return ChainHumanness(
        chain=chain,
        imgt_chain=imgt_chain,
        num_germline_residues=num_germline_residues,
        v_germline_names=[chain.name for chain in v_germline_chains],
        j_germline_names=[chain.name for chain in j_germline_chains],
        peptides=peptides,
        v_germline_family=v_germline_family,
        v_germline_suffix=v_germline_suffix,
        germline_family_residue_frequency=get_germline_family_residue_frequency(chain, imgt_chain, v_germline_family),
        chain_type_residue_frequency=get_chain_type_residue_frequency(chain, imgt_chain)
    )


def parse_peptide_humanness(peptide: str, oas_hits: pd.DataFrame, num_total_oas_subjects: int = None) -> PeptideHumanness:
    return PeptideHumanness(
        seq=peptide,
        num_oas_subjects=len(oas_hits),
        fraction_oas_subjects=len(oas_hits) / num_total_oas_subjects if num_total_oas_subjects else None,
        num_oas_occurrences=oas_hits['count'].sum()
    )


def get_oas_hits(peptides: Union[str, List[str]], engine: Engine, filter_chain=None):
    if not isinstance(peptides, list):
        peptides = list(peptides)
    filter_chain_statement = ""
    if filter_chain:
        assert filter_chain in ['Heavy', 'Light']
        filter_chain_statement = f"AND Complete{filter_chain}Seqs >= 10000"

    statement = "SELECT peptides.* FROM peptides " \
                "LEFT JOIN subjects ON peptides.subject=subjects.id " \
                "WHERE peptide IN (" + ",".join("?" * len(peptides)) + ") AND subjects.StudyPath <> 'Corcoran_2016' " \
                + filter_chain_statement
    return pd.read_sql(statement, params=peptides, con=engine)


FRACTION_SUBJECTS_THRESHOLDS = np.arange(0, 0.91, 0.01)
FRACTION_SUBJECTS_BINS = [f'{v:.0%}' for v in FRACTION_SUBJECTS_THRESHOLDS]


def get_fraction_subjects_bin(fraction):
    for bin_threshold, bin_label in zip(FRACTION_SUBJECTS_THRESHOLDS[::-1], FRACTION_SUBJECTS_BINS[::-1]):
        if bin_threshold is None or fraction >= bin_threshold:
            return bin_label
