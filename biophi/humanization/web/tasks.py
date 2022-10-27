from copy import deepcopy
from dataclasses import dataclass
from typing import List, Dict

from Bio.SeqRecord import SeqRecord

from biophi.common.utils.io import AntibodyInput
from biophi.common.web.tasks import celery
from biophi.humanization.methods.humanization import HumanizationParams, \
    AntibodyHumanization, humanize_antibody
from biophi.humanization.methods.humanness import OASisParams, \
    AntibodyHumanness, get_antibody_humanness, ChainHumanness, get_chain_humanness
from abnumber import Chain, Position
import pandas as pd


class HumanizeAntibodyTaskError(Exception):
    def __init__(self, message: str = None, input: AntibodyInput = None):
        super().__init__(message)
        self.input: AntibodyInput = input


@dataclass
class HumanizeAntibodyTaskResult:
    humanization: AntibodyHumanization
    parental_humanness: AntibodyHumanness
    humanized_humanness: AntibodyHumanness
    input: AntibodyInput
    humanization_params: HumanizationParams
    oasis_params: OASisParams

    def to_series(self):
        threshold = self.oasis_params.min_fraction_subjects_label
        return pd.Series({
            ('Info', 'Antibody'): self.input.name,
            ('OASis Percentile', 'Before'): self.parental_humanness.get_oasis_percentile(self.oasis_params.min_fraction_subjects),
            ('OASis Percentile', 'After'): self.humanized_humanness.get_oasis_percentile(self.oasis_params.min_fraction_subjects),
            ('OASis Identity', 'Before'): self.parental_humanness.get_oasis_identity(self.oasis_params.min_fraction_subjects),
            ('OASis Identity', 'After'): self.humanized_humanness.get_oasis_identity(self.oasis_params.min_fraction_subjects),
            ('Non-human Peptides', 'Before'): self.parental_humanness.get_num_nonhuman_peptides(self.oasis_params.min_fraction_subjects),
            ('Non-human Peptides', 'After'): self.humanized_humanness.get_num_nonhuman_peptides(self.oasis_params.min_fraction_subjects),
            ('Settings', 'Threshold'): threshold
        })

    @classmethod
    def to_overview_dataframe(cls, results: List['HumanizeAntibodyTaskResult']):
        return pd.DataFrame([r.to_series() for r in results])

    @classmethod
    def to_oasis_curve_dataframe(cls, results: List['HumanizeAntibodyTaskResult']):
        types = []
        index = []
        humanness_list = []
        for result in results:
            index.append(result.input.name)
            types.append('Parental')
            humanness_list.append(result.parental_humanness)
            index.append(result.input.name)
            types.append('Humanized')
            humanness_list.append(result.humanized_humanness)
        df = pd.DataFrame(
            [humanness.get_oasis_curve() for humanness in humanness_list]
        )
        df['Antibody'] = index
        df['Type'] = types
        df.set_index(['Antibody', 'Type'], inplace=True)
        return df

    @classmethod
    def to_sheets(cls, results: List['HumanizeAntibodyTaskResult'], full=True) -> Dict[str, pd.DataFrame]:
        sheets = {
            'Overview': cls.to_overview_dataframe(results)
        }
        if full:
            sheets['VH'] = ChainHumanness.to_sequence_dataframe(
                [r.humanized_humanness.vh for r in results if r.humanized_humanness.vh],
                species=False
            )
            sheets['VL'] = ChainHumanness.to_sequence_dataframe(
                [r.humanized_humanness.vl for r in results if r.humanized_humanness.vl],
                species=False
            )
            sheets['OASis Curves'] = cls.to_oasis_curve_dataframe(results)
            for result in results:
                peptides = result.humanized_humanness.to_peptide_dataframe()
                scores = result.humanization.to_score_dataframe()
                assert len(peptides) == len(scores), 'Scores and peptides not aligned in ' + result.input.name
                sheets[result.input.name] = pd.concat([peptides, scores], axis=1)
        return sheets

    def get_export_name(self, num_seqs=1):
        name = 'Humanized_'
        if num_seqs > 1:
            name += f'{num_seqs}seqs_'
        else:
            name += self.input.safe_name + '_'
        name += self.humanization_params.get_export_name()
        name += 'BioPhi'
        return name

    def get_humanized_records(self) -> List[SeqRecord]:
        records = []
        if self.humanization.vh:
            records.append(self.humanization.vh.humanized_chain.to_seq_record(
                description=f'VH ({self.get_export_name().replace("_"," ")})'
            ))
        if self.humanization.vl:
            records.append(self.humanization.vl.humanized_chain.to_seq_record(
                description=f'VL ({self.get_export_name().replace("_"," ")})'
            ))
        return records

    def get_parental_records(self) -> List[SeqRecord]:
        records = []
        if self.humanization.vh:
            records.append(self.humanization.vh.parental_chain.to_seq_record(description='VH'))
        if self.humanization.vl:
            records.append(self.humanization.vl.humanized_chain.to_seq_record(description='VL'))
        return records


@celery.task()
def humanize_antibody_task(input: AntibodyInput, humanization_params: HumanizationParams, oasis_params: OASisParams) -> HumanizeAntibodyTaskResult:
    vh_chain, vl_chain = None, None
    try:
        if input.heavy_protein_seq:
            vh_chain = Chain(
                input.heavy_protein_seq,
                scheme=humanization_params.scheme,
                cdr_definition=humanization_params.cdr_definition,
                name=input.name
            )

        if input.light_protein_seq:
            vl_chain = Chain(
                input.light_protein_seq,
                scheme=humanization_params.scheme,
                cdr_definition=humanization_params.cdr_definition,
                name=input.name
            )

        if vh_chain is not None and not vh_chain.is_heavy_chain():
            raise ValueError(f'Expected heavy chain, got {vh_chain.chain_type} chain: {vh_chain.seq}')

        if vl_chain is not None and not vl_chain.is_light_chain():
            raise ValueError(f'Expected light chain, got {vl_chain.chain_type} chain: {vl_chain.seq}')

        humanization = humanize_antibody(
            vh=vh_chain,
            vl=vl_chain,
            params=humanization_params
        )

        return HumanizeAntibodyTaskResult(
            humanization=humanization,
            parental_humanness=get_antibody_humanness(
                vh=vh_chain,
                vl=vl_chain,
                params=oasis_params
            ),
            humanized_humanness=get_antibody_humanness(
                vh=humanization.vh.humanized_chain if vh_chain else None,
                vl=humanization.vl.humanized_chain if vl_chain else None,
                params=oasis_params
            ),
            input=input,
            humanization_params=humanization_params,
            oasis_params=oasis_params
        )
    except Exception as e:
        raise HumanizeAntibodyTaskError(str(e), input=input) from e


@celery.task()
def mutate_humanized_antibody_task(result: HumanizeAntibodyTaskResult, pos: str, aa: str):
    assert len(aa) == 1, f'Invalid amino acid: {aa}'

    result = deepcopy(result)

    if pos.startswith('H'):
        chain = result.humanization.vh.humanized_chain
    elif pos.startswith('L'):
        chain = result.humanization.vl.humanized_chain
    else:
        raise ValueError(f'Position should start with H or L, got: "{pos}"')

    pos_object = chain._parse_position(pos)
    found = False
    for name, region in chain.regions.items():
        if pos_object in region:
            region[pos_object] = aa
            found = True
            break
    assert found, f'Position "{pos_object}" not in chain: \n{chain}'

    humanness = get_chain_humanness(chain, params=result.oasis_params)

    if pos.startswith('H'):
        result.humanized_humanness.vh = humanness
    elif pos.startswith('L'):
        result.humanized_humanness.vl = humanness
    else:
        raise ValueError()

    return result


class HumannessTaskError(Exception):
    def __init__(self, message: str = None, input: AntibodyInput = None):
        super().__init__(message)
        self.input: AntibodyInput = input


@dataclass
class HumannessTaskResult:
    input: AntibodyInput
    oasis_params: OASisParams
    humanness: AntibodyHumanness

    def to_series(self):
        threshold = self.oasis_params.min_fraction_subjects_label
        series = {
            'Antibody': self.input.name,
            'Threshold': threshold,
            'OASis Percentile': self.humanness.get_oasis_percentile(self.oasis_params.min_fraction_subjects),
            'OASis Identity': self.humanness.get_oasis_identity(self.oasis_params.min_fraction_subjects),
            'Germline Content': self.humanness.get_germline_content()
        }
        if self.humanness.vh:
            series['Heavy V Germline'] = self.humanness.vh.v_germline_names[0]
            series['Heavy J Germline'] = self.humanness.vh.j_germline_names[0]
            series['Heavy OASis Percentile'] = self.humanness.vh.get_oasis_percentile(self.oasis_params.min_fraction_subjects)
            series['Heavy OASis Identity'] = self.humanness.vh.get_oasis_identity(self.oasis_params.min_fraction_subjects)
            series['Heavy Non-human peptides'] = self.humanness.vh.get_num_nonhuman_peptides(self.oasis_params.min_fraction_subjects),
            series['Heavy Germline Content'] = self.humanness.vh.get_germline_content()
        if self.humanness.vl:
            series['Light V Germline'] = self.humanness.vl.v_germline_names[0]
            series['Light J Germline'] = self.humanness.vl.j_germline_names[0]
            series['Light OASis Percentile'] = self.humanness.vl.get_oasis_percentile(self.oasis_params.min_fraction_subjects)
            series['Light OASis Identity'] = self.humanness.vl.get_oasis_identity(self.oasis_params.min_fraction_subjects)
            series['Light Non-human peptides'] = self.humanness.vl.get_num_nonhuman_peptides(self.oasis_params.min_fraction_subjects)
            series['Light Germline Content'] = self.humanness.vl.get_germline_content()

        return pd.Series(series)

    @classmethod
    def to_overview_dataframe(cls, results: List['HumannessTaskResult']):
        return pd.DataFrame([r.to_series() for r in results]).set_index('Antibody')

    @classmethod
    def to_oasis_curve_dataframe(cls, results: List['HumannessTaskResult']) -> pd.DataFrame:
        df = pd.DataFrame(
            [result.humanness.get_oasis_curve() for result in results],
            index=[result.input.name for result in results]
        )
        df.index.name = 'Antibody'
        return df

    @classmethod
    def to_sheets(cls, results: List['HumannessTaskResult'], full=True):
        sheets = {
            'Overview': cls.to_overview_dataframe(results)
        }
        if full:
            sheets['VH'] = ChainHumanness.to_sequence_dataframe(
                [r.humanness.vh for r in results if r.humanness.vh],
                species=False
            )
            sheets['VL'] = ChainHumanness.to_sequence_dataframe(
                [r.humanness.vl for r in results if r.humanness.vl],
                species=False
            )
            sheets['OASis Curves'] = cls.to_oasis_curve_dataframe(results)
            for result in results:
                sheets[result.input.name] = result.humanness.to_peptide_dataframe()
        return sheets

    def get_records(self) -> List[SeqRecord]:
        records = []
        if self.humanness.vh:
            records.append(self.humanness.vh.chain.to_seq_record(description='VH'))
        if self.humanness.vl:
            records.append(self.humanness.vl.chain.to_seq_record(description='VL'))
        return records

@celery.task()
def humanness_task(input: AntibodyInput, oasis_params: OASisParams, scheme='kabat', cdr_definition='kabat'):
    vh_chain, vl_chain = None, None
    try:
        if input.heavy_protein_seq:
            vh_chain = Chain(input.heavy_protein_seq, scheme=scheme, cdr_definition=cdr_definition, name=input.name)

        if input.light_protein_seq:
            vl_chain = Chain(input.light_protein_seq, scheme=scheme, cdr_definition=cdr_definition, name=input.name)

        if vh_chain is not None and not vh_chain.is_heavy_chain():
            raise ValueError(f'Expected heavy chain, got {vh_chain.chain_type} chain: {vh_chain.seq}')

        if vl_chain is not None and not vl_chain.is_light_chain():
            raise ValueError(f'Expected light chain, got {vl_chain.chain_type} chain: {vl_chain.seq}')

        return HumannessTaskResult(
            input=input,
            oasis_params=oasis_params,
            humanness=get_antibody_humanness(
                vh=vh_chain,
                vl=vl_chain,
                params=oasis_params
            )
        )
    except Exception as e:
        raise HumannessTaskError(str(e), input=input) from e
