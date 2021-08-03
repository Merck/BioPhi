import sys
import time
from io import StringIO, BytesIO
from typing import Dict, Optional, List, Union
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.Chain import Chain as BioPdbChain
from Bio.PDB.Model import Model
from Bio.SeqUtils import seq1
from abnumber import Chain as AbNumberChain, ChainParseError
from flask import send_file, request, flash, current_app
from dataclasses import dataclass
from werkzeug.datastructures import FileStorage
import pandas as pd
from werkzeug.utils import redirect

from biophi.common.utils.formatting import get_valid_filename
from biophi.common.utils.seq import looks_like_antibody_heavy_chain, download_pdb, parse_plaintext_records, \
    looks_like_dna, is_valid_amino_acid_sequence
import os
from biophi.common.utils.seq import sanitize_sequence
from biophi.common.utils.stats import log_submission


@dataclass
class AntibodyPDB:
    pdb_data: str
    heavy_chain_name: str
    light_chain_name: str

    def get_model(self) -> Model:
        structure = PDBParser().get_structure('structure', StringIO(self.pdb_data))
        return structure[0]

    def get_heavy_pdb_chain(self) -> BioPdbChain:
        model = self.get_model()
        return model[self.heavy_chain_name]

    def get_light_pdb_chain(self) -> BioPdbChain:
        model = self.get_model()
        return model[self.light_chain_name]

    def get_heavy_sequence(self):
        return ''

    def get_light_sequence(self):
        return ''

    def convert_heavy_positions(self, ab_chain: AbNumberChain) -> Dict[str, str]:
        return convert_pdb_positions(self.get_heavy_pdb_chain(), ab_chain)

    def convert_light_positions(self, ab_chain: AbNumberChain) -> Dict[str, str]:
        return convert_pdb_positions(self.get_light_pdb_chain(), ab_chain)


@dataclass
class AntibodyInput:
    name: str
    heavy_protein_seq: str
    light_protein_seq: str
    pdb: AntibodyPDB = None

    def __post_init__(self):
        if self.heavy_protein_seq is not None:
            assert is_valid_amino_acid_sequence(self.heavy_protein_seq), \
                f'Invalid sequence for "{self.name}": "{self.heavy_protein_seq}"'
        if self.light_protein_seq is not None:
            assert is_valid_amino_acid_sequence(self.light_protein_seq), \
                f'Invalid sequence for "{self.name}": "{self.light_protein_seq}"'

    @property
    def safe_name(self):
        return get_valid_filename(self.name)

    @classmethod
    def from_pdb_id(cls, pdb_id):
        pdb_data = download_pdb(pdb_id).decode()
        return cls.from_pdb_data(pdb_data, name=pdb_id)

    @classmethod
    def from_pdb_data(cls, pdb_data, name) -> Optional['AntibodyInput']:
        heavy_record, light_record = None, None
        for record in SeqIO.parse(StringIO(pdb_data), 'pdb-seqres'):
            try:
                if looks_like_antibody_heavy_chain(record.seq):
                    if heavy_record is None:
                        heavy_record = record
                else:
                    if light_record is None:
                        light_record = record
            except ChainParseError:
                pass
            if heavy_record is not None and light_record is not None:
                break

        if heavy_record is None and light_record is None:
            return None

        return AntibodyInput(
            name=name,
            heavy_protein_seq=str(heavy_record.seq) if heavy_record else None,
            light_protein_seq=str(light_record.seq) if light_record else None,
            pdb=AntibodyPDB(
                pdb_data=pdb_data,
                heavy_chain_name=heavy_record.annotations['chain'] if heavy_record else None,
                light_chain_name=light_record.annotations['chain'] if light_record else None
            )
        )


def read_antibody_input_request():
    input_mode = request.form.get('input_mode', 'single')
    if input_mode == 'single':
        vh = sanitize_sequence(request.form['vh'])
        vl = sanitize_sequence(request.form['vl'])

        if not vh and not vl:
            flash('No input sequences were provided')
            return redirect(request.url)

        input = AntibodyInput(
            name=request.form['name'] or 'Antibody1',
            heavy_protein_seq=vh,
            light_protein_seq=vl
        )

        log_submission(
            antibody_inputs=[input],
            invalid_names=[],
            duplicate_names=[],
            unrecognized_files=[]
        )
        return [input]

    antibody_inputs, invalid_names, duplicate_names, unrecognized_files = parse_antibody_inputs(
        seq_string=request.form['sequence_text'],
        pdb_ids=[pdb_id.strip() for pdb_id in request.form['pdb_ids'].replace(',', ' ').split() if pdb_id.strip()],
        files=request.files.getlist("sequence_files[]")
    )

    log_submission(
        antibody_inputs=antibody_inputs,
        invalid_names=invalid_names,
        duplicate_names=duplicate_names,
        unrecognized_files=unrecognized_files
    )

    if invalid_names:
        flash(f'Skipped {len(invalid_names):,} invalid inputs: {", ".join(invalid_names[:3])}')

    if duplicate_names:
        flash(f'Skipped {len(duplicate_names):,} inputs with duplicate names: {", ".join(duplicate_names[:3])}')

    if unrecognized_files:
        flash(f'Skipped {len(unrecognized_files):,} unrecognized file types: {", ".join(unrecognized_files[:3])}')

    if not antibody_inputs:
        flash('No input sequences were provided')
        return redirect(request.url)

    max_inputs = current_app.config['MAX_INPUTS']
    if len(antibody_inputs) > max_inputs:
        flash(f'Maximum number of input sequences exceeded, processed only the first {max_inputs} inputs')
        antibody_inputs = antibody_inputs[:max_inputs]

    return antibody_inputs


def convert_pdb_positions(pdb_chain: BioPdbChain, ab_chain: AbNumberChain) -> Dict[str, str]:
    positions = {}
    chain_seq = ab_chain.seq
    pdb_seq = "".join(seq1(r.get_resname()) for r in pdb_chain.get_residues())[:len(chain_seq)]
    assert pdb_seq == chain_seq, f'PDB sequence does not match SEQRES field: "{pdb_seq}" != "{chain_seq}"'
    for (pos, aa), pdb_residue in zip(ab_chain, pdb_chain.get_residues()):
        hetflag, resseq, icode = pdb_residue.get_id()
        positions[str(pos)] = f'{resseq}^{icode.strip()}:{pdb_chain.get_id()}'
    return positions


def clean_extension(filename):
    """Return clean file extension

    - lowercase
    - without .gz
    - without . prefix

    examples:
    file.fa.gz -> fa
    file.pdb -> pdb
    nosuffix -> ''
    """
    filename = filename.lower()
    if filename.endswith('.gz'):
        filename = filename[:-3]

    _, extension = os.path.splitext(filename)
    return extension[1:]


def clean_antibody_name(name):
    """Get clean antibody name from FASTA identifier

    - Remove chain type suffix such as "_VH" or "_VL"
    """
    for suffix in ['_VH', '_VL', '_HC', '_LC']:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


def pair_antibody_records(records) -> (List[AntibodyInput], List[str], List[str]):
    invalid_names = []
    duplicate_names = []

    heavy_protein_seqs = {}
    light_protein_seqs = {}
    for protein_record in records:
        name = clean_antibody_name(protein_record.id)
        if looks_like_dna(protein_record.seq):
            protein_record = protein_record.translate()

        try:
            if looks_like_antibody_heavy_chain(protein_record.seq):
                protein_seqs = heavy_protein_seqs
            else:
                protein_seqs = light_protein_seqs
        except ChainParseError:
            invalid_names.append(name)
            # TODO log record
            continue
        if name in protein_seqs:
            duplicate_names.append(name)
            continue
        protein_seqs[name] = str(protein_record.seq)

    antibody_inputs = []
    for name in sorted(set(heavy_protein_seqs) | set(light_protein_seqs)):
        antibody_inputs.append(AntibodyInput(
            name=name,
            heavy_protein_seq=heavy_protein_seqs.get(name),
            light_protein_seq=light_protein_seqs.get(name),
        ))
    return antibody_inputs, invalid_names, duplicate_names


def parse_antibody_pdb_ids(pdb_ids):
    antibody_inputs = []
    invalid_names = []
    for pdb_id in pdb_ids:
        antibody_input = AntibodyInput.from_pdb_id(pdb_id)
        if antibody_input is None:
            invalid_names.append(pdb_id)
            continue
        antibody_inputs.append(antibody_input)
    return antibody_inputs, invalid_names


def parse_antibody_inputs(seq_string: str,
                          files: List[Union[str, FileStorage]] = None,
                          pdb_ids: List[str] = None
                          ) -> (List[AntibodyInput], List[str], List[str], List[str]):
    files = files or []
    pdb_ids = pdb_ids or []

    txt_antibody_inputs, txt_invalid_names, txt_duplicate_names = pair_antibody_records(
        parse_plaintext_records(seq_string)
    )
    pdb_antibody_inputs, pdb_invalid_names = parse_antibody_pdb_ids(
        sorted(set(pdb_ids))
    )
    file_antibody_inputs, file_invalid_names, file_duplicate_names, file_unrecognized_files = parse_antibody_files(
        files
    )

    return (
        txt_antibody_inputs + pdb_antibody_inputs + file_antibody_inputs,
        txt_invalid_names + pdb_invalid_names + file_invalid_names,
        txt_duplicate_names + file_duplicate_names,
        file_unrecognized_files
    )


def read_file_contents(file: Union[str, FileStorage]):
    if isinstance(file, str):
        with open(file, 'r') as fp:
            contents = fp.read()
        filename = os.path.basename(file)
    elif isinstance(file, FileStorage):
        contents = file.stream.read().decode()
        filename = file.filename
    else:
        raise NotImplementedError(f'Expected FileStorage or path string, got {type(file)}: {file}')
    return contents, filename


def parse_antibody_files(files: List[Union[str, FileStorage]]) -> (List[AntibodyInput], List[str], List[str], List[str]):
    names = set()
    pdb_antibody_inputs = []
    pdb_invalid_names = []
    pdb_duplicate_names = []
    unrecognized_files = []
    unpaired_records = []
    for file in files:
        if not file:
            continue
        contents, filename = read_file_contents(file)
        extension = clean_extension(filename)
        if extension in ['pdb']:
            name = filename.split('.', 1)[0]
            if name in names:
                pdb_duplicate_names.append(name)
                continue
            antibody_input = AntibodyInput.from_pdb_data(
                pdb_data=contents,
                name=name
            )
            if antibody_input is None:
                pdb_invalid_names.append(name)
                continue
            pdb_antibody_inputs.append(antibody_input)
            names.add(name)
        elif extension in ['fa', 'fna', 'faa', 'fasta']:
            unpaired_records += list(SeqIO.parse(StringIO(contents), 'fasta'))
        else:
            unrecognized_files.append(filename)

    seq_antibody_inputs, seq_invalid_names, seq_duplicate_names = pair_antibody_records(unpaired_records)

    return (
        pdb_antibody_inputs + seq_antibody_inputs,
        pdb_invalid_names + seq_invalid_names,
        seq_duplicate_names,
        unrecognized_files
    )


def chunk_list(lst, n) -> List[List]:
    """Yield successive n-sized chunks from lst."""
    if hasattr(lst, '__getitem__'):
        # Use slicing
        for i in range(0, len(lst), n):
            yield lst[i:i + n]
    else:
        # Use naive iteration
        batch = []
        for val in lst:
            batch.append(val)
            if len(batch) == n:
                yield batch
                batch = []
        if len(batch):
            yield batch


def send_text(text, name, extension='txt', timestamp=True):
    bytesio = BytesIO()
    bytesio.write(text.encode())
    bytesio.seek(0)

    if timestamp:
        name += time.strftime("_%Y-%m-%d_%H-%M-%S")

    return send_file(
        bytesio,
        attachment_filename=f'{name}.{extension}',
        as_attachment=True,
        cache_timeout=0
    )


def send_fasta(records, name, timestamp=True):
    stringio = StringIO()
    SeqIO.write(records, stringio, 'fasta')

    return send_text(stringio.getvalue(), name=name, extension='fa', timestamp=timestamp)


def write_sheet(df, writer, sheet_name='Sheet1', index=True, **kwargs):
    """
    Write df as an excel file to ExcelWriter, roughly similar to `df.to_excel` except that it handles
    `df` with MultiIndex columns and `index=False`.
    """

    # Avoid the "NotImplementedError: Writing to Excel with MultiIndex columns"
    # exception by temporarily changing the columns to a single-level index
    columns: pd.Index = df.columns
    df.columns = range(len(df.columns))
    df.to_excel(
        writer, startrow=columns.nlevels, header=False, sheet_name=sheet_name, index=index, **kwargs
    )
    df.columns = columns

    # Get the xlsxwriter workbook and worksheet objects.
    book = writer.book
    sheet = writer.sheets[sheet_name]

    # Add a header format.
    header_format = book.add_format({
        'bold': True,
        'text_wrap': True,
        'valign': 'top',
        'fg_color': '#e8e8f9',
        'border': 1
    })

    col_offset = (df.index.nlevels if index else 0)

    # Write the column headers with the defined format.
    for level in range(columns.nlevels):
        names = columns.get_level_values(level)
        prev_name = None
        last_name_change_i = 0
        for i, name in enumerate(names):
            if prev_name == name:
                sheet.merge_range(level, last_name_change_i + col_offset, level, i + col_offset, name, header_format)
                continue
            sheet.write(level, i + col_offset, name, header_format)
            last_name_change_i = i
            prev_name = name

    if index and any(df.index.names):
        for i, name in enumerate(df.index.names):
            sheet.write(df.columns.nlevels-1, i, name or '', header_format)
            sheet.write(df.columns.nlevels-1, i, name or '', header_format)


def write_sheets(df_dict: Dict[str, pd.DataFrame], fd_or_path):
    writer = pd.ExcelWriter(fd_or_path, engine='xlsxwriter')

    for sheet_name, df in df_dict.items():
        has_index_name = any(df.index.names)
        write_sheet(df, writer, sheet_name=sanitize_excel_sheet_name(sheet_name), index=has_index_name)

    writer.close()


def send_excel(df_dict: Dict[str, pd.DataFrame], name, timestamp=True):
    output = BytesIO()
    write_sheets(df_dict, output)
    output.seek(0)

    if timestamp:
        name += time.strftime("_%Y-%m-%d_%H-%M-%S")
    return send_file(
        output,
        attachment_filename=f'{name}.xlsx',
        as_attachment=True,
        cache_timeout=0
    )


def sanitize_excel_sheet_name(name):
    sanitized = name
    for char in '[]:*?/\\':
        sanitized = sanitized.replace(char, '_')
    if sanitized != name:
        print(f'Renaming excel sheet "{name}" to "{sanitized}"', file=sys.stderr)
    return sanitized
