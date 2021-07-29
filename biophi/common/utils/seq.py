import urllib
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools
import subprocess
import numpy as np
from io import StringIO
from typing import List, Tuple
from Bio import SeqIO
from abnumber import Chain
import re


AMINO_ACIDS_INCLUDING_X = frozenset(list('ACDEFGHIKLMNPQRSTVWXY'))


def sanitize_sequence(seq):
    # remove whitespace
    seq = ''.join(seq.split())
    # remove stop codon at the end
    if seq.endswith('*'):
        seq = seq[:-1]
    return seq


def is_valid_amino_acid_sequence(seq):
    return all(aa in AMINO_ACIDS_INCLUDING_X for aa in seq)


def parse_plaintext_records(text) -> List[SeqRecord]:
    text = text.strip()
    if not text:
        return []
    if text[0] != '>':
        text = f'> Unnamed({text[:10]}...)\n{text}'
    records = list(SeqIO.parse(StringIO(text), 'fasta'))
    for record in records:
        assert len(record.seq) != 0, f'Sequence "{record.id}" is empty'
    return records


def iterate_fasta_index(records, idx):
    for i, record in enumerate(records):
        if i == idx[0]:
            yield record
            idx = idx[1:]
            if not len(idx):
                break


def iterate_single_fasta(path, limit=None, random=False, random_sample_from_limit=None):
    records = SeqIO.parse(path, 'fasta')
    if limit == 0:
        return []
    if random:
        if limit is None:
            raise ValueError('Random can only be used together with limit')
        if random_sample_from_limit:
            size = random_sample_from_limit
        else:
            assert "'" not in path
            process = subprocess.Popen(f"grep '^>' '{path}' | wc -l", stdout=subprocess.PIPE, shell=True)
            size = int(process.communicate()[0].strip())
        limit = min(size, limit)
        random_idx = np.array(sorted(np.random.choice(size, limit, replace=False)))
        records = iterate_fasta_index(records, random_idx)
    elif limit:
        records = itertools.islice(records, limit)

    return records


def iterate_fasta(paths, limit=None, random=False, random_sample_from_limit=None):
    """
    Iterate through fasta sequence file(s), return generator of records
    """
    if isinstance(paths, str):
        paths = [paths]
    for path in paths:
        for record in iterate_single_fasta(path, limit=limit, random=random,
                                           random_sample_from_limit=random_sample_from_limit):
            yield record


def download_pdb(pdb_id):
    return urllib.request.urlopen(f'https://files.rcsb.org/download/{pdb_id.strip()}.pdb').read()


def looks_like_antibody_heavy_chain(seq):
    """Return True if sequence looks like an antibody heavy chain, False if it looks like a light chain

    Throws ChainParseError if sequence is not an antibody chain
    """
    if seq.endswith('TVSS'):
        return True
    if seq.endswith('EIK') or seq.endswith('LTVL') \
            or seq.endswith('DIQ') or seq.startswith('EIV') or seq.startswith('DIV') or seq.startswith('DVV'):
        return False
    return Chain(seq, scheme='imgt').chain_type == 'H'


DNA_SEQ_REGEX = re.compile(r'[ACTGUNactgun]+')


def looks_like_dna(seq):
    if not isinstance(seq, str) and not isinstance(seq, Seq):
        raise NotImplementedError(f'Expected string or Seq, got {type(seq)}: {seq}')
    return DNA_SEQ_REGEX.fullmatch(str(seq))


# B Asn or Asp
# Z Gln or Glu
# J Leu or Ile
# U Selenocysteine (UGA)
# O Pyrrolysine (UAG)
# X Unknown
PROTEIN_AMBIGUOUS_SEQ_REGEX = re.compile(r'[ABCDEFGHIJKLMNOPQRSTUVWXYZ]+')


def looks_like_protein(seq):
    if not isinstance(seq, str) and not isinstance(seq, Seq):
        raise NotImplementedError(f'Expected string or Seq, got {type(seq)}: {seq}')
    return PROTEIN_AMBIGUOUS_SEQ_REGEX.fullmatch(str(seq))


def validate_dna(seq):
    if seq is None:
        return
    if not looks_like_dna(seq):
        if seq == '':
            raise ValueError(f'DNA sequence can not be empty')
        raise ValueError(f'DNA sequence contains invalid characters: "{seq}"')


def validate_protein(seq):
    if seq is None:
        return
    if not looks_like_protein(seq):
        if seq == '':
            raise ValueError(f'Protein sequence can not be empty')
        if '*' in seq:
            raise ValueError(f'Protein sequence contains stop codons: "{seq}"')
        raise ValueError(f'Protein sequence contains invalid characters: "{seq}"')

