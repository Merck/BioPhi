import re

import click
from biophi import __version__


def spacer(length=60, **kwargs):
    click.echo('─' * length, **kwargs)


def logo(sublogo=''):
    logo = '''      __     ____  _       ____  _     _ 
  /| /  \   | __ )(_) ___ |  _ \| |__ (_)
 ( || [] )  |  _ \| |/ _ \| |_) | '_ \| |
  \_   _/   | |_) | | (_) |  __/| | | | |
    | |     |____/|_|\___/|_|   |_| |_|_|
    |_|              {}'''.format(f'version {__version__}'.rjust(20))

    if sublogo:
        assert sublogo.count('\n') == logo.count('\n'), 'Logos need to have same amount of lines'
        logo = '\n'.join([a + b for a, b in zip(logo.splitlines(), sublogo.splitlines())])
    click.echo(logo, err=True)
    spacer(length=len(logo.splitlines()[-1]) + 2, err=True)
    click.echo(err=True)


AA_NAMES = {
    'A': 'Alanine',
    'R': 'Arginine',
    'N': 'Asparagine',
    'D': 'Aspartic acid',
    'C': 'Cysteine',
    'Q': 'Glutamine',
    'E': 'Glutamic acid',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'L': 'Leucine',
    'K': 'Lysine',
    'M': 'Methionine',
    'F': 'Phenylalanine',
    'P': 'Proline',
    'S': 'Serine',
    'T': 'Threonine',
    'W': 'Tryptophan',
    'X': 'Unknown',
    'Y': 'Tyrosine',
    'V': 'Valine',
    '-': '(Gap)'
}


def aa_name(aa):
    return AA_NAMES[aa]


def get_valid_filename(s):
    s = str(s).strip().replace(' ', '_')
    return re.sub(r'[^-._a-zA-Z0-9]', '', s)


def human_size(bytes, units=(' bytes', 'kB', 'MB', 'GB', 'TB', 'PB', 'EB')):
    """ Returns a human readable string representation of bytes """
    return str(bytes) + units[0] if bytes < 1024 else human_size(bytes >> 10, units[1:])
