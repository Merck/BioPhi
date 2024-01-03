import click
from biophi.common.utils.formatting import logo
from biophi.common.utils.io import parse_antibody_files, write_sheets
from biophi.humanization.methods.humanness import OASisParams
from biophi.humanization.web.tasks import humanness_task, HumannessTaskResult
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial


@click.command()
@click.argument('inputs', required=True, nargs=-1)
@click.option('--output', required=False, help='Output XLSX report file path')
@click.option('--oasis-db', required=True, help='OAS peptide database connection string')
@click.option('--scheme', default='kabat', help='Numbering scheme (kabat, chothia, imgt, aho)')
@click.option('--cdr-definition', default='kabat', help='Numbering scheme (kabat, chothia, imgt, north)')
@click.option('--min-percent-subjects', default=10, type=float, help='Minimum percent of OAS subjects to consider peptide human')
@click.option('--summary', is_flag=True, default=False, type=bool, help='Only export a single summary sheet')
def oasis(inputs, output, oasis_db, scheme, cdr_definition, min_percent_subjects, summary):
    """OASis: Antibody humanness evaluation using 9-mer peptide search.

    OASis evaluates antibody humanness by searching all overlapping 9-mers
    in repertoires of more than 500 human subjects from the Observed Antibody Space database

    EXAMPLES:

        \b
        # Evaluate humanness from FASTA file(s), save OASis humanness report to directory
        biophi oasis input.fa --output ./report.xlsx \\
          --oasis-db path/to/OASis9mers.db

    INPUTS: Input FASTA file path(s)
    """
    logo('''     ___    _   ____  _
    / _ \  / \ / ___|(_)___
   | | | |/ _ \\\\___ \| / __|
   | |_| / ___ \___| | \__ \\
    \___/_/   \_\___/|_|___/
                           ''')

    assert 1 <= min_percent_subjects <= 90, '--min-percent-subjects should be between 1 and 90'

    if not output.endswith('.xlsx'):
        raise ValueError(f'The --output is a spreadsheet and should have an .xlsx extension')

    click.echo(f'Settings:', err=True)
    click.echo(f'- OASis database: {oasis_db}', err=True)
    click.echo(f'- Numbering scheme: {scheme}', err=True)
    click.echo(f'- CDR definition: {cdr_definition}', err=True)
    click.echo('', err=True)

    click.echo(f'Loading and numbering chains: {" ".join(inputs)}', err=True)

    antibody_inputs, invalid_names, duplicate_names, unrecognized_files = parse_antibody_files(
        files=inputs, verbose=True
    )

    if invalid_names:
        click.echo(f'Warning: Skipped {len(invalid_names):,} invalid inputs: {", ".join(invalid_names[:3])}', err=True)

    if duplicate_names:
        click.echo(f'Warning: Skipped {len(duplicate_names):,} inputs with duplicate names: {", ".join(duplicate_names[:3])}', err=True)

    if unrecognized_files:
        raise ValueError(f'Expected FASTA or PDB files, got: {unrecognized_files}')

    if not antibody_inputs:
        raise ValueError(f'No input sequences found')

    click.echo(f'Running OASis on {len(antibody_inputs)} antibodies...', err=True)

    show_unpaired_warning(antibody_inputs)

    if len(antibody_inputs) > 1 and not summary:
        click.echo('Use --summary to export only the summary sheet', err=True)

    oasis_params = OASisParams(
        oasis_db_path=oasis_db,
        min_fraction_subjects=min_percent_subjects/100
    )
    pool = Pool()
    iterator = pool.imap(partial(
        humanness_task_wrapper,
        oasis_params=oasis_params,
        scheme=scheme,
        cdr_definition=cdr_definition
    ), antibody_inputs)
    results = list(tqdm(iterator, desc='Generating results', total=len(antibody_inputs)))

    # Remove failed results
    results = [r for r in results if r is not None]

    click.echo('Saving report...', err=True)
    sheets = HumannessTaskResult.to_sheets(results, full=not summary)
    write_sheets(sheets, output)
    click.echo(f'Saved XLSX report to: {output}', err=True)


def humanness_task_wrapper(antibody_input, **kwargs):
    try:
        return humanness_task(antibody_input, **kwargs)
    except Exception as e:
        print('Error:', antibody_input.name, e, antibody_input.heavy_protein_seq, antibody_input.light_protein_seq)
        return None


def show_unpaired_warning(antibody_inputs):
    has_both_chain_types = any(i.heavy_protein_seq for i in antibody_inputs) \
                           and any(i.light_protein_seq for i in antibody_inputs)
    num_unpaired = sum((not i.heavy_protein_seq or not i.light_protein_seq) for i in antibody_inputs)

    if has_both_chain_types and num_unpaired:

        click.echo('Warning:', err=True)
        click.echo(f'Found {num_unpaired:,} unpaired ' + ('chains' if num_unpaired > 1 else 'chain'), err=True)
        click.echo('These will be processed and analyzed separately')
        click.echo('Please use the same fasta ID for heavy and light chain to report results for whole antibody'
                   '(no spaces, optionally with _VH and _VL suffix).', err=True)
        click.echo('', err=True)
