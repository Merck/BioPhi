from functools import partial
from multiprocessing import Pool

import click
from Bio import SeqIO
import pandas as pd
from biophi.common.utils.formatting import logo
from biophi.common.utils.io import parse_antibody_files, write_sheets
from biophi.common.utils.seq import iterate_fasta
from biophi.humanization.cli.oasis import show_unpaired_warning
from biophi.humanization.methods.humanness import OASisParams
from biophi.humanization.methods.humanization import humanize_chain, SapiensHumanizationParams, HumanizationParams
from abnumber import Chain, ChainParseError, SUPPORTED_CDR_DEFINITIONS, SUPPORTED_SCHEMES
import os
import sys
from tqdm import tqdm
from biophi.humanization.web.tasks import humanize_antibody_task, HumanizeAntibodyTaskResult


@click.command()
@click.argument('inputs', required=False, nargs=-1)
@click.option('--output', required=False, help='Output directory path. With --fasta-only, output FASTA file path.')
@click.option('--fasta-only', is_flag=True, default=False, type=bool, help='Output only a FASTA file with humanized sequences (speeds up processing)')
@click.option('--scores-only', is_flag=True, default=False, type=bool, help='Output only a CSV file with Sapiens position*residue scores')
@click.option('--mean-score-only', is_flag=True, default=False, type=bool, help='Output only a CSV file with one Sapiens score per sequence')
@click.option('--oas-db', required=False, help='OAS peptide database connection string (required to run OASis)')
@click.option('--version', default='latest', help='Sapiens trained model name')
@click.option('--iterations', type=int, default=1, help='Run Sapiens given number of times to discover more humanizing mutations')
@click.option('--scheme', default=HumanizationParams.cdr_definition, help=f'Numbering scheme: one of {", ".join(SUPPORTED_SCHEMES)}')
@click.option('--cdr-definition', default=HumanizationParams.cdr_definition, help=f'CDR definition: one of {", ".join(SUPPORTED_CDR_DEFINITIONS)}')
@click.option('--humanize-cdrs', is_flag=True, default=False, type=bool, help='Allow humanizing mutations in CDRs')
@click.option('--limit', required=False, metavar='N', type=int, help='Process only first N records')
def sapiens(inputs, output, fasta_only, scores_only, mean_score_only, version, iterations, scheme, cdr_definition, humanize_cdrs, limit, oas_db):
    """Sapiens: Antibody humanization using deep learning.

     Sapiens is trained on 20 million natural antibody sequences
     from antibody repertoires of more than 500 human subjects
     found in the Observed Antibody Space database

    EXAMPLES:

        \b
        # Humanize FASTA file(s), print to standard output
        biophi sapiens input.fa

        \b
        # Humanize FASTA file(s) without running OASis (fast)
        biophi sapiens input.fa --fasta-only --output humanized.fa

        \b
        # Humanize FASTA file(s), save to directory along with OASis humanness report
        biophi sapiens input.fa --output ./report/ \\
          --oas-db sqlite:////Absolute/path/to/oas_human_subject_9mers_2019_11.db

    INPUTS: Input FASTA file path(s). If not provided, creates an interactive session.
    """

    logo('''    ____              _                
   / ___|  __ _ _ __ (_) ___ _ __  ___ 
   \___ \ / _` | '_ \| |/ _ \ '_ \/ __|
    ___| | |_| | |_| | |  __/ | | \__ \\
   |____/ \__,_|  __/|_|\___|_| |_|___/
               |_| ''')

    click.echo(f'Settings:', err=True)
    click.echo(f'- Predicting using Sapiens model: {version}', err=True)
    if scores_only or mean_score_only:
        click.echo(f'- Producing Sapiens scores only', err=True)
        if fasta_only:
            raise ValueError('Cannot use --fasta-only together with --scores-only')
        if iterations != 1:
            raise ValueError('Iterations cannot be used with --scores-only')
        iterations = 0
    else:
        click.echo('- Mutations in CDR regions are enabled' if humanize_cdrs else '- Ignoring mutations in CDR regions', err=True)
        click.echo(f'- Humanizing using {iterations} {"iteration" if iterations == 1 else "iterations"}', err=True)
    click.echo(err=True)

    humanization_params = SapiensHumanizationParams(
        model_version=version,
        humanize_cdrs=humanize_cdrs,
        scheme=scheme,
        cdr_definition=cdr_definition,
        iterations=iterations
    )
    oasis_params = OASisParams(
        oasis_db_path=oas_db,
        min_fraction_subjects=0.10
    ) if oas_db else None

    if inputs:
        if scores_only or mean_score_only:
            return sapiens_scores_only(
                inputs,
                output,
                limit=limit,
                humanization_params=humanization_params,
                mean=mean_score_only
            )
        elif fasta_only:
            return sapiens_fasta_only(
                inputs,
                output,
                limit=limit,
                humanization_params=humanization_params
            )
        else:
            from biophi.common.web.views import app
            with app.app_context():
                return sapiens_full(
                    inputs,
                    output,
                    limit=limit,
                    humanization_params=humanization_params,
                    oasis_params=oasis_params
                )
    else:
        if output:
            click.echo('Warning! Ignoring --output and --report parameter in interactive session (no inputs provided)', err=True)
        return sapiens_interactive(scheme=scheme, humanization_params=humanization_params)


def sapiens_interactive(scheme, humanization_params):
    while True:
        seq = input('Enter sequence:\n')
        try:
            chain = Chain(seq, scheme=scheme, cdr_definition=humanization_params.cdr_definition)
        except ChainParseError as e:
            print(e)
            continue
        humanization = humanize_chain(chain, params=humanization_params)
        click.echo('Parental -> Humanized:')
        click.echo(humanization.alignment.format())


def sapiens_scores_only(inputs, output_csv, humanization_params, limit=None, mean=False):
    chains = [Chain(
                record.seq,
                name=record.id,
                scheme=humanization_params.scheme,
                cdr_definition=humanization_params.cdr_definition
              ) for record in iterate_fasta(inputs, limit=limit)]

    if output_csv:
        click.echo(f'Writing scores CSV to file: {output_csv}', err=True)
        f = open(output_csv, 'wt')
    else:
        click.echo('Writing to stdout', err=True)
        f = None

    headers = None
    for chain in tqdm(chains, disable=(not output_csv)):
        humanization = humanize_chain(chain, params=humanization_params)
        scores = humanization.to_score_dataframe()
        seq = chain.seq
        chain_label = 'H' if chain.is_heavy_chain() else 'L'

        if mean:
            scores_by_pos_and_aa = scores.melt(ignore_index=False).set_index('variable', append=True)['value']
            seq_scores = scores_by_pos_and_aa.loc[list(enumerate(seq))]
            mean_score = seq_scores.mean()
            # Return single row with mean score value for given chain
            scores = pd.DataFrame([{'id': chain.name, 'chain': chain_label, 'sapiens_score': mean_score}]).set_index('id')
        else:
            scores.insert(0, 'id', chain.name)
            scores.insert(1, 'chain', chain_label)
            scores.insert(2, 'input_aa', list(seq))

        if headers is None:
            scores.to_csv(f if output_csv else sys.stdout, header=True)
            headers = scores.columns
        else:
            assert headers.equals(scores.columns), f'Score dataframe columns do not agree: ' \
                                                   f'{headers.values} VS {scores.columns.values}'
            scores.to_csv(f if output_csv else sys.stdout, header=False)

    if output_csv:
        click.echo(f'Saved {len(chains):,} humanized chains to: {output_csv}', err=True)
        click.echo('Done.', err=True)


def sapiens_fasta_only(inputs, output_fasta, humanization_params, limit=None):
    chains = [Chain(
                record.seq,
                name=record.id,
                scheme=humanization_params.scheme,
                cdr_definition=humanization_params.cdr_definition
              ) for record in iterate_fasta(inputs, limit=limit)]

    if output_fasta:
        click.echo(f'Writing FASTA to file: {output_fasta}', err=True)
        f = open(output_fasta, 'wt')
    else:
        click.echo('Writing to stdout', err=True)
        f = None

    for chain in tqdm(chains, disable=(not output_fasta)):
        humanization = humanize_chain(chain, params=humanization_params)
        chain_label = 'VH' if humanization.humanized_chain.is_heavy_chain() else 'VL'
        Chain.to_fasta(
            humanization.humanized_chain,
            f if output_fasta else sys.stdout,
            description=f'{chain_label} {humanization_params.get_export_name().replace("_"," ")}'
        )

    if output_fasta:
        click.echo(f'Saved {len(chains):,} humanized chains to: {output_fasta}', err=True)
        click.echo('Done.', err=True)


def sapiens_full(inputs, output_dir, humanization_params, oasis_params, limit=None):
    if oasis_params is None:
        raise ValueError('OASis params need to be provided for full output, or consider using --fasta-only')
    if not os.path.exists(output_dir) or not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if len(os.listdir(output_dir)):
        raise FileExistsError(f'Output directory exists and is not empty: {output_dir}')

    try:
        antibody_inputs, invalid_names, duplicate_names, unrecognized_files = parse_antibody_files(
            files=inputs
        )

        if invalid_names:
            click.echo(f'Warning: Skipped {len(invalid_names):,} invalid inputs: {", ".join(invalid_names[:3])}', err=True)

        if duplicate_names:
            click.echo(f'Warning: Skipped {len(duplicate_names):,} inputs with duplicate names: {", ".join(duplicate_names[:3])}', err=True)

        if unrecognized_files:
            raise ValueError(f'Expected FASTA or PDB files, got: {unrecognized_files}')

        if not antibody_inputs:
            raise ValueError(f'No input sequences found')

        if limit:
            antibody_inputs = antibody_inputs[:limit]

        click.echo(f'Running Sapiens humanization on {len(antibody_inputs)} antibodies...', err=True)

        show_unpaired_warning(antibody_inputs)

        pool = Pool()
        iterator = pool.imap(
            partial(humanize_antibody_task, humanization_params=humanization_params, oasis_params=oasis_params),
            antibody_inputs
        )
        results = list(tqdm(iterator, total=len(antibody_inputs)))

        click.echo('Saving report...', err=True)
        alignments_path = os.path.join(output_dir, 'alignments.txt')
        fasta_path = os.path.join(output_dir, 'humanized.fa')
        excel_path = os.path.join(output_dir, 'Sapiens.xlsx')

        alignments = [result.humanization.get_alignment_string() for result in results]
        with open(alignments_path, 'wt') as fd:
            fd.write('\n\n'.join(alignments))
        click.echo(f'Saved alignments to: {alignments_path}', err=True)

        records = [record for result in results for record in result.get_humanized_records()]
        SeqIO.write(records, fasta_path, 'fasta')
        click.echo(f'Saved humanized chains to: {fasta_path}', err=True)

        sheets = HumanizeAntibodyTaskResult.to_sheets(results, full=True)
        write_sheets(sheets, excel_path)
        click.echo(f'Saved XLSX report to: {excel_path}', err=True)

        click.echo('Done.', err=True)

    except:
        # remove empty output directory
        if os.path.exists(output_dir) and os.path.isdir(output_dir) and not len(os.listdir(output_dir)):
            os.rmdir(output_dir)
        raise

