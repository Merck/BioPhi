from typing import List

from abnumber.germlines import get_germline_v_families, get_germline_v_genes

from biophi.common.utils.io import (send_excel, send_fasta, send_text, read_antibody_input_request)
from biophi.common.utils.scheduler import scheduler
from biophi.humanization.methods.humanization import SapiensHumanizationParams, \
    CDRGraftingHumanizationParams, ManualHumanizationParams
from biophi.humanization.methods.humanness import (
    OASIS_MIN_SUBJECTS_THRESHOLDS, OASisParams)
from biophi.humanization.web.tasks import (HumanizeAntibodyTaskResult,
                                           HumannessTaskResult,
                                           humanize_antibody_task,
                                           humanness_task, mutate_humanized_antibody_task)
from flask import (Blueprint, current_app, redirect, render_template,
                   request, url_for)
from celery.exceptions import TimeoutError
from werkzeug import Response

biophi_humanization = Blueprint(
    'biophi_humanization',
    __name__,
    template_folder='templates',
    static_folder='static'
)


def _get_germline_lists():
    heavy_germlines = sorted(get_germline_v_families('H') + get_germline_v_genes('H'))
    light_germlines = sorted(get_germline_v_families('K') + get_germline_v_genes('K') + ['IGKV']
                             + get_germline_v_families('L') + get_germline_v_genes('L') + ['IGLV'])
    heavy_germlines = [n for n in heavy_germlines if '/' not in n]
    light_germlines = [n for n in light_germlines if '/' not in n]
    return heavy_germlines, light_germlines


@biophi_humanization.route('/humanize/', methods=['GET'])
def humanize_get():
    heavy_germlines, light_germlines = _get_germline_lists()

    humanization_method = 'sapiens'
    sapiens_humanization_params = SapiensHumanizationParams()
    cdr_grafting_humanization_params = CDRGraftingHumanizationParams()
    sequence_text = ''
    humanness_params = None
    scheme = None
    cdr_definition = None
    input_mode = 'single'

    humanize_task_id = request.args.get('humanize_task_id')
    if humanize_task_id:
        humanization_results: List[HumanizeAntibodyTaskResult] = scheduler.get_results(humanize_task_id)
        if 'index' in request.args:
            humanization_results = [humanization_results[int(request.args['index'])-1]]
        parental_records = [record for result in humanization_results for record in result.get_parental_records()]
        sequence_text = ''.join(r.format('fasta-2line') for r in parental_records)
        humanization_params = humanization_results[0].humanization_params
        humanness_params = humanization_results[0].oasis_params
        scheme = humanization_params.scheme
        cdr_definition = humanization_params.cdr_definition
        input_mode = 'bulk'
        if isinstance(humanization_params, SapiensHumanizationParams):
            sapiens_humanization_params = humanization_params
            humanization_method = 'sapiens'
        elif isinstance(humanization_params, CDRGraftingHumanizationParams):
            cdr_grafting_humanization_params = humanization_params
            humanization_method = 'cdr_grafting'
        else:
            raise NotImplementedError(str(type(humanization_params)))

    humanness_task_id = request.args.get('humanness_task_id')
    if humanness_task_id:
        humanness_results: List[HumannessTaskResult] = scheduler.get_results(humanness_task_id)
        if 'index' in request.args:
            humanness_results = [humanness_results[int(request.args['index'])-1]]
        records = [record for result in humanness_results for record in result.get_records()]
        sequence_text = ''.join(r.format('fasta-2line') for r in records)
        input_mode = 'bulk'

    return render_template(
        'humanization/humanize_input.html',
        input_mode=input_mode,
        scheme=scheme,
        cdr_definition=cdr_definition,
        sequence_text=sequence_text,
        humanness_params=humanness_params,
        humanization_method=humanization_method,
        sapiens_humanization_params=sapiens_humanization_params,
        cdr_grafting_humanization_params=cdr_grafting_humanization_params,
        heavy_germlines=heavy_germlines,
        light_germlines=light_germlines,
        OASIS_MIN_SUBJECTS_THRESHOLDS=OASIS_MIN_SUBJECTS_THRESHOLDS
    )


@biophi_humanization.route('/humanize/', methods=['POST'])
def humanize_post():
    antibody_inputs = read_antibody_input_request()

    if isinstance(antibody_inputs, Response):
        return antibody_inputs

    scheme = request.form['scheme']
    humanizing_cdr_definition = request.form['cdr_definition']
    if humanizing_cdr_definition == 'kabat_vernier':
        cdr_definition = 'kabat'
        backmutate_vernier = True
    else:
        cdr_definition = humanizing_cdr_definition
        backmutate_vernier = False

    method = request.form.get('method', 'manual')
    if method == 'manual':
        humanization_params = ManualHumanizationParams()
    elif method == 'sapiens':
        iterations = int(request.form['iterations'])
        assert iterations <= 10, 'Sapiens iteration limit exceeded'
        humanization_params = SapiensHumanizationParams(
            scheme=scheme,
            cdr_definition=cdr_definition,
            backmutate_vernier=backmutate_vernier,
            humanize_cdrs='humanize_cdrs' in request.form,
            model_version='latest',
            iterations=iterations
        )
    elif method == 'cdr_grafting':
        humanization_params = CDRGraftingHumanizationParams(
            scheme=scheme,
            cdr_definition=cdr_definition,
            backmutate_vernier=backmutate_vernier,
            heavy_v_germline=request.form['heavy_v_germline'],
            light_v_germline=request.form['light_v_germline'],
            sapiens_iterations=(1 if 'sapiens_final_pass' in request.form else 0)
        )
    else:
        raise ValueError(request.form['method'])

    min_fraction_subjects = OASIS_MIN_SUBJECTS_THRESHOLDS.get(
        request.form['min_subjects'],
        float(request.form.get('custom_min_subjects', 10)) / 100
    )

    oasis_params = OASisParams(
        oasis_db_path=current_app.config['OASIS_DB_PATH'],
        min_fraction_subjects=min_fraction_subjects
    )

    task_id = scheduler.schedule_tasks(
        humanize_antibody_task,
        inputs=[dict(
                    input=input,
                    humanization_params=humanization_params,
                    oasis_params=oasis_params
                ) for input in antibody_inputs]
    )

    if request.args.get('designer'):
        return redirect("/humanization/designer/" + scheduler.get_result_task_id(task_id, 0))

    return redirect("/humanization/humanize/results/" + task_id)


@biophi_humanization.route('/humanize/results/<task_id>/', methods=['GET'])
def humanize_results_get(task_id):
    show = request.args.get('show', 'table')
    limit = int(request.args.get('limit', 10))

    if not scheduler.are_results_ready(task_id):
        running, completed, total = scheduler.get_results_progress(task_id)
        return render_template(
            'loading.html',
            active_nav='humanize',
            running=running,
            completed=completed,
            total=total
        )

    results: List[HumanizeAntibodyTaskResult] = scheduler.get_results(task_id)

    successful_results = [r for r in results if not isinstance(r, Exception)]
    has_both_chain_types = any(r.humanization.vh for r in successful_results) \
                           and any(r.humanization.vl for r in successful_results)
    num_unpaired = sum((not r.humanization.vh or not r.humanization.vl) for r in successful_results)

    page = max(int(request.args.get('page', 1)), 1)
    has_prev_page = page > 1
    has_next_page = len(results) > page*limit
    offset = (page-1)*limit
    results = results[offset: offset+limit]

    return render_template(
        'humanization/humanize_results.html',
        task_id=task_id,
        results=results,
        has_both_chain_types=has_both_chain_types,
        num_unpaired=num_unpaired,
        page=page,
        show=show,
        has_prev_page=has_prev_page,
        has_next_page=has_next_page,
        offset=offset
    )


@biophi_humanization.route('/humanize/results/<task_id>/<index>', methods=['GET'])
def humanize_detail_get(task_id, index):
    result_task_id = scheduler.get_result_task_id(task_id, index)
    result: HumanizeAntibodyTaskResult = scheduler.get_result(task_id, index)
    return render_template('humanization/humanize_detail.html',
                           task_id=task_id, result_task_id=result_task_id, result=result, result_index=index)


@biophi_humanization.route('/humanize/report/<task_id>/humanized.fa', methods=['GET'])
def humanize_detail_export_humanized_fasta(task_id):
    index = request.args.get('index')
    result: HumanizeAntibodyTaskResult = scheduler.get_result(task_id, index)
    return send_fasta(
        result.get_humanized_records(),
        name=result.get_export_name(),
    )


@biophi_humanization.route('/humanize/report/<task_id>/alignment.txt', methods=['GET'])
def humanize_detail_export_alignment(task_id):
    index = request.args.get('index')
    result: HumanizeAntibodyTaskResult = scheduler.get_result(task_id, index)
    return send_text(
        result.humanization.get_alignment_string(),
        name=result.get_export_name(),
    )


@biophi_humanization.route('/humanize/report/<task_id>/oasis.xls', methods=['GET'])
def humanize_detail_export_oasis_table(task_id):
    index = request.args.get('index')
    result: HumanizeAntibodyTaskResult = scheduler.get_result(task_id, index)
    return send_excel(
        {
            'Humanized': result.humanized_humanness.to_peptide_dataframe(),
            'Parental': result.parental_humanness.to_peptide_dataframe()
        },
        name=result.get_export_name()
    )


@biophi_humanization.route('/humanize/report/batch/<task_id>/humanized.fa', methods=['GET'])
def humanize_batch_export_humanized_fasta(task_id):
    results: List[HumanizeAntibodyTaskResult] = scheduler.get_results(task_id)
    results = [r for r in results if not isinstance(r, Exception)]
    records = [record for result in results for record in result.get_humanized_records()]
    return send_fasta(records, name=results[0].get_export_name(num_seqs=len(results)))


@biophi_humanization.route('/humanize/report/batch/<task_id>/alignments.txt', methods=['GET'])
def humanize_batch_export_alignments(task_id):
    results: List[HumanizeAntibodyTaskResult] = scheduler.get_results(task_id)
    results = [r for r in results if not isinstance(r, Exception)]
    text = '\n\n'.join([result.humanization.get_alignment_string() for result in results])
    name = results[0].input.safe_name if len(results) == 1 else results[0].get_export_name(num_seqs=len(results))
    return send_text(text, name=name)


@biophi_humanization.route('/humanize/report/batch/<task_id>/sapiens.xlsx', methods=['GET'])
def humanize_batch_export_table(task_id):
    results: List[HumanizeAntibodyTaskResult] = scheduler.get_results(task_id)
    results = [r for r in results if not isinstance(r, Exception)]
    full = bool(request.args.get('full'))
    sheets = HumanizeAntibodyTaskResult.to_sheets(results, full=full)
    return send_excel(
        sheets,
        name=results[0].get_export_name(num_seqs=len(results))
    )


@biophi_humanization.route('/designer/', methods=['GET'])
def designer_get():
    humanness_task_id = request.args.get('humanness_task_id')
    if humanness_task_id:
        result: HumannessTaskResult = scheduler.get_result(humanness_task_id, request.args['index'])

        task_id = scheduler.schedule_task(
            humanize_antibody_task,
            input=result.input,
            oasis_params=result.oasis_params,
            humanization_params=ManualHumanizationParams()
        )
        return redirect(url_for('biophi_humanization.designer_detail_get', task_id=task_id))

    return render_template(
        'humanization/designer_input.html',
        input_mode='single',
        OASIS_MIN_SUBJECTS_THRESHOLDS=OASIS_MIN_SUBJECTS_THRESHOLDS
    )

@biophi_humanization.route('/designer/<task_id>', methods=['GET'])
def designer_detail_get(task_id):
    tab = request.args.get('tab')
    pos = request.args.get('pos')
    aa = request.args.get('aa')

    try:
        timeout = 15 if pos and aa else 5
        result: HumanizeAntibodyTaskResult = scheduler.get_result(task_id, timeout=timeout)
    except TimeoutError:
        return render_template('loading.html', active_nav='designer')

    if pos and aa:
        new_task_id = scheduler.schedule_task(mutate_humanized_antibody_task, result=result, pos=pos, aa=aa)
        return redirect(url_for('biophi_humanization.designer_detail_get', task_id=new_task_id, tab=tab))

    return render_template('humanization/designer_detail.html',
                           tab=tab,
                           task_id=task_id, result=result)


@biophi_humanization.route('/humanness/', methods=['GET'])
def humanness_get():
    return render_template('humanization/humanness_report_input.html',
                           input_mode='single',
                           OASIS_MIN_SUBJECTS_THRESHOLDS=OASIS_MIN_SUBJECTS_THRESHOLDS)


@biophi_humanization.route('/humanness/', methods=['POST'])
def humanness_post():
    antibody_inputs = read_antibody_input_request()
    if isinstance(antibody_inputs, Response):
        return antibody_inputs

    scheme = request.form['scheme']
    cdr_definition = request.form['cdr_definition']

    min_fraction_subjects = OASIS_MIN_SUBJECTS_THRESHOLDS.get(
        request.form['min_subjects'],
        float(request.form.get('custom_min_subjects', 10)) / 100
    )

    oasis_params = OASisParams(
        oasis_db_path=current_app.config['OASIS_DB_PATH'],
        min_fraction_subjects=min_fraction_subjects
    )

    task_id = scheduler.schedule_tasks(
        humanness_task,
        inputs=[dict(
                    input=input,
                    oasis_params=oasis_params,
                    scheme=scheme,
                    cdr_definition=cdr_definition
                ) for input in antibody_inputs]
    )

    return redirect("/humanization/humanness/report/" + task_id)


@biophi_humanization.route('/humanness/report/<task_id>/', methods=['GET'])
def humanness_report_table_get(task_id):
    limit = int(request.args.get('limit', 10))

    if not scheduler.are_results_ready(task_id):
        running, completed, total = scheduler.get_results_progress(task_id)
        return render_template(
            'loading.html',
            active_nav='humanness',
            running=running,
            completed=completed,
            total=total
        )

    results: List[HumannessTaskResult] = scheduler.get_results(task_id)

    page = max(int(request.args.get('page', 1)), 1)
    has_prev_page = page > 1
    has_next_page = len(results) > page*limit
    offset = (page-1)*limit
    results = results[offset: offset+limit]

    return render_template(
        'humanization/humanness_report_table.html',
        task_id=task_id,
        results=results,
        page=page,
        has_prev_page=has_prev_page,
        has_next_page=has_next_page,
        offset=offset
   )


@biophi_humanization.route('/humanness/report/<task_id>/detail/<int:result_index>/', methods=['GET'])
def humanness_report_detail_get(task_id, result_index):
    result: HumannessTaskResult = scheduler.get_result(task_id, result_index)
    return render_template('humanization/humanness_report_detail.html',
                           task_id=task_id, result_index=result_index, result=result)


@biophi_humanization.route('/humanness/report/<task_id>/oasis.xls', methods=['GET'])
def humanness_export_oasis_table(task_id):
    results: List[HumannessTaskResult] = scheduler.get_results(task_id)
    results = [r for r in results if not isinstance(r, Exception)]
    full = bool(request.args.get('full'))
    sheets = HumannessTaskResult.to_sheets(results, full=full)
    return send_excel(sheets, name='OASis_Results' if full else 'OASis_Summary')


@biophi_humanization.route('/humanness/report/<task_id>/detail/<int:result_index>/oasis.xls', methods=['GET'])
def humanness_detail_export_oasis_table(task_id, result_index):
    result: HumannessTaskResult = scheduler.get_result(task_id, result_index)
    return send_excel(
        {
            'Peptides': result.humanness.to_peptide_dataframe(),
        },
        name=result.input.safe_name
    )
