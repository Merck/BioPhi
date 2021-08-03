import hashlib

from flask import current_app, request
from sqlalchemy import create_engine
import os
import pandas as pd
import datetime
import time

from sqlalchemy.exc import OperationalError


def get_stats(table='access_log', past_days=7):
    e = get_engine()
    from_date = (datetime.datetime.now() - datetime.timedelta(days=past_days)).date()
    from_timestamp = time.mktime(from_date.timetuple())
    try:
        stats = pd.read_sql(f'SELECT * FROM {table} WHERE timestamp >= {from_timestamp}', e, parse_dates=['timestamp'])
    except OperationalError as e:
        raise ValueError('Stats tables are not filled, try submitting a task first.', e)
    stats = stats.set_index('timestamp')
    return stats


def log_submission(antibody_inputs, invalid_names, duplicate_names, unrecognized_files):
    if not current_app.config['STATS_DB_PATH']:
        return
    num_inputs_unpaired = sum((not i.heavy_protein_seq or not i.light_protein_seq) for i in antibody_inputs)
    unrecognized_extensions = set([f.split('.')[-1] for f in unrecognized_files])
    log_data(
        data=dict(
            input_hash=hashlib.md5(','.join([str(i) for i in antibody_inputs]).encode()).hexdigest(),
            num_inputs=len(antibody_inputs),
            num_inputs_unpaired=num_inputs_unpaired,
            num_inputs_heavy=sum(bool(i.heavy_protein_seq) for i in antibody_inputs),
            num_inputs_light=sum(bool(i.light_protein_seq) for i in antibody_inputs),
            num_invalid_seqs=len(invalid_names),
            num_duplicate_names=len(duplicate_names),
            unrecognized_extensions=','.join(unrecognized_extensions) if unrecognized_files else None
        ),
        table='submission_log'
    )


def log_task_result(running_seconds=None, exception=None):
    if not current_app.config['STATS_DB_PATH']:
        return
    log_data(
        data=dict(
            running_seconds=running_seconds,
            exception=str(exception) if exception else None
        ),
        table='task_log'
    )


def log_access(exception=None):
    if not current_app.config['STATS_DB_PATH']:
        return
    log_data(
        data=dict(
            method=request.method,
            path=request.path,
            exception=str(exception) if exception else None
        ),
        table='access_log'
    )


def log_data(data, table):
    e = get_engine()
    try:
        request_data = dict(
            referer=request.headers.get('Referer'),
            browser_hash=hashlib.md5(
                '{}_{}'.format(request.headers.get('User-Agent'), request.remote_addr).encode()
            ).hexdigest(),
            ip=request.remote_addr,
            endpoint_name=request.endpoint.split('.')[-1] if request.endpoint else None
        )
    except RuntimeError:
        request_data = {}

    pd.DataFrame([dict(
        timestamp=time.time(),
        **request_data,
        **data
    )]).set_index('timestamp').to_sql(table, e, if_exists='append')


def get_engine():
    db_path = current_app.config['STATS_DB_PATH']
    if db_path is None:
        raise ValueError('Configure STATS_DB_PATH env var to use statistics')
    return create_engine('sqlite:///' + os.path.abspath(db_path))

