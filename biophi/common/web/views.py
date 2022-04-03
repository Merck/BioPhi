import datetime
import json

import humanize
from Bio.SeqUtils import seq3
import pandas as pd
from markupsafe import Markup
from werkzeug.urls import url_encode
from biophi.common.utils.formatting import aa_name
from flask import Flask, url_for, request, render_template
import os
from biophi.common.utils.scheduler import use_scheduler, TaskNotFoundError
from biophi.common.utils.stats import get_stats, log_access
from biophi.humanization.web.views import biophi_humanization


app = Flask(__name__)

app.secret_key = 'gijo080)Q@%0h8q808th0018020ahofijvvi018a-b8n8881244o09g-fff221111ttgj09s'

app.config.from_object('biophi.common.web.app_config')

app.jinja_env.globals.update(aa_name=aa_name)
app.jinja_env.globals.update(seq3=seq3)
app.jinja_env.globals.update(python_set=set)
app.jinja_env.globals.update(zip=zip)
app.jinja_env.globals.update(sum=sum)
app.jinja_env.globals.update(enumerate=enumerate)
app.jinja_env.globals.update(isinstance=isinstance)
app.jinja_env.globals.update(Exception=Exception)
app.jinja_env.globals.update(list=list)
app.jinja_env.globals.update(int=int)
app.jinja_env.globals.update(sorted=sorted)
app.jinja_env.globals.update(min=min)
app.jinja_env.globals.update(max=max)
app.jinja_env.globals.update(BANNER_HTML=app.config['BANNER_HTML'])

app.register_blueprint(biophi_humanization, url_prefix='/humanization')

use_scheduler('celery')

#
# Configure handlers
#

if app.config['STATS_DB_PATH']:
    @app.teardown_request
    def teardown_request_log_stats(exception):
        if request.endpoint and 'static' in request.endpoint:
            return
        if request.endpoint in ['stats']:
            return
        log_access(exception=exception)


@app.errorhandler(404)
def error_404(e):
    return render_template('error_404.html'), 404


if not app.config['DEBUG']:
    @app.errorhandler(Exception)
    def error(e):
        return render_template('error.html', e=e), 500


@app.errorhandler(TaskNotFoundError)
def error_task_not_found(e):
    return render_template('error_task_not_found.html'), 404

#
# Define basic endpoints
#


@app.route('/')
def index():
    web_path = os.path.dirname(__file__)
    news_path = os.path.join(web_path, 'static', 'news.json')
    with open(news_path) as f:
        news = json.load(f)
    for item in news:
        item['date'] = datetime.datetime.fromisoformat(item['date'])
    oasis_enabled = app.config['OASIS_DB_PATH'] is not None
    mailchimp_newsletter = app.config['MAILCHIMP_NEWSLETTER']
    return render_template('index.html', news=news,
                           oasis_enabled=oasis_enabled, mailchimp_newsletter=mailchimp_newsletter)


@app.route('/stats')
def stats():
    if not app.config['STATS_DB_PATH']:
        raise ValueError('Set STATS_DB_PATH env var to enable stats')
    past_days = 7
    return render_template(
        'stats.html',
        access_log=get_stats('access_log', past_days=past_days),
        submission_log=get_stats('submission_log', past_days=past_days),
        task_log=get_stats('task_log', past_days=past_days),
        pd=pd
    )

#
# Define template functions and filters
#

@app.template_global()
def icon(name, s=16):
    href = url_for('static', filename='bootstrap-icons-1.0.0/bootstrap-icons.svg')
    return Markup(f'<svg class="bi" width="{s}" height="{s}" fill="currentColor"><use xlink:href="{href}#{name}"/></svg>')


@app.template_global()
def info_icon(text, filled=False, secondary=True, delay=100, s=16):
    name = 'info-circle-fill' if filled else 'info-circle'
    return Markup(f'<span class="{"text-secondary" if secondary else ""}" data-tooltip-classes="tooltip-wide" '
                  f'data-bs-delay="{delay}" data-bs-toggle="tooltip" title="{text}">{icon(name, s=s)}</span>')


@app.template_filter()
def naturalday(v):
    return humanize.naturalday(v)


@app.template_filter()
def naturaltime(v):
    return humanize.naturaltime(v)


@app.template_global()
def big_number_format(num, precision=1):
    if pd.isna(num) or num is None:
        return '?'
    if abs(num) >= 1_000_000:
        return '{:,.{precision}f}M'.format(num / 1_000_000, precision=precision)
    if abs(num) >= 1_000:
        return '{:.{precision}f}k'.format(num / 1_000, precision=precision)
    return '{:,}'.format(num)


@app.template_global()
def url_for_arg(relative=False, **new_values):
    """Add query arguments to current URL"""
    args = request.args.copy()

    for key, value in new_values.items():
        args[key] = value

    if relative:
        return '?{}'.format(url_encode(args))
    return '{}?{}'.format(request.path, url_encode(args))


@app.template_filter('autoversion')
def autoversion_filter(filename):
    # determining fullpath might be project specific
    bin_dir = os.path.dirname(os.path.realpath(__file__))
    fullpath = os.path.join(bin_dir, filename.lstrip('/'))
    try:
        timestamp = str(os.path.getmtime(fullpath))
    except OSError:
        return filename
    newfilename = "{0}?v={1}".format(filename, timestamp)
    return newfilename
