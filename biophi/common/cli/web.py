import click
from biophi.common.utils.formatting import logo
from biophi.common.utils.scheduler import use_scheduler
from biophi.common.web.views import app
import time

@click.command()
@click.option('--host', default='localhost', help='Server port')
@click.option('--port', type=int, default=5000, help='Server port')
def web(host, port):
    """Run the BioPhi web interface in simplified mode (without a processing queue)"""
    logo()
    click.echo('Note! This is a simplified BioPhi server!')
    click.echo('      See BioPhi documentation on how to run BioPhi with a processing queue.')
    click.echo('')
    time.sleep(1)
    click.echo('Starting BioPhi...')

    use_scheduler('simple')

    app.run(debug=True, host=host, port=port, processes=1, use_reloader=False)
