import click
import traceback
from biophi.common.utils.formatting import spacer
from biophi.humanization.cli.sapiens import sapiens
from biophi.humanization.cli.oasis import oasis
from biophi.common.cli.web import web

class MainGroup(click.Group):
    def __call__(self, *args, **kwargs):
        try:
            return self.main(*args, **kwargs)
        except Exception as e:
            message = e.args[0] if e.args else ''
            traceback.print_exc()
            spacer(err=True)
            click.echo(f'BioPhi failed with {type(e).__name__}' + (f': {message}' if message else ''))
            if len(e.args) > 1:
                spacer(err=True)
                for arg in e.args[1:]:
                    click.echo(arg, err=True)
                spacer(err=True)
            exit(1)


@click.group(cls=MainGroup)
def main():
    pass

# Register all commands
# Humanization & Humanness
main.add_command(sapiens)
main.add_command(oasis)
main.add_command(web)


if __name__ == '__main__':
    main()
