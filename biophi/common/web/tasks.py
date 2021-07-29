from celery import Celery
from celery.signals import after_task_publish

celery = Celery('biophi')
config_obj = 'biophi.common.web.celery_config'
celery.config_from_object(config_obj)

# Register tasks
import biophi.humanization.web.tasks

class ContextTask(celery.Task):
    def __call__(self, *args, **kwargs):
        from biophi.common.web.views import app
        with app.app_context():
            return self.run(*args, **kwargs)

celery.Task = ContextTask


@after_task_publish.connect
def update_sent_state(sender=None, headers=None, **kwargs):
    task = celery.tasks.get(sender)
    backend = task.backend if task else celery.backend
    backend.store_result(headers['id'], None, "SENT")
