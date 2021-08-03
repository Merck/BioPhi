import time

from celery import Celery
from celery.signals import after_task_publish, task_prerun, task_postrun, task_failure

from biophi.common.utils.stats import log_task_result

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


@task_prerun.connect
def update_running_state(task_id, task, **kwargs):
    backend = task.backend if task else celery.backend
    backend.store_result(task_id, None, "RUNNING")
    task.start_time = time.time()


@task_postrun.connect
def log_task_postrun(task_id, task, **kwargs):
    from biophi.common.web.views import app
    queued_seconds = 0
    running_seconds = time.time() - task.start_time
    with app.app_context():
        log_task_result(queued_seconds=queued_seconds, running_seconds=running_seconds)


@task_failure.connect
def log_task_failure(task_id, exception, **kwargs):
    from biophi.common.web.views import app
    with app.app_context():
        log_task_result(exception=exception)

