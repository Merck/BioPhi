from celery import group
from abc import ABC, abstractmethod
from biophi.common.web.tasks import celery
import uuid

class TaskNotFoundError(Exception):
    pass


class Scheduler(ABC):

    @abstractmethod
    def get_result(self, task_id, index=None, timeout=15):
        pass

    @abstractmethod
    def get_result_task_id(self, task_id, index):
        pass

    @abstractmethod
    def get_results(self, task_id, timeout=15):
        pass

    @abstractmethod
    def are_results_ready(self, task_id):
        pass

    @abstractmethod
    def get_results_progress(self, task_id):
        pass

    @abstractmethod
    def schedule_task(self, fun, *args, **kwargs):
        pass

    @abstractmethod
    def schedule_tasks(self, fun, inputs):
        pass


class CeleryScheduler(Scheduler):

    def get_celery_group_result(self, task_id):
        group_result = celery.GroupResult.restore(task_id)
        if group_result is None:
            raise TaskNotFoundError(task_id)
        return group_result

    def get_result(self, task_id, index=None, timeout=15):
        if index is not None:
            assert int(index) > 0, f'Index is 1-indexed, got {index}'
            return self.get_celery_group_result(task_id)[int(index) - 1].get(timeout=timeout)
        async_result = celery.AsyncResult(task_id)
        if async_result.state == 'PENDING':
            # We update status of all created tasks to SENT, so no task should be PENDING
            # PENDING tasks are actually missing tasks
            raise TaskNotFoundError(task_id)
        return async_result.get(timeout=timeout)

    def get_result_task_id(self, task_id, index):
        return self.get_celery_group_result(task_id)[int(index) - 1].id

    def get_results(self, task_id, timeout=15):
        return self.get_celery_group_result(task_id).join(timeout=timeout, propagate=False)

    def are_results_ready(self, task_id):
        return self.get_celery_group_result(task_id).ready()

    def get_results_progress(self, task_id):
        group_result = self.get_celery_group_result(task_id)
        return group_result.completed_count(), len(group_result)

    def schedule_task(self, fun, *args, **kwargs):
        async_result = fun.delay(*args, **kwargs)
        return async_result.id

    def schedule_tasks(self, fun, inputs):
        group_result = group([fun.s(**kwargs) for kwargs in inputs]).apply_async()
        group_result.save()
        return group_result.id


class SimpleInMemoryScheduler(Scheduler):
    results = {}

    def save_result(self, result):
        task_id = None
        while task_id is None or task_id in self.results:
            task_id = uuid.uuid4().hex
        self.results[task_id] = result
        return task_id

    def get_result(self, task_id, index=None, timeout=15):
        if task_id not in self.results:
            raise TaskNotFoundError(task_id)
        if index is not None:
            assert int(index) > 0, f'Index is 1-indexed, got {index}'
            return self.get_result(self.results[task_id][int(index)-1])
        return self.results[task_id]

    def get_result_task_id(self, task_id, index):
        task_ids = self.get_result(task_id)
        return task_ids[int(index)-1]

    def get_results(self, task_id, timeout=15):
        task_ids = self.get_result(task_id)
        return [self.get_result(i) for i in task_ids]

    def are_results_ready(self, task_id):
        return task_id in self.results

    def get_results_progress(self, task_id):
        return None, None

    def schedule_task(self, fun, *args, **kwargs):
        return self.save_result(fun(*args, **kwargs))

    def schedule_tasks(self, fun, inputs):
        task_ids = []
        for kwargs in inputs:
            task_ids.append(self.save_result(fun(**kwargs)))
        return self.save_result(task_ids)


class NotInitializedScheduler(Scheduler):
    def throw(self):
        raise NotImplementedError('Scheduler not initialized')

    def get_result_task_id(self, task_id, index):
        self.throw()

    def get_results(self, task_id, timeout=15):
        self.throw()

    def are_results_ready(self, task_id):
        self.throw()

    def get_results_progress(self, task_id):
        self.throw()

    def schedule_task(self, fun, *args, **kwargs):
        self.throw()

    def schedule_tasks(self, fun, inputs):
        self.throw()

    def get_result(self, task_id, index=None, timeout=15):
        self.throw()


class SchedulerProxy:
    wrapped = NotInitializedScheduler()

    def __getattr__(self, attr):
        return getattr(self.wrapped, attr)


scheduler: Scheduler = SchedulerProxy()


def use_scheduler(name):
    global scheduler
    if name == 'celery':
        scheduler.wrapped = CeleryScheduler()
    elif name == 'simple':
        scheduler.wrapped = SimpleInMemoryScheduler()
    else:
        raise ValueError(f"Unsupported scheduler: {name}")