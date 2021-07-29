import os

broker_url = os.environ.get('CELERY_BROKER_URL', 'redis://localhost:6379/0'),
result_backend = os.environ.get('CELERY_RESULT_BACKEND', 'redis://localhost:6379/0')
accept_content = ['pickle']
task_serializer = 'pickle'
result_serializer = 'pickle'
