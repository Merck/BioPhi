import os

# Max file upload size in bytes
MAX_CONTENT_LENGTH = int(os.environ.get('MAX_CONTENT_LENGTH', 1 * 1024 * 1024))

# Maximum number of input antibodies to be processed
MAX_INPUTS = int(os.environ.get('MAX_INPUTS', 5000))

# Path to BioPhi usage statistics sqlite database
STATS_DB_PATH = os.environ.get('STATS_DB_PATH')

# Path to OASis sqlite database
OASIS_DB_PATH = os.environ.get('OASIS_DB_PATH')

# Destroy context on exception even in Debug mode, so that we can log the exception to the stats DB
# This makes sure that teardown_request is called in debug mode when getting an exception in the endpoint
PRESERVE_CONTEXT_ON_EXCEPTION = False

# Show newsletter popup at the bottom of landing page
# using provided user/newsletter ID (something like c7bcd0367a4cdbbfe2ef413ff/c5b4d513538dcdb730f72a9db)
MAILCHIMP_NEWSLETTER = os.environ.get('MAILCHIMP_NEWSLETTER')

# Optional banner HTML (will not be sanitized!)
BANNER_HTML = os.environ.get('BANNER_HTML')
