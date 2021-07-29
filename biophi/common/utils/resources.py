import os
import biophi

def get_resource_path(name, module):
    return os.path.join(os.path.dirname(biophi.__file__), module, 'resources', name)