"""
WSGI config for CAPE_project.

It exposes the WSGI callable as a module-level variable named ``application``.
This is used by any WSGI-compatible server to serve your project.

For more information, see:
https://docs.djangoproject.com/en/5.2/howto/deployment/wsgi/
"""

import os
from django.core.wsgi import get_wsgi_application

# Set the default Django settings module for the project
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'CAPE_project.settings')

# WSGI application callable used by WSGI servers
application = get_wsgi_application()