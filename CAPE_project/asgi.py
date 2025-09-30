"""
ASGI config for CAPE_project.

It exposes the ASGI callable as a module-level variable named ``application``.
This is used by any ASGI-compatible server to serve your project.

For more information, see:
https://docs.djangoproject.com/en/5.2/howto/deployment/asgi/
"""

import os
from django.core.asgi import get_asgi_application

# Set the default Django settings module for the project
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'CAPE_project.settings')

# ASGI application callable used by ASGI servers
application = get_asgi_application()