"""
Application configuration for the 'cape_app' Django application.

This class is used by Django to configure application-specific settings
and metadata, such as the default primary key field type.
"""

from django.apps import AppConfig


class CapeAppConfig(AppConfig):
    # Default primary key field type for models in this app
    default_auto_field = 'django.db.models.BigAutoField'

    # Name of the application
    name = 'cape_app'