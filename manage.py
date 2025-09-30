#!/usr/bin/env python
"""
Django's command-line utility for administrative tasks.

This file acts as the entry point for running Django management commands,
such as:
    - runserver
    - migrate
    - createsuperuser
    - shell

It sets up the Django environment and delegates the commands to Django's
command-line execution system.
"""

import os
import sys


def main():
    """
    Main entry point for administrative tasks.

    Sets the default Django settings module and runs commands passed
    via the command line.
    """
    # Set the default settings module for the project
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'CAPE_project.settings')

    try:
        # Import Django's command-line execution function
        from django.core.management import execute_from_command_line
    
    except ImportError as exc:
        # If Django is not installed or environment not activated, raise informative error
        raise ImportError(
            "Couldn't import Django. Are you sure it's installed and "
            "available on your PYTHONPATH environment variable? Did you "
            "forget to activate a virtual environment?"
        ) from exc

    # Pass the command-line arguments to Django's execution system
    execute_from_command_line(sys.argv)


# Execute main() only if this file is run directly
if __name__ == '__main__':
    main()
