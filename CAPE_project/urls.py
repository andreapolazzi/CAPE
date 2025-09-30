"""
URL configuration for CAPE_project.

The `urlpatterns` list routes URLs to views.

For more information, see:
https://docs.djangoproject.com/en/5.2/topics/http/urls/
"""

from django.contrib import admin
from django.urls import path, include

urlpatterns = [
    # Admin site URLs
    path('admin/', admin.site.urls),

    # Include all URLs defined in the 'cape_app' application
    path('', include('cape_app.urls')),
]