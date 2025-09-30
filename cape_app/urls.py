"""
URL configuration for the 'cape_app' application.

This module defines all the URL patterns for the application. Each URL is mapped
to a corresponding view function in `views.py`, which handles the request and
returns a response (rendering templates, processing data, or providing downloads).

Use this module to manage navigation between different pages of the application.
"""

from django.urls import path
from . import views

# List of URL patterns
urlpatterns = [
    # Homepage of the application
    path('', views.homepage, name='homepage'),

    # Pathways analysis page: compute and display metabolic pathways
    path('pathways', views.pathways, name='pathways'),

    # EC choice page: user selects EC numbers for each step of a pathway
    path('ec_choice', views.ec_choice, name='ec_choice'),

    # Orthologs page: process selected ECs and show available ortholog sequences
    path('orthologs', views.orthologs, name='orthologs'),

    # Summary page: display selected ortholog sequences for final review
    path('summary', views.summary, name='summary'),

    # Codon optimization page: perform codon optimization on selected protein sequences
    path('codon_optimization', views.codon_optimization, name='codon_optimization'),

    # Download FASTA page: allow downloading codon-optimized sequences
    path('download_fasta', views.download_fasta, name='download_fasta'),

    # Results page: additional results view
    path('results', views.results, name='results'),

    # Help page: provide user guidance and documentation
    path('help', views.help, name='help'),
]
