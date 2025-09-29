from django.urls import path
from . import views

'''List of url patterns'''
urlpatterns = [
    path('', views.homepage, name='homepage'),

    path('pathways', views.pathways, name='pathways'),
    path('ec_choice', views.ec_choice, name='ec_choice'),

    path('orthologs', views.orthologs, name='orthologs'),
    path('summary', views.summary, name='summary'),

    path('codon_optimization', views.codon_optimization, name='codon_optimization'),
    path('download_fasta', views.download_fasta, name='download_fasta'),
    
    path('plasmid', views.plasmid, name='plasmid'),
    path('results', views.results, name='results'),
    path('help', views.help, name='help'),
]