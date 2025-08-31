from django.urls import path
from . import views

#Define a list of url patterns
urlpatterns = [
    path('', views.homepage, name='homepage'),
    path('help', views.help, name='help'),
    path('pathways', views.pathways, name='pathways'),
    path('orthologs', views.orthologs, name='orthologs'),
    path('summary', views.summary, name='summary'),
    path('confirm', views.confirm, name='confirm')
]