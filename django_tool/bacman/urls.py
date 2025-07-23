from django.urls import path
from . import views

#Define a list of url patterns
urlpatterns = [
    path('', views.homepage, name='homepage'),
    path('contact', views.contact, name='contact'),
    path('help', views.help, name='help'),
    path('pathways', views.pathways, name='pathways'),
    #path('orthologs', views.orthologs, name='orthologs'),
    path('aaa', views.aaa, name='aaa'),                 #PROVA
]