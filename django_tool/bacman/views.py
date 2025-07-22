from django.shortcuts import render
import json     #mi serve per salvare EC che non printo a schermo
from .algorithms.step_one import primo
from .algorithms.step_two import secondo

# Create your views here.
def homepage(request):
    return render(request, 'bacman/bacman-webpage.html')







def contact(request):
    return render(request, 'bacman/contactBacman.html')

def help(request):
    return render(request, 'bacman/helpBacman.html')