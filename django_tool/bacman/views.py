from django.shortcuts import render
import json     #mi serve per salvare EC che non printo a schermo
from .algorithms.step_one import primo
from .algorithms.step_two import secondo

# Create your views here.
def homepage(request):
    return render(request, 'bacman/bacman-webpage.html')

def pathways(request):
    if request.method == "POST":
        pollutant = request.POST.get("pollutant")
        final = request.POST.get("final")

    paths, ec = primo(pollutant, final)
    request.session['ec_list'] = json.dumps(ec)

    return render(request, 'bacman/pathways.html', { 'first_path': paths[0], 'second_path': paths[1], 'third_path': paths[2]})







def contact(request):
    return render(request, 'bacman/contactBacman.html')

def help(request):
    return render(request, 'bacman/helpBacman.html')