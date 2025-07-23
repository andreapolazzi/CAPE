from django.shortcuts import render
import json     #mi serve per salvare EC che non printo a schermo
from .algorithms.step_one import primo
from .algorithms.step_two import secondo

#HOMEPAGE
def homepage(request):
    return render(request, 'bacman/bacman-webpage.html')

#PATHWAYS PAGE: the 3 best pathways are presented, user can choose one
def pathways(request):
    if request.method == "POST":
        pollutant = request.POST.get("pollutant")
        final = request.POST.get("final")

    paths, ec = primo(pollutant, final)
    request.session['ec_list'] = json.dumps(ec)

    return render(request, 'bacman/pathways.html', { 'first_path': paths[0], 'second_path': paths[1], 'third_path': paths[2]})

#ORTHOLOGS PAGE: fasta files of all orthologs are presented, user chooses
def aaa(request):
    if request.method == "POST":
        selected_index = int(request.POST.get("selected_pathway"))
        ec = json.loads(request.session.get('ec_list', '[]'))

        selected_ec = ec[selected_index]
        df = secondo(selected_ec)

        #da capire come mettere altra colonna

        #table = df.to_html(index=True)

        return render(request, "bacman/aaa.html", {"table": df})       #PROVA SOLO PER VEDERE SE VA




def contact(request):
    return render(request, 'bacman/contactBacman.html')

def help(request):
    return render(request, 'bacman/helpBacman.html')