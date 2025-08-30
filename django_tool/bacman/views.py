from django.shortcuts import render
import json     #mi serve per salvare EC che non printo a schermo

from collections import defaultdict

from .algorithms.pathways.main import main
from .algorithms.orthologs.processing import process_ec_list

#HOMEPAGE
def homepage(request):
    return render(request, 'bacman/HomeCAPE.html')


#PATHWAYS PAGE: the 3 best pathways are presented, user can choose one
def pathways(request):
    if request.method == "POST":
        pollutant = request.POST.get("pollutant")
        final = request.POST.get("final")
        inchi = request.POST.get("inchi")

    paths_old, ec = main(pollutant, final, inchi)
    request.session['ec_list'] = json.dumps(ec)

    paths = []
    for item in paths_old:
        paths.append(item.replace('\n', '<br>'))

    return render(request, 'bacman/pathways.html', { 'first_path': paths[0], 'second_path': paths[1], 'third_path': paths[2], 'pollutant': pollutant, 'sink':final})


#ORTHOLOGS PAGE: fasta files of all orthologs are presented, user chooses
def orthologs(request):
    if request.method == "POST":
        selected_index = int(request.POST.get("selected_pathway"))
        ec = json.loads(request.session.get('ec_list', '[]'))
        selected_ec = ec[selected_index]

        raw_result = process_ec_list(selected_ec)

        cleaned_result = []
        for item in raw_result:
            ec_number = item.get('ec_number')
            summary = item.get('summary')
            data = item.get('data')

            if hasattr(data, "to_dict"):
                data = data.to_dict(orient='records')

            cleaned_result.append({
                "ec_number": ec_number,
                "summary": summary,
                "data": data,
            })

        request.session['orthologs'] = json.dumps(cleaned_result)

        return render(request, 'bacman/orthologs.html', { 'results': cleaned_result })


def summary(request):
    if request.method == "POST":
        selections = request.POST.getlist('selections')
        parsed_data = defaultdict(list)

        for line in selections:
            ec_number, organism, header, sequence, length = line.split("||")
            parsed_data[ec_number].append({
                "ec": ec_number,
                "organism": organism,
                "header": header,
                "sequence": sequence,
                "length": length,
            })

        return render(request, 'bacman/summary.html', {'grouped_selections': dict(parsed_data)})

def confirm(request):
    return render(request, 'bacman/confirm.html')


def help(request):
    return render(request, 'bacman/HelpCAPE.html')