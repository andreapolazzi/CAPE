from django.shortcuts import render
import json
from collections import defaultdict
import re

from .algorithms.pathways.main import main
from .algorithms.orthologs.processing import process_ec_list

#HOMEPAGE
def homepage(request):
    return render(request, 'bacman/HomeCAPE.html')


#PATHWAYS PAGE: the 3 best pathways are presented, user can choose one
def pathways(request):
    if request.method == "POST":
        pollutant = request.POST.get("pollutant")
        p_inchi = request.POST.get("p_inchi")
        final = request.POST.get("final")
        f_inchi = request.POST.get("f_inchi")

    paths_old, ec = main(pollutant, final, p_inchi, f_inchi)
    request.session['ec_list'] = json.dumps(ec)

    color_map = {"KEGG green edge": "#21b57c",
                 "KEGG white edge": "#d6cb2d",
                 "RP predicted edge": "#f2784c"}
    
    paths = []
    
    for item in paths_old:
        x = ''
        for line in item.split("\n"):
            def color_edge(match):
                line_color = "#000000"
                for key, color in color_map.items():
                    if key in line:
                        line_color = color
                return f'<span style="color:{line_color}">{match.group(0)}</span>'
            
            cleaned_line = re.sub(r'\([^)]*\)', '', line)
            colored_line = re.sub(r'\[.*?\]', color_edge, cleaned_line)

            x += colored_line + '<br>'
        
        paths.append(x)

    return render(request, 'bacman/pathways.html', { 'first_path': paths[0], 'second_path': paths[1], 'third_path': paths[2], 'pollutant': pollutant, 'sink':final})

#EC CHOICE PAGE: user can customize ec
def ec_choice(request):
    return render(request, 'bacman/ec_choice.html' )

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

#SUMMARY PAGE: to confirm sequences chosen
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

        all_ec = json.loads(request.session.get('orthologs', '[]'))

        summary_ordered = []
        for item in all_ec:
            ec_num = item.get('ec_number')
            rows = parsed_data.get(ec_num, [])
            summary_ordered.append({'ec_number': ec_num, 'rows': rows})

        return render(request, 'bacman/summary.html', {'summary_ordered': summary_ordered})

#CODON OPTIMIZATION: choose max/weight and restriction
def codon_optimization(request):
    return render(request, 'bacman/codon_optimization.html')

#PLASMID: build the plasmid
def plasmid(request):
    return render(request, 'bacman/plasmid.html')

#RESULTS: final results, plasmid completed, final file
def results(request):
    return render(request, 'bacman/results.html')

def help(request):
    return render(request, 'bacman/HelpCAPE.html')