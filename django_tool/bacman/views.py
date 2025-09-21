from django.shortcuts import render
from django.utils.safestring import mark_safe
import json
from collections import defaultdict
import re

from .algorithms.pathways.main import main
from .algorithms.orthologs.processing import process_ec_list


COLOR_MAP = {'green': '#21b57c',
             'white': '#f2c94c',
             'RetroPath': '#f2784c'}

#HOMEPAGE
def homepage(request):
    return render(request, 'bacman/HomeCAPE.html')


#PATHWAYS PAGE: the 3 best pathways are presented, user can choose one
def pathways(request):
    if request.method == "POST":
        pollutant = request.POST.get("pollutant")
        p_inchi = request.POST.get("p_inchi")
        sink = request.POST.get("sink")
        f_inchi = request.POST.get("f_inchi")

        pathways, db = main(pollutant, sink, p_inchi, f_inchi)

        request.session['pathways'] = pathways
        request.session['database'] = db
        request.session['pollutant'] = pollutant
        request.session['sink'] = sink
        
        color_map = COLOR_MAP
        
        path_text = []
        for path in pathways:
            score, steps = path[0], path[1]
            
            summary = f'SCORE: {score:.2f}<br>'
            
            for i, step in enumerate(steps, 1):
                start_compound, end_compound = step[0], step[2]
                enzyme_ec, enzyme_color = step[1][0][0], step[1][0][1]

                coloring = color_map.get(enzyme_color, '#333333')
                colored_enzyme = f'<span style = "color:{coloring}"> {enzyme_ec} </span>'
                
                summary += f"{i}. {start_compound} → {colored_enzyme} → {end_compound}<br>"
                
            path_text.append(mark_safe(summary))
        
        context = {'pathways': enumerate(path_text),
                'pollutant': pollutant,
                'sink': sink}

        return render(request, 'bacman/pathways.html', context )

#EC CHOICE PAGE: user can customize ec
def ec_choice(request):
    if request.method == 'POST':
        selected_index = int(request.POST.get('selected_pathway'))

        pathways = request.session.get('pathways')
        db = request.session.get('db')
        pollutant = request.session.get('pollutant')
        sink = request.session.get('sink')

        selected_pathway = pathways[selected_index]
        steps = selected_pathway[1]

        color_map = COLOR_MAP


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