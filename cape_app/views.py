from django.shortcuts import render
from django.utils.safestring import mark_safe
from django.http import HttpResponse
from collections import defaultdict
import random

from .algorithms.pathways.main import main
from .algorithms.orthologs.processing import process_ec_list
from .algorithms.codon_opt.backtranslate_pd630 import backtranslate_protein


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
        num = request.POST.get("number")

        #print(pollutant, p_inchi, sink, f_inchi)

        pathways, db = main(pollutant, sink, p_inchi, f_inchi)

        request.session['pathways'] = pathways
        request.session['database'] = db
        request.session['pollutant'] = pollutant
        request.session['sink'] = sink
        
        path_text = []
        for path in pathways:
            weight, steps = path[0], path[1]
            
            summary = f'WEIGHT: {weight:.2f}<br>'
            
            for i, step in enumerate(steps, 1):
                start_compound, end_compound = step[0], step[2]
                enzyme_ec, enzyme_color = step[1][0][0], step[1][0][1]

                coloring = COLOR_MAP.get(enzyme_color, '#333333')
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
        db = request.session.get('database')
        pollutant = request.session.get('pollutant')
        sink = request.session.get('sink')

        selected_pathway = pathways[selected_index]
        steps = selected_pathway[1]

        '''
        noec_steps = []
        for i, step in enumerate(steps):
            first_ec = step[1][0][0] if step[1] else None
            if first_ec == 'noec':
                noec_steps.append([first_ec, i]) #da capire cosa mi passa samu
        '''

        

        detailed_steps = []
        for i, step in enumerate(steps):
            start_compound, end_compound = step[0], step[2]

            available_ecs = []
            for ec in step[1]:
                enzyme_id, enzyme_color = ec[0], ec[1]
                enzyme_info = db.get(enzyme_id, {})

                coloring = COLOR_MAP.get(enzyme_color, '#333333')

                available_ecs.append({'ec_number': enzyme_id,
                                      'color': enzyme_color,
                                      'enzyme_name': enzyme_info.get('enzyme', 'Unknown enzyme'),
                                      'reactions': enzyme_info.get('reactions', []),
                                      'coloring': coloring})
            
            start_info = db.get(start_compound.lower(), {})
            end_info = db.get(end_compound.lower(), {})

            detailed_steps.append({'step_number': i,
                                   'start_compound': start_compound,
                                   'end_compound': end_compound,
                                   's_info': start_info,
                                   'e_info': end_info,
                                   'ecs': available_ecs})
            
        request.session['selected_steps'] = detailed_steps

        context = {'steps': detailed_steps,
                   'pollutant': pollutant,
                   'sink': sink}

        return render(request, 'bacman/ec_choice.html', context)

#ORTHOLOGS PAGE: fasta files of all orthologs are presented, user chooses
def orthologs(request):
    if request.method == "POST":
        detailed_steps = request.session.get('selected_steps')
        pollutant = request.session.get('pollutant')
        sink = request.session.get('sink')

        non_green_ecs = []

        for step in detailed_steps:
            step_n = step['step_number']

            selected_key = f'selected_ec_{step_n}'
            selected_value = request.POST.get(selected_key, '').strip()

            if selected_value:
                if selected_value == 'manual':
                    manual_key = f'manual_ec_{step_n}'
                    manual_ec = request.POST.get(manual_key, '').strip()

                    if manual_ec:
                        non_green_ecs.append(manual_ec)
                
                else:
                    matching_ec = next((ec for ec in step['ecs'] if ec['ec_number'] == selected_value), None)

                    if matching_ec and matching_ec['color'] != 'green':
                         non_green_ecs.append(selected_value)
                                       
        
        raw_result = process_ec_list(non_green_ecs)

        cleaned_result = []
        for item in raw_result:
            ec_number = item.get('ec_number')
            summary = item.get('summary')
            data = item.get('data')
            hidden = item.get('hidden_data')

            if hasattr(data, 'to_dict'):
                data = data.to_dict(orient = 'records')

            if hasattr(hidden, 'to_dict'):
                hidden = hidden.to_dict(orient = 'records')

            cleaned_result.append({'ec_number': ec_number,
                                   'summary': summary,
                                   'data': data,
                                   'hidden': hidden})
            
        context = {'pollutant': pollutant,
                   'sink': sink,
                   'results': cleaned_result,
                   'ecs': non_green_ecs}
        
        request.session['orthologs'] = cleaned_result

        return render(request, 'bacman/orthologs.html', context)

#SUMMARY PAGE: to confirm sequences chosen
def summary(request):
    if request.method == "POST":
        selections = request.POST.getlist('selections')
        all_ec = request.session.get('orthologs', [])

        parsed_data = defaultdict(list)
        

        for line in selections:
            ec_number, organism, header, sequence, length = line.split('||')
            parsed_data[ec_number].append({'ec': ec_number,
                                           'organism': organism,
                                           'header': header,
                                           'sequence': sequence,
                                           'length': length})

        summary_ordered = []
        for item in all_ec:
            ec_num = item.get('ec_number')
            rows = parsed_data.get(ec_num, [])
            summary_ordered.append({'ec_number': ec_num, 'rows': rows})

        selected_sequences = []
        for ec in summary_ordered:
            selected_sequences.extend(ec['rows'])
        request.session['selected_sequences'] = selected_sequences

        context = {'summary_ordered': summary_ordered}

        return render(request, 'bacman/summary.html', context)

#CODON OPTIMIZATION: choose max/weight and restriction
def codon_optimization(request):
    fasta_output = None
    optimization_method = ''
    final_restriction_sites = []
    
    if request.method == 'POST':
        # Recupero dei parametri dall'utente
        optimization_method = request.POST.get('optimization_method')
        restriction_sites = request.POST.getlist('restriction_sites')
        custom_sites = request.POST.get('custom_restriction_site', '').strip()
        seed = request.POST.get('seed')

        rng = random.Random(seed)

        final_restriction_sites = []
        
        for site in restriction_sites:
            if ',' in site:
                final_restriction_sites.extend([s.strip() for s in site.split(',') if s.strip()])
            else:
                final_restriction_sites.append(site.strip())

        if custom_sites:
            custom_sites_list = [site.strip() for site in custom_sites.split(',') if site.strip()]
            final_restriction_sites.extend(custom_sites_list)

        
        # Recupero delle sequenze selezionate dalla sessione o POST
        selected_sequences = request.session.get('selected_sequences', [])
        
        fasta_lines = []
        for seq_info in selected_sequences:
            seq = backtranslate_protein(seq_info['sequence'], optimization_method,  rng, final_restriction_sites)
            header = f">{seq_info['ec']}|{seq_info['organism']}|{seq_info['header']}"
            fasta_lines.append(header)
            for i in range(0, len(seq), 60):
                    fasta_lines.append(seq[i:i+60])
        
        fasta_output = "\n".join(fasta_lines)
        request.session['fasta_output'] = fasta_output
    
    context = {'fasta_output': fasta_output,
               'method': optimization_method,
               'restriction_sites': final_restriction_sites}
    
    return render(request, 'codon_optimization.html', context)

#To download fasta   
def download_fasta(request):
    fasta_output = request.session.get('fasta_output')

    if not fasta_output:
        return HttpResponse("No FASTA available. Please run codon optimization first.")

    response = HttpResponse(fasta_output, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="optimized_sequences.fasta"'
    return response


#PLASMID: build the plasmid
def plasmid(request):
    return render(request, 'bacman/plasmid.html')

#RESULTS: final results, plasmid completed, final file
def results(request):
    return render(request, 'bacman/results.html')

def help(request):
    return render(request, 'bacman/HelpCAPE.html')