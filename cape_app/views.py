"""
Views for the 'cape_app' application.

This module contains all view functions responsible for handling HTTP requests,
processing data, and rendering templates for the application. Each function
corresponds to a page or action within the CAPE workflow, from pathway computation
to codon optimization and results.
"""

# Standard library
import random
from collections import defaultdict
import os
import glob


# Django
from django.http import HttpResponse
from django.shortcuts import render
from django.utils.safestring import mark_safe
from django.http import JsonResponse, FileResponse, HttpResponse
from django.views.decorators.http import require_http_methods
from django.conf import settings

# Local algorithms
from .algorithms.pathways.main import main
from .algorithms.orthologs.processing import process_ec_list
from .algorithms.codon_opt.backtranslate_pd630 import backtranslate_protein, write_fasta
from .algorithms.selenzyme.noec_api import query_smarts

# Color map needed to paint ec numbers to identify their origin
COLOR_MAP = {'green': '#21b57c',
             'white': '#f2c94c',
             'RetroPath': '#f2784c',
             'NOEC': '#c726b9'}


def homepage(request):
    """
    Render the homepage of the application.

    Supported method: GET (default browser request)

    Output:
        HttpResponse rendering the HomeCAPE.html
    """

    # Render the template
    return render(request, 'cape/HomeCAPE.html')


def pathways(request):
    """
    Display the best n metabolic pathways from pollutant to compound.
    Calls the algorithm 'main' and formats results for visualization.
    
    Supported method: POST

    Input parameters - via POST:
        pollutant (str): starting compound to degrade
        sink (str): final compound to reach
        p_inchi (str): InChi code of the pollutant
        f_inchi (str): InChi code of the sink compound
        max_paths (int): number of pathways to return

    Session Data - write:
        pathways: raw pathways data
        database: dictionary containing information about compounds and enzymes
        pollutant: starting compound to degrade
        sink: final compound to reach
        requested_paths: number of pathways the user asked for

    Output - HttpResponse rendering the template pathways_html:
        context:
            pathways: formatted pathways ready to be displayed
            pollutant: starting compound to degrade
            sink: final compound to reach
    """

    if request.method == "POST":
        
        # Retrieve data from the submitted form
        pollutant = request.POST.get("pollutant")
        p_inchi = request.POST.get("p_inchi").strip()
        sink = request.POST.get("sink")
        f_inchi = request.POST.get("f_inchi")

        if not (p_inchi.startswith("InChI=") and f_inchi.startswith("InChI=")):
            return HttpResponse("Error: InChI are requested in the form 'InChI=...", status=400)

        # Determine how many pathways the user asked for, default to 3 on invalid input
        n = request.POST.get("max_paths")
        try:
            max_paths = int(n)
            if max_paths <= 0:
                raise ValueError
        except (TypeError, ValueError):
            max_paths = 3

        # Compute pathways using the function 'main'
        pathways = main(p_inchi, f_inchi, pollutant, sink, max_paths)
        
        # Save the data in session
        request.session['pathways'] = pathways
        request.session['pollutant'] = pollutant
        request.session['sink'] = sink
        request.session['requested_paths'] = max_paths
        
        # Build HTML representation of each pathway in the form:
        path_text = []
        for weight, steps in pathways:
            summary = f'WEIGHT: {weight:.2f}<br>'
            
            for i, step in enumerate(steps, 1):
                start_compound, end_compound = step[0], step[2]
                
                default_enzyme = step[1][0] if step[1] else {}
                enzyme_ec = default_enzyme.get('EC')
                
                # Map enzyme color
                if enzyme_ec == 'NOEC':
                    coloring = COLOR_MAP['NOEC']
                else:
                    coloring = COLOR_MAP.get(default_enzyme.get('color'), '#333333')
                
                colored_enzyme = f'<span style="color:{coloring}">{enzyme_ec}</span>'
                
                # Construct the step line
                summary += f"{i}. {start_compound} → {colored_enzyme} → {end_compound}<br>"
                
            path_text.append(mark_safe(summary))
        
        # Prepare context for rendering the template
        context = {'pathways': enumerate(path_text),
                   'pollutant': pollutant,
                   'sink': sink,
                   'requested_paths': max_paths}

        # Render the template with the data 'context'
        return render(request, 'cape/pathways.html', context )


def ec_choice(request):
    """
    Display the EC choice page where the user can customize EC assignments
    for each step of the selected metabolic pathway.
    
    Supported method: POST
    
    Input parameters - via POST:
        selected_pathway (int): index of the pathway chosen by the user
        
    Session data - read:
        pathways: raw pathways data
        database: dictionary containing information about compounds and enzymes
        pollutant: starting compound to degrade
        sink: final compound to reach

    Session data - write:
        selected_steps: detailed information about each step, including compounds,
                        EC numbers, enzyme info, and coloring
    
    Output - HttpResponse rendering the template ec_choice.html:
        context:
            steps: list of steps with detailed compound/enzyme information
            pollutant: starting compound to degrade
            sink: final compound to reach
    """

    if request.method == 'POST':

        # Retrieve the selected pathway index from the form
        selected_index = int(request.POST.get('selected_pathway'))

        # Load data from session
        pathways = request.session.get('pathways')
        pollutant = request.session.get('pollutant')
        sink = request.session.get('sink')

        # Extract the selected pathway and its steps
        _, steps = pathways[selected_index]

        detailed_steps = []
        non_green = len(steps)
        for i, step in enumerate(steps):
            start_compound, enzyme_list, end_compound = step

            # Collect all possible EC number for this step
            available_ecs = []
            step_type = None
            
            for enzyme in enzyme_list:
                enzyme_ec = enzyme.get('EC', '').strip().upper()
                enzyme_color = enzyme.get('color')
                enzyme_info = {k: v for k, v in enzyme.items() if k not in ["EC", "color"]}

                # Map enzyme color
                if enzyme_ec == 'NOEC':
                    coloring = COLOR_MAP['NOEC']
                else:
                    coloring = COLOR_MAP.get(enzyme_color, '#333333')

                colored_ec = f'<span style="color:{coloring}">{enzyme_ec}</span>'

                enzyme_entry = {'colored_ec': colored_ec,
                                'color': enzyme_color,
                                'EC': enzyme_ec,
                                **enzyme_info}
                
                if enzyme_ec == 'NOEC' and 'Reaction SMILES' in enzyme_info:
                    smarts = str(enzyme_info['Reaction SMILES'])
                    out = query_smarts(smarts)
                    data = out[smarts]

                    enzyme_entry['data'] = data
                    enzyme_entry['has_data'] = not data.empty
                    
                    enzyme_entry['ec_selenzyme'] = data['EC Number'].dropna().unique().tolist()

                    step_type = 'selenzyme'
                else:
                    enzyme_entry['has_data'] = False


                available_ecs.append(enzyme_entry)
            
            
            if step_type is None and available_ecs:
                if available_ecs[0]['color'] == 'green':
                    step_type = 'green'
                    non_green -= 1
                    available_ecs = [available_ecs[0]]
                else:
                    step_type = 'choice'

            
            # Build the detailed step representation
            detailed_steps.append({'step_number': i,
                                   'start_compound': start_compound,
                                   'end_compound': end_compound,
                                   'ecs': available_ecs,
                                   'step_type': step_type})

        for step in detailed_steps:
            for enzyme_entry in step['ecs']:
                if hasattr(enzyme_entry.get('data'), 'to_dict'):
                    enzyme_entry['data'] = enzyme_entry['data'].to_dict(orient='records')

        # Save detailed steps in session   
        request.session['selected_steps'] = detailed_steps

        # Prepare context for rendering the template
        context = {'steps': detailed_steps,
                   'pollutant': pollutant,
                   'sink': sink,
                   'non_green': non_green}
        
        # Render the template with the data 'context'
        return render(request, 'cape/ec_choice.html', context)


def orthologs(request):
    """
    Process orthologs for selected EC numbers and display results.
    ECs not already marked as 'green' are analyzed to provide available ortholog sequences.

    Supported method: POST

    Input parameters - via POST:
        selected_ec_<step_number> (str): EC number chosen for each step
        manual_ec_<step_number> (str, optional): manually entered EC number if "manual" option selected

    Session data - read:
        selected_steps: list of detailed steps with compounds and ECs
        pollutant: starting compound to degrade
        sink: final compound to reach

    Session data - write:
        orthologs: list of ortholog search results for selected ECs, including summary, sequence data, and hidden metadata

    Output - HttpResponse rendering the template orthologs.html:
        context:
            pollutant: starting compound to degrade
            sink: final compound to reach
            results: cleaned list of ortholog search results
            ecs: list of non-green ECs
    """

    if request.method == "POST":

        # Retrieve detailed steps and context from session
        detailed_steps = request.session.get('selected_steps')
        pollutant = request.session.get('pollutant')
        sink = request.session.get('sink')

        # List to collect EC numbers that require ortholog search
        non_green_ecs = []

        # Iterate through all steps to collect selected EC numbers
        for step in detailed_steps:
            step_number = step['step_number']
            selected_value = request.POST.get(f'selected_ec_{step_number}', '').strip()

            if selected_value:

                # User provided a manual EC number
                if selected_value == 'manual':
                    manual_ec = request.POST.get(f'manual_ec_{step_number}', '').strip()
                    if manual_ec:
                        non_green_ecs.append(manual_ec)
                
                else:
                    matching_ec = None
                    for ec_entry in step['ecs']:
                        if ec_entry.get('EC') == selected_value:
                            matching_ec = ec_entry

                        if 'ec_selenzyme' in ec_entry and selected_value in ec_entry['ec_selenzyme']:
                            matching_ec = ec_entry
                    
                    # Append EC number only if not already marked as 'green'
                    if matching_ec and matching_ec.get('color') != 'green':
                        non_green_ecs.append(selected_value)

        # Run the ortholog search algorithm
        raw_result = process_ec_list(non_green_ecs)

        # Clean results and prepare them for rendering
        cleaned_result = []
        for item in raw_result:
            ec_number = item.get('ec_number')
            summary = item.get('summary')
            data = item.get('data')
            hidden = item.get('hidden_data')

            # Convert pandas DataFrames to plain dictionaries (if present)
            if hasattr(data, 'to_dict'):
                data = data.to_dict(orient = 'records')

            if hasattr(hidden, 'to_dict'):
                hidden = hidden.to_dict(orient = 'records')

            cleaned_result.append({'ec_number': ec_number,
                                   'summary': summary,
                                   'data': data,
                                   'hidden': hidden})

        # Prepare context for rendering the template    
        context = {'pollutant': pollutant,
                   'sink': sink,
                   'results': cleaned_result,
                   'ecs': non_green_ecs}
        
        # Save ortholog results in session
        request.session['orthologs'] = cleaned_result

        # Render the template with the data 'context'
        return render(request, 'cape/orthologs.html', context)


def summary(request):
    """
    Display a summary of the selected ortholog sequences for each EC number.
    Collects user selections and stores them in the session for downstream processing.

    Supported method: POST

    Input parameters - via POST:
        selections (list of str): each selection in the format "ec_number||organism||header||sequence||length"

    Session data - read:
        orthologs: list of ortholog search results from previous step

    Session data - write:
        selected_sequences: list of dictionaries with detailed sequence information

    Output - HttpResponse rendering the template summary.html:
        context:
            summary_ordered: list of EC numbers with corresponding selected sequences
    """

    if request.method == "POST":

        # Retrieve the list of selected sequences submitted by the user
        selections = request.POST.getlist('selections')

        # Load all orthologs from session
        all_ec = request.session.get('orthologs', [])

        # Parse selections into structured dictionaries grouped by EC number
        parsed_data = defaultdict(list)
        for line in selections:
            # Each selection is expected in the format: ec_number||organism||header||sequence||length
            ec_number, organism, header, sequence, length = line.split('||')
            parsed_data[ec_number].append({'ec': ec_number,
                                           'organism': organism,
                                           'header': header,
                                           'sequence': sequence,
                                           'length': length})

        # Order summary according to orthologs previously computed
        summary_ordered = []
        for item in all_ec:
            ec_num = item.get('ec_number')
            rows = parsed_data.get(ec_num, [])
            summary_ordered.append({'ec_number': ec_num, 'rows': rows})

        # Flatten selected sequences for use in downstream steps
        selected_sequences = []
        for ec in summary_ordered:
            selected_sequences.extend(ec['rows'])

        fasta_entries = []
        for seq in selected_sequences:
            header = f">{seq['ec']}|{seq['organism']}|{seq['header']}"
            sequence = seq['sequence']
            wrapped_sequence = "\n".join(sequence[i:i+60] for i in range(0, len(sequence), 60))
            fasta_entry = f"{header}\n{wrapped_sequence}\n"
            fasta_entries.append(fasta_entry)

        fasta_output = "\n".join(fasta_entries)

        # Save selected sequences in session for later steps
        request.session['selected_sequences'] = selected_sequences
        request.session['faa_output'] = fasta_output

        # Prepare context for rendering the template 
        context = {'summary_ordered': summary_ordered}

        # Render the template with the data 'context'
        return render(request, 'cape/summary.html', context)


def codon_optimization(request):
    """
    Perform codon optimization for the selected protein sequences.
    Allows the user to choose an optimization method, provide restriction sites, and optionally specify a random seed for reproducibility.

    Supported method: POST

    Input parameters - via POST:
        optimization_method (str): method to use for codon optimization
        restriction_sites (list of str): predefined restriction sites to avoid
        custom_restriction_site (str, optional): custom restriction sites, comma-separated
        seed (int, optional): random seed for codon selection

    Session data - read:
        selected_sequences: list of dictionaries containing protein sequences and metadata

    Session data - write:
        fasta_output: codon-optimized sequences in FASTA format

    Output - HttpResponse rendering the template codon_optimization.html:
        context:
            fasta_output: optimized sequences in FASTA format
            method: chosen optimization method
            restriction_sites: final list of restriction sites used
    """

    fasta_output = None
    optimization_method = ''
    final_restriction_sites = []
    
    if request.method == 'POST':

        # Retrieve user parameters from POST
        optimization_method = request.POST.get('optimization_method')
        restriction_sites = request.POST.getlist('restriction_sites')
        custom_sites = request.POST.get('custom_restriction_site', '').strip()
        seed = request.POST.get('seed')

        # Initialize random number generator for reproducibility
        rng = random.Random(seed)

        # Build the final list of restriction sites
        final_restriction_sites = []
        for site in restriction_sites:
            # Split comma-separated entries into individual sites
            if ',' in site:
                final_restriction_sites.extend([s.strip() for s in site.split(',') if s.strip()])
            else:
                final_restriction_sites.append(site.strip())

        # Add any custom restriction sites provided by the user
        if custom_sites:
            custom_sites_list = [site.strip().upper() for site in custom_sites.split(',') if site.strip()]
            
            for site in custom_sites_list:
                if not all(c in 'ACGT' for c in site):
                    return HttpResponse(f"Error: Invalid restriction site '{site}'. Must contain only A,C,G,T characters.", status=400)
    
            final_restriction_sites.extend(custom_sites_list)

        # Retrieve the selected protein sequences from session
        selected_sequences = request.session.get('selected_sequences', [])
        
        fasta_entries = []
        for seq_info in selected_sequences:
            # Optimize the protein sequence into a DNA sequence
            seq = backtranslate_protein(seq_info['sequence'], optimization_method,  rng, final_restriction_sites)
            header = f"{seq_info['ec']}|{seq_info['organism']}|{seq_info['header']}"

            fasta_entries.append((header, seq))

        write_fasta(fasta_entries, './cape_app/pipeline_results/codon_opt_results/opt_genes_pd630.fna')

        with open('./cape_app/pipeline_results/codon_opt_results/opt_genes_pd630.fna', 'r') as f:
            fasta_output = f.read()

        # Save the optimized sequences in session for downstream use
        request.session['fasta_output'] = fasta_output
    
    # Prepare context for rendering the template
    context = {'fasta_output': fasta_output,
               'method': optimization_method,
               'restriction_sites': final_restriction_sites}
    
    # Render the template with the data 'context'
    return render(request, 'codon_optimization.html', context)

  
def download_fasta(request):
    """
    Allow the user to download the codon-optimized sequences as a FASTA file.

    Supported method: GET (or POST)

    Session data - read:
        fasta_output: codon-optimized sequences in FASTA (.fna) format

    Output - HttpResponse:
        - If FASTA is available: a downloadable file named 'opt_genes_pd630.fna'
        - If FASTA is not available: a plain text message indicating the user should run codon optimization first
    """
    file_type = request.GET.get('type', 'dna')

    if file_type == 'protein':
        fasta_output = request.session.get('faa_output')
        filename = "selected_proteins.faa"
    else:
        fasta_output = request.session.get('fasta_output')
        filename = "opt_genes_pd630.fna"

    # If no FASTA is available, inform the user
    if not fasta_output:
        return HttpResponse("No FASTA available. Please run codon optimization first.")

    # Return FASTA as a downloadable file
    response = HttpResponse(fasta_output, content_type='text/plain')
    response['Content-Disposition'] = f'attachment; filename="{filename}"'
    return response

@require_http_methods(["GET"])
def get_optimized_sequences(request):
    """Read and return the optimized sequences"""
    try:
        opt_file = os.path.join(settings.BASE_DIR, 'cape_app', 'pipeline_results', 'codon_opt_results', 'opt_genes_pd630.fna')
        if os.path.exists(opt_file):
            with open(opt_file, 'r') as f:
                content = f.read()
            return JsonResponse({'success': True, 'content': content})
        else:
            return JsonResponse({'success': False, 'error': 'Optimized sequences file not found'}, status=404)
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)}, status=500)

@require_http_methods(["GET"])
def get_hero_parts(request):
    """Read and return all HERO toolkit parts"""
    try:
        hero_parts = []
        
        # Read from rhodo_toolkit folder
        rhodo_path = os.path.join(settings.BASE_DIR, 'export', 'rhodo_toolkit')
        if os.path.exists(rhodo_path):
            for fasta_file in sorted(glob.glob(os.path.join(rhodo_path, '*.fasta'))):
                with open(fasta_file, 'r') as f:
                    hero_parts.append(f.read())
        
        # Read from constitutive_prom folder
        prom_path = os.path.join(settings.BASE_DIR, 'export', 'constitutive_prom')
        if os.path.exists(prom_path):
            for fasta_file in sorted(glob.glob(os.path.join(prom_path, '*.fasta'))):
                with open(fasta_file, 'r') as f:
                    hero_parts.append(f.read())

        # Read from primers folder
        primers_path = os.path.join(settings.BASE_DIR, 'export', 'primers')
        if os.path.exists(primers_path):
            for fasta_file in sorted(glob.glob(os.path.join(primers_path, '*.fasta'))):
                with open(fasta_file, 'r') as f:
                    hero_parts.append(f.read())

        # Read from reporter_gene folder
        reporter_path = os.path.join(settings.BASE_DIR, 'export', 'reporter_gene')
        if os.path.exists(reporter_path):
            for fasta_file in sorted(glob.glob(os.path.join(reporter_path, '*.fasta'))):
                with open(fasta_file, 'r') as f:
                    hero_parts.append(f.read())
        
        content = '\n\n'.join(hero_parts)
        return JsonResponse({'success': True, 'content': content})
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)}, status=500)

@require_http_methods(["GET"])
def get_backbone(request):
    """Read and return the pLoxship backbone sequence"""
    try:
        backbone_file = os.path.join(settings.BASE_DIR, 'export', 'composite_part', 'pLoxship.fa')
        if os.path.exists(backbone_file):
            with open(backbone_file, 'r') as f:
                content = f.read()
            return JsonResponse({'success': True, 'content': content})
        else:
            return JsonResponse({'success': False, 'error': 'Backbone file not found'}, status=404)
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)}, status=500)

@require_http_methods(["GET"])
def download_optimized_sequences(request):
    """Download the optimized sequences file"""
    try:
        opt_file = os.path.join(settings.BASE_DIR, 'cape_app', 'pipeline_results', 'codon_opt_results', 'opt_genes_pd630.fna')
        if os.path.exists(opt_file):
            response = FileResponse(open(opt_file, 'rb'), as_attachment=True, filename='opt_genes_pd630.fna')
            return response
        else:
            return JsonResponse({'error': 'File not found'}, status=404)
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)

@require_http_methods(["GET"])
def download_hero_parts(request):
    """Download HERO toolkit parts as a single file"""
    try:
        hero_parts = []
        
        # Read from rhodo_toolkit folder
        rhodo_path = os.path.join(settings.BASE_DIR, 'export', 'rhodo_toolkit')
        if os.path.exists(rhodo_path):
            for fasta_file in sorted(glob.glob(os.path.join(rhodo_path, '*.fasta'))):
                with open(fasta_file, 'r') as f:
                    hero_parts.append(f.read())
        
        # Read from constitutive_prom folder
        prom_path = os.path.join(settings.BASE_DIR, 'export', 'constitutive_prom')
        if os.path.exists(prom_path):
            for fasta_file in sorted(glob.glob(os.path.join(prom_path, '*.fasta'))):
                with open(fasta_file, 'r') as f:
                    hero_parts.append(f.read())

        # Read from primers folder
        primers_path = os.path.join(settings.BASE_DIR, 'export', 'primers')
        if os.path.exists(primers_path):
            for fasta_file in sorted(glob.glob(os.path.join(primers_path, '*.fasta'))):
                with open(fasta_file, 'r') as f:
                    hero_parts.append(f.read())

        # Read from reporter_gene folder
        reporter_path = os.path.join(settings.BASE_DIR, 'export', 'reporter_gene')
        if os.path.exists(reporter_path):
            for fasta_file in sorted(glob.glob(os.path.join(reporter_path, '*.fasta'))):
                with open(fasta_file, 'r') as f:
                    hero_parts.append(f.read())
        
        content = '\n\n'.join(hero_parts)
        response = HttpResponse(content, content_type='text/plain')
        response['Content-Disposition'] = 'attachment; filename="HERO_toolkit_parts.fna"'
        return response
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)

@require_http_methods(["GET"])
def download_backbone(request):
    """Download the pLoxship backbone"""
    try:
        backbone_file = os.path.join(settings.BASE_DIR, 'export', 'composite_part', 'pLoxship.fa')
        if os.path.exists(backbone_file):
            response = FileResponse(open(backbone_file, 'rb'), as_attachment=True, filename='pLoxship_backbone.fa')
            return response
        else:
            return JsonResponse({'error': 'File not found'}, status=404)
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)
    
def export(request):
    return render(request, 'cape/export.html')

def help(request):
    return render(request, 'cape/HelpCAPE.html')
