import pandas as pd

def parse_reaction_pathways(file_path='./results/results.csv'):
    """
    Parse the reaction pathways from a CSV file, identifying paths from
    source compounds at iteration 0 to sink compounds.
    
    Args:
        file_path (str): Path to the CSV file containing reaction data
        
    Returns:
        pd.DataFrame: A dataframe containing all the pathways
    """
    # Read the CSV data
    df = pd.read_csv(file_path)
    
    # Dictionary to track all compounds and their relationships
    compounds = {}
    
    # Process all reactions to build the network
    for _, row in df.iterrows():
        substrate = row['Substrate SMILES']
        substrate_inchi = row['Substrate InChI']
        products = row['Product SMILES'].split('.')
        product_inchis = row['Product InChI'].split('.')
        transformation = row['Transformation ID']
        iteration = int(row['Iteration'])
        in_sink = row['In Sink'] == 1
        sink_name = row['Sink name'] if in_sink else None
        rule_id = row['Rule ID']
        ec_number = row['EC number']
        score = row['Score']
        
        # Add substrate to compounds dictionary
        if substrate not in compounds:
            compounds[substrate] = {
                'inchi': substrate_inchi,
                'iterations': {iteration},
                'precursors': [],
                'products': [],
                'in_sink': False,
                'sink_name': None
            }
        else:
            compounds[substrate]['iterations'].add(iteration)
        
        # Process each product
        for i, product in enumerate(products):
            product_inchi = product_inchis[i] if i < len(product_inchis) else ""
            
            # Add product to compounds dictionary
            if product not in compounds:
                compounds[product] = {
                    'inchi': product_inchi,
                    'iterations': {iteration + 1},  # Product appears in next iteration
                    'precursors': [(substrate, transformation, rule_id, ec_number, score)],
                    'products': [],
                    'in_sink': in_sink,
                    'sink_name': sink_name
                }
            else:
                compounds[product]['iterations'].add(iteration + 1)
                if (substrate, transformation, rule_id, ec_number, score) not in compounds[product]['precursors']:
                    compounds[product]['precursors'].append((substrate, transformation, rule_id, ec_number, score))
                if in_sink:
                    compounds[product]['in_sink'] = True
                    compounds[product]['sink_name'] = sink_name
            
            # Update substrate's products list
            compounds[substrate]['products'].append((product, transformation, rule_id, ec_number, score))
    
    # Find sink compounds
    sink_compounds = [(smiles, info['sink_name']) 
                     for smiles, info in compounds.items() 
                     if info['in_sink']]
    
    # Function to trace pathway from sink back to iteration 0
    def trace_pathway(current, visited=None, path=None):
        if visited is None:
            visited = set()
        if path is None:
            path = []
        
        # Avoid cycles
        if current in visited:
            return []
        
        visited.add(current)
        current_path = path + [current]
        
        # If compound is in iteration 0, we found a complete path
        if 0 in compounds[current]['iterations']:
            return [current_path]
        
        # Otherwise, continue tracing back through precursors
        paths = []
        for precursor, _, _, _, _ in compounds[current]['precursors']:
            new_paths = trace_pathway(precursor, visited.copy(), current_path)
            paths.extend(new_paths)
        
        return paths
    
    # Collect all pathways
    all_paths = []
    path_id = 1
    
    for sink_smiles, sink_name in sink_compounds:
        # Trace all pathways from this sink compound back to iteration 0
        pathways = trace_pathway(sink_smiles)
        
        for pathway in pathways:
            # Reverse pathway to go from iteration 0 to sink
            pathway.reverse()
            
            # Create rows for this pathway
            for i in range(len(pathway) - 1):  # For each step in the pathway
                current = pathway[i]
                next_compound = pathway[i+1]
                
                # Find the transformation between these compounds
                transformation = "Unknown"
                rule_id = "Unknown"
                ec_number = "Unknown"
                score = "Unknown"
                for prod, trans, rule, ec, sc in compounds[current]['products']:
                    if prod == next_compound:
                        transformation = trans
                        rule_id = rule
                        ec_number = ec
                        score = sc
                        break
                
                # Add this step to all_paths
                all_paths.append({
                    'Path_ID': f"Path_{path_id}",
                    'Sink_Name': sink_name,
                    'Step': i,
                    'Source_Compound': current,
                    'Source_InChI': compounds[current]['inchi'],
                    'Target_Compound': next_compound,
                    'Target_InChI': compounds[next_compound]['inchi'],
                    'Transformation': transformation,
                    'Rule_ID': rule_id,
                    'EC_Number': ec_number,
                    'Score': score,
                    'Source_Iteration': min(compounds[current]['iterations']),
                    'Target_Iteration': min(compounds[next_compound]['iterations'])
                })
            
            path_id += 1
    
    # Create the final dataframe
    pathway_df = pd.DataFrame(all_paths)
    
    return pathway_df

# Run the parser and return the dataframe
if __name__ == "__main__":
    try:
        pathways_df = parse_reaction_pathways()
        print(f"Successfully identified {pathways_df['Path_ID'].nunique()} distinct pathways")
        print("\nSummary of pathways dataframe:")
        print(f"Shape: {pathways_df.shape}")
        print("\nFirst few rows with all columns:")
        pd.set_option('display.max_columns', None)
        print(pathways_df.head())
        print("\nColumn names:")
        print(pathways_df.columns)
    except Exception as e:
        print(f"Error parsing pathways: {e}")