import pandas as pd

def results_to_df(results_path: str) -> pd.DataFrame:
    df = pd.read_csv(results_path)
    return df

"""
def extract_row_data(df, col_name, query):
    target_rows = df[df[col_name] == query]
    #print(target_row)

    if not target_rows.empty:
        results = []
        for _, row in target_rows.iterrows():
            results.append((
                row['Substrate SMILES'],
                row['Product SMILES'],
                row['EC number'],
                row['Iteration']
            ))
        return results
    else:
        return []
"""
 
def extract_row_data(df, queries):
    """
    queries: ((col_name1, cond1), (col_name2, cond2) ...)
    only checks for == cond
    """
    target_rows = df
    for query in queries:
        #print(query)
        col_name, cond = query
        target_rows = target_rows[target_rows[col_name] == cond]
        #print(target_rows)
    target_rows.drop_duplicates(inplace=True)
    if not target_rows.empty:
            results = []
            for _, row in target_rows.iterrows():
                results.append({
                    'Substrate SMILES': row['Substrate SMILES'],
                    'Product SMILES': row['Product SMILES'],
                    'EC number': row['EC number'],
                    'Iteration': row['Iteration']
                })
            return results    


def reactions_given_prodSMILES(df, prod_smiles, iteration):
    return extract_row_data(df, (('Iteration', iteration), ('Product SMILES', prod_smiles),))


def extend_path(path, df):
    #print(path)
    return [path + [reaction] for reaction in reactions_given_prodSMILES(df, path[-1]['Substrate SMILES'], path[-1]['Iteration']-1) if reaction not in path]

def find_pathways(df, sink_name):
    for sink in extract_row_data(df, (('Sink name', f'[{sink_name}]'),)):
        #print(sink)
        path_length = sink['Iteration']
        pathways = [[sink] + [reaction] for reaction in reactions_given_prodSMILES(df, sink['Substrate SMILES'], path_length-1)]
        for i in range(path_length-1):
            for j in range(len(pathways)):
                ext_paths = extend_path(pathways[0], df)
                del pathways[0]
                pathways.extend(ext_paths)
    return pathways



if __name__ == "__main__":
    # Example usage
    df = results_to_df("results/results.csv")
    #print(df.head())

    #data = extract_row_data(df, 'Sink name', '[beta-Carboxy-cis,cis-muconate]')
    #print(data)
    """
    paths = find_pathways(df, 'beta-Carboxy-cis,cis-muconate')
    for path in paths:
        print(path[::-1], "\n")
    """
    sinks_df = df[ df["In Sink"] != 0 ]
    sinks_df.drop_duplicates(subset=["Sink name"], inplace=True)
    sinks = sinks_df["Sink name"].tolist()
    for i in range(len(sinks)):
        sinks[i] = sinks[i].replace("[", "").replace("]", "")

    for sink in sinks:
        print(f"Sink: {sink}")
        paths = find_pathways(df, sink)
        for path in paths:
            print(path[::-1], "\n")
        print("\n\n")