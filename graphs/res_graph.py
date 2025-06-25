import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from pubchempy import get_compounds
from networkx.drawing.nx_agraph import graphviz_layout


def load_res(file_path):
    """
    Load the results from a CSV file into a pandas DataFrame.
    """
    return pd.read_csv(file_path)

def add_names(df):
    """
    Add 'Substrate Name' and 'Product Name' columns to the DataFrame based on InChI IDs
    Queries pubchem to get names.
    """

    def get_name(inchi):
        if inchi == 'InChI=1S':
            return "H+"
        compounds = get_compounds(inchi, 'inchi')
        if compounds:
            if compounds[0].synonyms:
                return compounds[0].synonyms[0]
            elif compounds[0].iupac_name:
                return compounds[0].iupac_name
            else:
                return "Unknown"
        return "Unknown"

    df['Substrate Name'] = df['Substrate InChI'].apply(get_name)
    df['Product Name'] = df['Product InChI'].apply(get_name)
    
    return df



if __name__ == "__main__":

    #file_path = "/work/igem/iGEM/graphs/res_cases/ota_MET_scope.csv"
    #file_path = "/work/igem/iGEM/graphs/res_cases/ota_results.csv"
    file_path = "/work/igem/iGEM/graphs/res_cases/otc_MET_scope.csv"
    df = load_res(file_path)
    print("Initial df\n", df)
    df = add_names(df)
    print("df with names\n", df)
    
    print("Substrates\n", df['Substrate Name'].unique().tolist())
    print("Products\n", df["Product Name"].unique().tolist())


    G = nx.from_pandas_edgelist(
                df=df, 
                source='Substrate Name', 
                target='Product Name', 
                edge_attr=True, 
                create_using=nx.DiGraph()
            )
    
    # Use Graphviz 'dot' layout for hierarchical layering
    pos = graphviz_layout(G, prog='dot')


    # Draw the graph
    #plt.figure(figsize=(16, 10))
    plt.figure(figsize=(32, 20))
    """
    nx.draw(
        G, 
        pos,
        with_labels=True, 
        node_color='lightblue', node_size=500, font_size=10,
        edge_color='gray', 
        arrows=True, arrowsize=20
        )
    """
    # Prepare edge colors and widths based on "In Sink" attribute
    edge_colors = []
    edge_widths = []
    for u, v, d in G.edges(data=True):
        if d.get("In Sink", 0) == 1:
            edge_colors.append("green")
            edge_widths.append(4)
        else:
            edge_colors.append("gray")
            edge_widths.append(1)

    nx.draw(
        G, 
        pos,
        with_labels=True, 
        node_color='lightblue', node_size=500, font_size=10,
        edge_color=edge_colors,
        width=edge_widths,
        arrows=True, arrowsize=20
    )
    
    # Draw EC number as edge labels
    edge_labels = nx.get_edge_attributes(G, 'EC number')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='red', font_size=8)


    plt.title("Directed Graph from Results")
    #plt.show()
    #plt.savefig("/work/igem/iGEM/results/ota_scope.png", dpi=400)
    #plt.savefig("/work/igem/iGEM/results/ota_results.png", dpi=400)
    plt.savefig("/work/igem/iGEM/results/otc_scope.png", dpi=400)