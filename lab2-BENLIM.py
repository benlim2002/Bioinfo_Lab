import requests
import pandas as pd
import networkx as nx
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(layout="wide")

# BioGRID's Functions
def retrieve_ppi_biogrid(target_protein):
    biogrid_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": "1000170679ebd9bab92680ec1bd55699",
        "format": "json",
        "searchNames": True,
        "geneList": target_protein,
        "organism": 9606,
        "searchbiogridids": True,
        "includeInteractors": True
    }
    response = requests.get(biogrid_url, params=params)
    return response.json()

# Process data into dataframe then create NetworkX Graph for visual
def process_biogrid_data(data):
    dataframe = pd.DataFrame.from_dict(data, orient='index')
    dataframe['OFFICIAL_SYMBOL_A'] = dataframe['OFFICIAL_SYMBOL_A'].str.upper()
    dataframe['OFFICIAL_SYMBOL_B'] = dataframe['OFFICIAL_SYMBOL_B'].str.upper()
    network_graph = nx.from_pandas_edgelist(dataframe, "OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B")
    return dataframe, network_graph



# STRING's Functions
def retrieve_ppi_string(target_protein):
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": 9606
    }
    response = requests.get(string_url, params=params)
    return response.json()

# Process data into dataframe then create NetworkX graph for visual
def process_string_data(data):
    dataframe = pd.json_normalize(data)
    network_graph = nx.from_pandas_edgelist(dataframe, "preferredName_A", "preferredName_B")
    return dataframe, network_graph



# Function to calculate and display centrality metrics
def get_centralities(network_graph):
    degree_cent = nx.degree_centrality(network_graph)
    betweenness_cent = nx.betweenness_centrality(network_graph)
    closeness_cent = nx.closeness_centrality(network_graph)
    eigenvector_cent = nx.eigenvector_centrality(network_graph, max_iter=100000000) #set max_iterations bcs some protein wont display centrality if iteration is too mcuh
    pagerank_cent = nx.pagerank(network_graph, max_iter=1000000000)
    
    def get_top25(centrality_data):
        return dict(sorted(centrality_data.items(), key=lambda x: x[1], reverse=True)[:25])#func to get top 25 only, since ther will be too many
    
    return {
        "Degree Centrality": get_top25(degree_cent),
        "Betweenness Centrality": get_top25(betweenness_cent),
        "Closeness Centrality": get_top25(closeness_cent),
        "Eigenvector Centrality": get_top25(eigenvector_cent),
        "PageRank Centrality": get_top25(pagerank_cent)
    }
    


# Function for network graph and struct    
# Function to visualize the network graph
def generate_network(network_graph, centrality_data):
    fig, ax = plt.subplots(figsize=(8, 6))
    slayout = nx.spring_layout(network_graph, seed=42)

    high_centrality_nodes = []
    for measure, top25 in centrality_data.items():
        high_centrality_nodes.extend(top25.keys())
        
    high_centrality_nodes = list(set(high_centrality_nodes))

    nx.draw(network_graph, slayout, with_labels=True, node_color='skyblue', node_size=1000, edge_color="gray", font_size=8, ax=ax)
    nx.draw_networkx_nodes(network_graph, slayout, nodelist=high_centrality_nodes, node_size=800, node_color='pink')
    
    return fig


# Function to describe the network structure
def describe_network(network_graph):
    num_edges = len(network_graph.edges)
    num_nodes = len(network_graph.nodes)
    description = f"Number of interactions in network: {num_edges} \n\nNumber of proteins in network: {num_nodes}"
    return description


st.title('Lab 2 - Ben Lim Choong Chuen')
st.subheader("Retrieve Human Protein-Protein Interactions from STRING or BioGRID")

# Input 
box, sel = st.columns([5,1])

with box:
    target_protein = st.text_input('Enter Protein Symbol')

# Selection buttons for BioGRID and STRING
with sel:
    if 'selected_source' not in st.session_state:
        st.session_state.selected_source = None

    if st.button("BioGRID"):
        st.session_state.selected_source = "BioGRID"
    
    if st.button("STRING"):
        st.session_state.selected_source = "STRING"

if st.button("Retrieve"):
    network_data = None
    if st.session_state.selected_source == "BioGRID":
            network_data = retrieve_ppi_biogrid(target_protein)
            dataframe, network_graph = process_biogrid_data(network_data)
            
    elif st.session_state.selected_source == "STRING":
            network_data = retrieve_ppi_string(target_protein)
            dataframe, network_graph = process_string_data(network_data)

    else:
        st.warning("Please select either BioGRID or STRING before retrieving data.")

    # display data
    if network_data:

        col1, col2 = st.columns([1, 1])  

        with col1:
            target_protein_text = f"<span style='color:orange;'>{target_protein}</span>"
            st.markdown(f"### {st.session_state.selected_source} PPI data for protein: {target_protein_text}", unsafe_allow_html=True)
            st.write(f"Network data size: {dataframe.shape[0]} rows and {dataframe.shape[1]} columns")
            st.dataframe(dataframe) 
            
            st.write("### Description of the network")
            st.markdown(describe_network(network_graph))
            st.write("***The top 25 are in pink***")

            centrality_data = get_centralities(network_graph)  # Get centrality data
            fig = generate_network(network_graph, centrality_data)  # Pass both arguments
            st.pyplot(fig)

        with col2:

            st.subheader('Centrality Measures')
            st.write("*Only top 25 for each centrality*")
            for measure, data in centrality_data.items():
                st.write(f"**{measure}:**")
                formatted_data = ", ".join([f"({node}, {value:.4f})" for node, value in data.items()])
                st.markdown(f"<p style='color:orange;'>{formatted_data}</p>", unsafe_allow_html=True)
                
    else:
        st.warning("No data found for the specified protein.")
