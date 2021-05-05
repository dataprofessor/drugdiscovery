"""
Docker Container: https://hub.docker.com/r/continuumio/anaconda3
RDKit Installation: https://www.rdkit.org/docs/Install.html
"""
import mols2grid
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

st.title("Filter FDA Approved Drugs by Molecular Weight with Streamlit")


@st.cache(allow_output_mutation=True)
def download_dataset():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv(
        "https://www.cureffi.org/wp-content/uploads/2013/10/drugs.txt", sep="\t"
    ).dropna()
    return df


def calc_mw(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the molecular weight"""
    mol = Chem.MolFromSmiles(smiles_string)
    return ExactMolWt(mol)


# Copy the dataset so any changes are not applied to the original cached version
df = download_dataset().copy()
df["mol_weight"] = df.apply(lambda x: calc_mw(x["smiles"]), axis=1)

weight_cutoff = st.slider(
    label="Show compounds that weigh below:",
    min_value=0,
    max_value=500,
    value=150,
    step=10,
)


df_result = df[df["mol_weight"] < weight_cutoff]
st.write(df_result)


raw_html = mols2grid.display(df_result, mapping={"smiles": "SMILES"})._repr_html_()
components.html(raw_html, width=900, height=900, scrolling=True)
