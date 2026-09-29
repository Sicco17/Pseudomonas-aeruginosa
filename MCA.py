'''

Code for the execution of the Multiple Correspondence Anlaysis
on the excel file containing all MIC values and ST labels for each isolate
and generation of 3 graphs:
    one for High risk clones
    one for emerging rick clones
    one for highly represented STs


'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  #dataframe generation
import prince  #execute MCA
import seaborn as sns #produce the graph

'''

Transform the different MIC values to ordinal categorical values based on the rank
This is done to include the left (< x) and right (> x) censored values 

'''
def mic_to_ordinal_str(val):
    if pd.isna(val):          #avoid errors for missing values
        return np.nan
    val_str = str(val).strip().replace(",", ".")

    if val_str == "" or val_str.lower() in ["nan", "null", "-", ""]:
        return np.nan

    if val_str.startswith(">"):    #replace right censored data with its rank
        try:
            num_part = float(val_str.replace(">", ""))
            return f"Rango_{int(np.log2(num_part)) + 1}(>)"
        except ValueError:
            return val_str
    elif val_str.startswith("<") or val_str.startswith("≤"):   #replace left censored data with its rank
        try:
            num_part = float(
                val_str.replace("<", "").replace("≤", "").replace('"', "")
            )
            return f"Rango_{int(np.log2(num_part))}(≤)"
        except ValueError:
            return val_str
    else:
        try:
            num_part = float(val_str.replace('"', ""))
            if num_part > 0:
                return f"Rango_{int(np.floor(np.log2(num_part)))}"
            return "Rango_0"
        except ValueError:
            return val_str

'''

INITIAL LINES: read the excel file

'''

file_excel = "onlyMIC.xlsx"  #name of the file containing MIC values
df_raw = pd.read_excel(file_excel)
df_raw.columns = df_raw.columns.str.strip()

'''

extract columns

'''

colonne_intervallo = list(df_raw.columns[2:20])
colonne_antibiotici = [
    c
    for c in colonne_intervallo
    if "Unnamed" not in str(c) and str(c).strip() != ""
]

'''

Bild the corresponding dataframe and apply the ordinal transformation 

'''

df_processed = pd.DataFrame(index=df_raw.index)
for col in colonne_antibiotici:
    df_processed[col] = df_raw[col].apply(mic_to_ordinal_str)

df_processed["ST_Metadata"] = df_raw.iloc[:, 1].astype(str).str.strip()
df_processed.index = df_raw.iloc[:, 0].astype(str).str.strip()

'''

Dataframe processing and metadata cleaning

'''

df_mca_pre = df_processed.dropna(subset=colonne_antibiotici, how="all")
righe_da_escludere = ["isolate id", "grand total", "row labels", "total", "nan"]
df_mca_pre = df_mca_pre[
    ~df_mca_pre.index.str.lower().str.strip().isin(righe_da_escludere)
]

st_labels = df_mca_pre["ST_Metadata"]
df_mca = df_mca_pre.drop(columns=["ST_Metadata"])

for col in df_mca.columns:
    if df_mca[col].isnull().all():
        df_mca = df_mca.drop(columns=[col]) 
        continue
    if df_mca[col].isnull().any():
        df_mca[col] = df_mca[col].fillna(df_mca[col].mode()[0])
    df_mca[col] = df_mca[col].astype("category")


print(
    f"Dataset validated! Samples: {df_mca.shape[0]}, Antibiotics: {df_mca.shape[1]}" 
)

#shape check
if df_mca.shape[0] == 0 or df_mca.shape[1] == 0:
    raise ValueError(
        "Error: The dataset is empty. Check your data alignment codes."
    )

'''

MCA execution and coordinates recovery for graph building

'''

mca = prince.MCA(n_components=2, random_state=42) #random_state as fixed set of random numbers for reproducibility
mca = mca.fit(df_mca)
coordinate_ceppi = mca.row_coordinates(df_mca)

# dimensions calculated based on prince parameter of native inertia

var_dim1 = float(mca.explained_inertia_[0]) * 100
var_dim2 = float(mca.explained_inertia_[1]) * 100

df_plot = pd.DataFrame(
    {
        "Dim1": coordinate_ceppi.iloc[:, 0],
        "Dim2": coordinate_ceppi.iloc[:, 1],
        "ST": st_labels.values,
    },
    index=coordinate_ceppi.index,
)

'''

Divide the different STs in 3 groups

'''

gruppi_st = [
    ["235", "111", "175", "654", "308", "357", "244", "233"], #high risk clones
    ["309", "773"],   #emerging high risk clones
    ["253", "2592"],  #remaining most represented STs 
]


'''

graph production based on the different defined groups

'''

sns.set_theme(style="whitegrid") # graph backgound

for idx, st_attivi in enumerate(gruppi_st, start=1): #one group at a time 
    plt.figure(figsize=(14, 10))
    df_attivi = df_plot[df_plot["ST"].isin(st_attivi)].copy()
    df_esclusi = df_plot[~df_plot["ST"].isin(st_attivi)].copy()
    if not df_esclusi.empty: #groups not considered represented as white circles
        sns.scatterplot(
            data=df_esclusi,
            x="Dim1",
            y="Dim2",
            color="white",
            edgecolor="black",
            linewidth=1.2,
            s=130,
            alpha=0.6,
            zorder=3,
            label="Other STs (Not considered)",
        )
        
    if not df_attivi.empty:   #considered groups are differentiated based on colour palette
        df_attivi["ST_Label"] = "ST " + df_attivi["ST"]
        ordine_legenda_attivi = [f"ST {st}" for st in st_attivi]
        palette_scelta = "tab10" if len(st_attivi) <= 10 else "tab20"
        colori_generati = plt.cm.get_cmap(palette_scelta)(
            np.linspace(0, 1, len(st_attivi))
        )
        
        mappa_colori = {
            f"ST {st}": colori_generati[i] for i, st in enumerate(st_attivi)
        }
        
        sns.scatterplot(
            data=df_attivi,
            x="Dim1",
            y="Dim2",
            hue="ST_Label",
            hue_order=ordine_legenda_attivi,
            palette=mappa_colori,
            s=150,
            edgecolor="black",
            linewidth=0.8,
            alpha=0.9,
            zorder=5,
        )

    # Central geometric grid axis lines
    plt.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    plt.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

    #Labels
    plt.xlabel(f"Dimension 1 ({var_dim1:.1f}%)", weight="bold", fontsize=11)
    plt.ylabel(f"Dimension 2 ({var_dim2:.1f}%)", weight="bold", fontsize=11)

    #Graph Titles
    lista_st_titolo = ", ".join(st_attivi)
    plt.title(
        f"MCA Mapping - [GRAPH {idx}]\nHighlighted Lineages: {lista_st_titolo}",
        fontsize=13,
        pad=15,
        weight="bold",
    )

    #Legend
    plt.legend(
        title="Sequence Type (ST)",
        title_fontproperties={"weight": "bold"},
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=True,
        facecolor="white",
    )

    plt.axis("equal")
    plt.tight_layout()

# Show all three graphs at the same time in the same window
plt.show()
