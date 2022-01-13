import numpy as np
import pandas as pd
import math
from itertools import product

def affichage_resultat(dataframe):
    """
    Affiche les résultats sous forme de dataframe.
    Permet de filtrer les résultats en indiquant si on connait déjà le nombre d'un atome dans la molécule.
    """
    atome_ou_pas = input("Y-a-t-il un atome pour lequel sa quantité est connue dans la molécule? (o/n) :")
    if(atome_ou_pas=='o'):
        #Si on sait, par exemple, que l'ensemble de la molécule que l'on cherche ne possède que 14 atomes d'azotes
        #On ne va afficher que les lignes qui satisfont ce résultat
        nom_atome = input("Nom de l'atome (C, H, O, N, S, Fe, Al) : ")
        nb_atome = input(f"Combien y-a-t-il de {nom_atome} : ")
        if(nom_atome in dataframe.columns):
            #On vérifie quand même que l'atome est présent dans la dataframe
            dataframe[dataframe[nom_atome]== nb_atome]
        else:
            print("Cet atome n'est pas considéré.")
    elif(atome_ou_pas=='n'):
        #Sinon on affiche l'ensemble du dataframe, si aucun atome connu
        dataframe
        
    else:
        #L'utilisateur n'a ni tapé 'o' ni tapé 'n'
        print("Ce n'était pas une des options attendues.")

def trouver_molecule(masse_exp):
    """
    Trouve les molécules possibles pour la masse expérimentale donnée.
    Renvoie les molécules possibles dans un dataframe.
    """
    molecules = pd.read_excel('Pyoverdine Germain Final.xlsx', header=2, sheet_name="Molécules")
    molecules = molecules.drop(molecules.columns[16:], axis=1)
    acides_amines = molecules.iloc[:19]
    cyclisation = molecules.iloc[19:21]
    chromophores = molecules.iloc[21:27]
    side_chain = molecules.iloc[27:-4]
    metaux = molecules.iloc[-2:]
    resultats = pd.read_excel('Pyoverdine Germain Final.xlsx', header=3, sheet_name="Charges")
    resultats = resultats.drop([1,3]) # On enlève les lignes vides du tableau

    #Création d'un dictionnaire contenant la masse théorique de chacune des molécules
    masses_acides = {}
    masses_cycle = {}
    masses_chromo = {}
    masses_side_chain = {}
    masses_metaux = {}
    ranges_ = [range(2)] * 19
    for index, row in acides_amines.iterrows():
        masses_acides[row["ID"]]=row["Masse théo"]
    for index, row in cyclisation.iterrows():
        masses_cycle[row["ID"]]=row["Masse théo"]
    for index, row in chromophores.iterrows():
        masses_chromo[row["ID"]]=row["Masse théo"]
    for index, row in side_chain.iterrows():
        masses_side_chain[row["ID"]]=row["Masse théo"]
    for index, row in metaux.iterrows():
        masses_metaux[row["ID"]]=row["Masse théo"]

    candidates = []
    #On considère les candidats pour 2 mêmes molécules max
    ranges_familles = [range(3)]*5

    for nuplets in product(*ranges_):
        # On calcule la masse théorique des acides aminés
        masse_theo_acides = sum([nuplets[i] * list(masses_acides.values())[i] for i in range(len(nuplets))])
        
        # Puis, on regarde pour chaque chromophore, side_chain, métal et cycles si ça correspond
        # (les conditions de type "Si Ch3 alors pas de side" ne sont pas encore implémentées)
        for chromo in masses_chromo:
            for side in masses_side_chain:
                for met in masses_metaux:
                    for nuplet_famille in product(*ranges_familles):
                        masses_familles = [masses_cycle["Cyclisation -H2O"], masses_cycle["OH Terminal"], masses_chromo[chromo], masses_side_chain[side], masses_metaux[met]]
                        masse_theo_totale = masse_theo_acides + sum([nuplet_famille[i] * masses_familles[i] for i in range(len(nuplet_famille))])
                        if(math.isclose(masse_exp, masse_theo_totale, rel_tol=1e-05)):
                            # Si la masse théorique est à peu près égale à l'expérimentale,
                            # alors on récupère le nombre de chaque molécule présente.
                            tuple_cycles = (nuplet_famille[0], nuplet_famille[1])
                            tuple_chromo = (0,) * pd.Index(chromophores['ID']).get_loc(chromo) + (nuplet_famille[2],) + (0,) * (len(chromophores['ID']) - 1 - pd.Index(chromophores['ID']).get_loc(chromo))
                            tuple_side = (0,) * pd.Index(side_chain['ID']).get_loc(side) + (nuplet_famille[3],) + (0,) * (len(side_chain['ID']) - 1 - pd.Index(side_chain['ID']).get_loc(side))
                            tuple_metaux = (0,) * pd.Index(metaux['ID']).get_loc(met) + (nuplet_famille[4],) + (0,) * (len(metaux['ID']) - 1 - pd.Index(metaux['ID']).get_loc(met))
                            resultat = nuplets + tuple_cycles + tuple_chromo + tuple_side + tuple_metaux
                            #Si COH=1 et OH terminal 0 ou COH=0 et on vérifie pas pour OH terminal
                            if((nuplets[14]>=1 and tuple_cycles[1]==0) or nuplets[14]==0):
                                #Si Chr3 ou Chr4 >= 1 et pas de side chain ou alors Chr3 et Chr4 = 0
                                if(((tuple_chromo[1]>=1 or tuple_chromo[2]>=1) and sum(list(tuple_side))==0)) or (tuple_chromo[1]==0 and tuple_chromo[2]==0):
                                    #1 chromophore et 1 side chain au plus
                                    if(sum(list(tuple_chromo))<=1 and sum(list(tuple_side))<=1):
                                        candidates.append(resultat)

    # On crée un DataFrame pour rendre les résultats plus lisibles
    candidates_df = resultats.iloc[4:,4:]

    pd.options.display.max_columns = None
    pd.options.display.max_rows = None

    # Pour chaque combinaison possible, on va calculer le nombre de chaque atome présent
    for candidate in candidates:
        candidate2 = candidate + (0,)*(len(candidates_df.columns)-len(candidate))
        serie = pd.Series(candidate2, index=candidates_df.columns)
        atomes = {"C":0, "H":0, "N":0, "O":0, "S":0, "Fe":0, "Al":0}
        for index, value in serie.items():
            if(index in molecules["ID"].values):
                atomes["C"] += molecules.fillna(0)[molecules["ID"]==index]["NbrC"].item() * value
                if(index == "Fe(III) - 3 H+" or index == "Al(III) - 3 H+"):
                    atomes["H"] -= molecules.fillna(0)[molecules["ID"]==index]["NbrH"].item() * value
                else:
                    atomes["H"] += molecules.fillna(0)[molecules["ID"]==index]["NbrH"].item() * value
                atomes["N"] += molecules.fillna(0)[molecules["ID"]==index]["NbrN"].item() * value
                atomes["O"] += molecules.fillna(0)[molecules["ID"]==index]["NbrO"].item() * value
                atomes["S"] += molecules.fillna(0)[molecules["ID"]==index]["NbrS"].item() * value
                atomes["Fe"] += molecules.fillna(0)[molecules["ID"]==index]["NbrFe"].item() * value
                atomes["Al"] += molecules.fillna(0)[molecules["ID"]==index]["NbrAl"].item() * value
        serie = pd.Series(candidate + (0,) + tuple(atomes.values()), index=candidates_df.columns)
        candidates_df = candidates_df.append(serie,ignore_index=True)
    return(candidates_df)
