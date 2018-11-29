import numpy as np
import pandas as pd
import os
import sys



### FUNCTIONS ###

# Va chercher les .fasta du benchmark
# Dictionnaire de query + nom
#go into each folder to retrieve each query file.fasta

Chemin="/home/roselyne/Bureau/BIM-info/MEETU/Benchmark building/test_dataset"

def Tuple_of_query(Chemin_Repertoire):
    Liste_Of_Tuple=[]
    for nom in os.listdir(Chemin_Repertoire):
        temporary=[nom]
        QueryFolder=Chemin_Repertoire+"/"+nom
        adress=QueryFolder+"/"+nom+".fasta"
        fastaFile=open(adress,"r")
        lines=fastaFile.readlines()
        temporary.append(lines[1])
        temporary.append(adress)
        Liste_Of_Tuple.append(tuple(temporary))
    return Liste_Of_Tuple


# Results
Benchmark_Fold={"UBQ":"protg Fold", "DEP":"myb_DNA-binding", "igvar-h":"Desulfoferrodox",
                "Rib_hydrolayse":"TIR", "Lum_binding":"EFTU_C", "DnaJ":"ATP-synt_DE_N Fold",
                "Cohesin":"Fimbrial Fold", "PAC":"GAF", "Agglutinin":"intb", "hemery":"SRP54", "LRR":"Recep_L_domain",
                "SSB":"TIMP", "TBCA": "BAG"}
Benchmark_SF={"His_biosynth":"OMPdecase", "svmp":"Astacin", "Lipoprotein_4":"mofe", "FAD-oxidase_NC": "MurB_C",
              "ETF_alpha":"Arginosuc_synth", "histone": "Arch_histone", "GA":"B", "PCNA":"DNA_PPF"} 


#Scoring function
def score(List_Benchmark_Fold,  List_Benchmark_SF, Results, threshold):
    accuracy=0
    for i in range(threshold):
        for query in Results[i].keys():
            r=Results[i].get(query)
            #Fold (accorder un poids)
            f=r.get('Fold')
            if (f==List_Benchmark_Fold.get(query)) and (f!= None) and (List_Benchmark_Fold.get(query)!= None):
                accuracy+=1
            #SF (accorder un poids)
            sf=r.get('SF')
            if (sf==List_Benchmark_SF.get(query)) and (sf!= None) and (List_Benchmark_SF.get(query)!= None):
                accuracy+=1
    print("Pour les "+str(threshold)+" meilleurs scores on trouve une accuracy de "+str(accuracy))   
    return accuracy

#Plot
import matplotlib.pyplot as plt
def Affichage_Accuracy(List_Benchmark_Fold,  List_Benchmark_SF, Results, threshold):
    acc=[]
    accuracy=0
    for i in range(threshold):
        for query in Results[i].keys():
            r=Results[i].get(query)
            #Fold (accorder un poids)
            f=r.get('Fold')
            if (f==List_Benchmark_Fold.get(query)) and (f!= None) and (List_Benchmark_Fold.get(query)!= None):
                accuracy+=1
            #SF (accorder un poids)
            sf=r.get('SF')
            if (sf==List_Benchmark_SF.get(query)) and (sf!= None) and (List_Benchmark_SF.get(query)!= None):
                accuracy+=1  
        acc.append(accuracy)
    plt.plot(range(threshold), acc, 'ro')
    plt.xlabel("Seuil")
    plt.ylabel("Accuracy")
    plt.show()
    return
    



Benchmark_as_tuple=Tuple_of_query(Chemin)

#Score donne 1/3 
#Score 1 & 2 : 2/3
# tous : 3/3 
Small_Results=[{"SSB":{"Fold":"Taratata"}, "PAC":{"SF":"Toto"}, "Lipoprotein_4":{"SF":"mofe"}}, 
               {"SSB":{"Fold":"TIMP"}, "PAC":{"SF":"GAF"}, "Lipoprotein_4":{"SF":"histone"}},
               {"SSB":{"SF":"blabla"}, "PAC":{"Fold":"GAF"}, "Lipoprotein_4":{"Fold":"BAG"}}]


Affichage_Accuracy(Benchmark_Fold, Benchmark_SF, Small_Results, 3)


