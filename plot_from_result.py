import numpy as np

import os, shutil
import warnings
import sys, getopt
import matplotlib.pyplot as plt


Benchmark_Fold={"UBQ":"protg", "DEP":"myb_DNA-binding", "igvar-h":"Desulfoferrodox",
                "Rib_hydrolayse":"TIR", "Lum_binding":"EFTU_C", "DnaJ":"ATP-synt_DE_N",
                "Cohesin":"Fimbrial", "PAC":"GAF", "Agglutinin":"intb", "hemery":"SRP54", "LRR":"Recep_L_domain",
                "SSB":"TIMP", "TBCA": "BAG"}
Benchmark_SF={"His_biosynth":"OMPdecase", "svmp":"Astacin", "Lipoprotein_4":"mofe", "FAD-oxidase_NC": "MurB_C",
              "ETF_alpha":"Arginosuc_synth", "histone": "Arch_histone", "GA":"B", "PCNA":"DNA_PPF"} 

#Plot
def Affichage_Accuracy(List_Benchmark_Fold,  List_Benchmark_SF, Results, threshold, name):
    
    #stocke les accuracy
    acc=[0]*threshold
    
    
    for query in Results.keys():
        r=Results.get(query)
        rank_query=sorted(r, key=r.__getitem__) #Range les resultats du plus petit au plus grand
        rank_query.reverse() #list type
        if threshold>len(rank_query):
            threshold=len(rank_query)
            print("Pas assez de resultats possibles pour afficher le threshold demande dans le cas de "+query)
        print(List_Benchmark_Fold.get(query))
        print(List_Benchmark_SF.get(query))
        accuracy_inter = 0  
        for i in range(threshold):     
            f=rank_query[i]
            #Fold (accorder un poids)
            if (f==List_Benchmark_Fold.get(query)) and (f!= None) and (List_Benchmark_Fold.get(query)!= None):
                accuracy_inter+=1
            #SF (accorder un poids)
            if (f==List_Benchmark_SF.get(query)) and (f!= None) and (List_Benchmark_SF.get(query)!= None):
                accuracy_inter+=1  
            acc[i]=acc[i]+accuracy_inter

    
    nb_queries=len(Results.keys())        
    for k in range(len(acc)):
        acc[k]=acc[k]/nb_queries

    #plot    
    plt.figure()    
    plt.plot(range(threshold), acc, 'ro')
    y = np.linspace(0, 1, threshold, endpoint=False)
    plt.plot(range(threshold),y, linestyle='-.')
    plt.xlabel("Threshold")
    plt.ylabel("Semi-ROC")
    plt.title(name)
    plt.savefig(name+".png")
    plt.show()
    
    return acc


def plot_score(List_Benchmark_Fold,  List_Benchmark_SF, Results):
    list_query=Results.keys()
    nombre_abscisse=len(list_query)
    plt.figure()
    x=range(nombre_abscisse-1)
    Top_one=[]
    for query in list_query:
        results_for_the_query=Results.get(query)
        for res in results_for_the_query.keys():
            if res==List_Benchmark_Fold.get(query) or List_Benchmark_SF.get(query)==res: 
                Top_one.append(results_for_the_query[res])
        
    plt.plot(x, Top_one, '.', color='r', marker='s')
    
    compteur=0
    for query in list_query:
        results_for_the_query=Results.get(query)
        not_matched=[]
        for res in results_for_the_query.keys():
            not_matched.append(results_for_the_query[res])
        plt.scatter([compteur]*len(not_matched), not_matched, alpha=0.6, marker='.', color='b')
        compteur=compteur+1
    plt.xticks(x, list_query, rotation='vertical')
    plt.show()


r1=np.load("checkComparisonResults2.npy").item()
Affichage_Accuracy(Benchmark_Fold, Benchmark_SF, r1, 405, "control check")
plot_score(Benchmark_Fold, Benchmark_SF, r1)