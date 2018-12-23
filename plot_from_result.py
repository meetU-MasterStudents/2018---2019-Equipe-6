import numpy as np

import os, shutil
import warnings
import sys, getopt
import matplotlib.pyplot as plt


Benchmark_Fold={"UBQ":"protg", "DEP":"myb_DNA-binding", "igvar-h":"Desulfoferrodox",
                "Rib_hydrolayse":"TIR", "Lum_binding":"EFTU_C", "DnaJ":"ATP-synt_DE_N",
                "Cohesin":"Fimbrial", "PAC":"GAF", "Agglutinin":"intb", "hemery":"SRP54", "LRR":"Recep_L_domain",
                "SSB":"TIMP", "TBCA": "BAG", "RRF":"RRF"}
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


    #Area under curve
    Aire=np.trapz(acc,range(threshold), dx=1)
    Aire=round(Aire/threshold,2)
    #plot    
    fig=plt.figure()    
    #plt.plot(range(threshold), acc, 'ro')
    lw = 2
    plt.plot(range(threshold), acc, color='darkorange',
         lw=lw, label='Semi-ROC curve area '+str(Aire),linewidth=3)
    y = np.linspace(0, 1, threshold, endpoint=False)
    plt.plot(range(threshold),y, linestyle='--')
    plt.ylim((-0.01,1.01))
    plt.xlim((-1,405))
    plt.xlabel("Threshold")
    plt.ylabel("Number of positives")
    plt.title(name)
    plt.legend(loc="lower right")
    fig.savefig(name+".png")
    plt.show()
    
    return acc

def Accuracy_according_to_Fold_and_SF(List_Benchmark_Fold,  List_Benchmark_SF, Results, threshold):
        #stocke les accuracy
    accSF=[0]*threshold
    accF=[0]*threshold
    
    
    for query in Results.keys():
        r=Results.get(query)
        rank_query=sorted(r, key=r.__getitem__) #Range les resultats du plus petit au plus grand
        rank_query.reverse() #list type
        if threshold>len(rank_query):
            threshold=len(rank_query)
            print("Pas assez de resultats possibles pour afficher le threshold demande dans le cas de "+query)
        accuracy_SF = 0  
        accuracy_Fold = 0
        for i in range(threshold):     
            f=rank_query[i]
            #Fold (accorder un poids)
            if (f==List_Benchmark_Fold.get(query)) and (f!= None) and (List_Benchmark_Fold.get(query)!= None):
                accuracy_Fold+=1
            #SF (accorder un poids)
            if (f==List_Benchmark_SF.get(query)) and (f!= None) and (List_Benchmark_SF.get(query)!= None):
                accuracy_SF+=1  
            accSF[i]=accSF[i]+accuracy_SF
            accF[i]=accF[i]+accuracy_Fold

    
    nb_queries_Fold=len(List_Benchmark_Fold.keys()) 
    nb_queries_SF=len(List_Benchmark_SF.keys())   
    for k in range(len(accSF)):
        accSF[k]=accSF[k]/nb_queries_SF
    for k in range(len(accF)):       
        accF[k]=accF[k]/nb_queries_Fold
        

    #plot     
    fig=plt.figure()  
    plt.subplot(1, 2, 1)
   
    #Area under curve
    AireSF=np.trapz(accSF,range(threshold), dx=1)
    AireSF=round(AireSF/threshold,2)
         
    lw = 2
    plt.plot(range(threshold), accSF, color='darkorange',
         lw=lw, label='Semi-ROC curve area '+str(AireSF),linewidth=3)
    y = np.linspace(0, 1, threshold, endpoint=False)
    plt.plot(range(threshold),y, linestyle='--')
    plt.ylim((-0.01,1.01))
    plt.xlim((-1,405))
    plt.xlabel("Threshold")
    plt.ylabel("Number of positives")
    plt.title("Super-Family", fontsize=15)
    plt.legend(loc="lower right")

    plt.subplot(1, 2, 2)
    #Area under curve
    AireF=np.trapz(accF,range(threshold), dx=1)
    AireF=round(AireF/threshold,2)
         
    lw = 2
    plt.plot(range(threshold), accF, color='darkorange',
         lw=lw, label='Semi-ROC curve area '+str(AireF),linewidth=3)
    y = np.linspace(0, 1, threshold, endpoint=False)
    plt.plot(range(threshold),y, linestyle='--')
    plt.ylim((-0.01,1.01))
    plt.xlim((-1,405))
    plt.xlabel("Threshold")
    plt.ylabel("Number of positives")
    plt.title("Fold", fontsize=15)
    plt.legend(loc="lower right")
    
    
    plt.tight_layout()
    plt.show()   
    fig.savefig('Fold_and_SuperFamily.jpg')
    
    return

import seaborn as sns
def plot_score(List_Benchmark_Fold,  List_Benchmark_SF, Results):
    list_query=Results.keys()
    nombre_abscisse=len(list_query)

  
    x=range(nombre_abscisse)
    Top_one=[]
    for query in list_query:
        results_for_the_query=Results.get(query)
        for res in results_for_the_query.keys():
            if res==List_Benchmark_Fold.get(query) or List_Benchmark_SF.get(query)==res: 
                Top_one.append(results_for_the_query[res])
        
    fig=plt.figure()
    plt.plot(x, Top_one, '.', color='r', marker='s')
    
    compteur=0

    for query in list_query:
        results_for_the_query=Results.get(query)
        not_matched=[]
        for res in results_for_the_query.keys():
            not_matched.append(results_for_the_query[res])
        c = [0.0, float(np.std(not_matched))/float(80), float(80-np.std(not_matched))/float(80)] #R,G,B 
        print(np.std(not_matched))
        plt.scatter([compteur]*len(not_matched), not_matched, alpha=0.6, marker='.', color=c)
        compteur=compteur+1
    plt.xticks(x, list_query, rotation='vertical')
    plt.show()
    fig.savefig("details_on_each_queries")


r1=np.load("check6.npy").item()
Affichage_Accuracy(Benchmark_Fold, Benchmark_SF, r1, 405, "control check")
Accuracy_according_to_Fold_and_SF(Benchmark_Fold, Benchmark_SF, r1, 405)
plot_score(Benchmark_Fold, Benchmark_SF, r1)