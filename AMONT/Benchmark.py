import numpy as np
import pandas as pd
import os, shutil
import warnings
from ThreadManager import *

### FUNCTIONS ###

#Put HELP here
# Va chercher les .fasta du benchmark
# Dictionnaire de query + nom
#go into each folder to retrieve each query file.fasta

warnings.filterwarnings("ignore")
benchmarkPath="../2018---2019-partage-master_old/Data/test_dataset"
homstradPath="..//2018---2019-partage-master_old/Data/HOMSTRAD"
if os.path.exists('QueryResults'):
    shutil.rmtree('QueryResults')
os.mkdir('QueryResults')
homstradProfilesPath = "HomstradResults"

processHomstrad = False
 
def ReadDatabase(Chemin_Repertoire):
    Liste_Of_Tuple=[]
    for nom in os.listdir(Chemin_Repertoire):
        temporary=[nom]
        QueryFolder=Chemin_Repertoire+"/"+nom
        adress=QueryFolder+"/"+nom+".fasta"
        fastaFile=open(adress,"r")
        lines=fastaFile.readlines()
        temporary.append(lines[1].strip())
        temporary.append(adress)
        Liste_Of_Tuple.append(tuple(temporary))
    return Liste_Of_Tuple

def GetHomstradProfiles(repertoire):
    data = []
    for nom in os.listdir(repertoire):
        File=open(repertoire+"/"+nom,'r')
        file=pd.read_csv(File, sep="\t| " , header=None, engine="python") 
        n,m=file.shape
        Fil=file.iloc[:n-1,1:m-1]
        fil=Fil.as_matrix()
        data.append((nom.replace('_PSSMProfile',''),fil))
        File.close()
    return data

homstrad = ReadDatabase(homstradPath)
homstradDict = {}
for i in range(len(homstrad)):
    homstradDict[homstrad[i][0]] = homstrad[i][1]
queries = ReadDatabase(benchmarkPath)

if(processHomstrad):
    if os.path.exists(homstradProfilesPath):
        shutil.rmtree(homstradProfilesPath)
    os.mkdir(homstradProfilesPath)
    os.system('./Profile -t '+homstradPath)

dataProfileHomstrad = GetHomstradProfiles(homstradProfilesPath)

#parameter loop here
evalue = "1e-8"
database = "swissprot"
#Scores = MultiThreadQuery(queries,homstradDict,dataProfileHomstrad,evalue,database)
Scores = MultiQuery(queries,homstradDict,dataProfileHomstrad,evalue,database)

# Results
Benchmark_Fold={"UBQ":"protg", "DEP":"myb_DNA-binding", "igvar-h":"Desulfoferrodox",
                "Rib_hydrolayse":"TIR", "Lum_binding":"EFTU_C", "DnaJ":"ATP-synt_DE_N",
                "Cohesin":"Fimbrial", "PAC":"GAF", "Agglutinin":"intb", "hemery":"SRP54", "LRR":"Recep_L_domain",
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
        accuracy /= len(Results[i].keys()) 
        acc.append(accuracy)
    plt.plot(range(threshold), acc, 'ro')
    plt.xlabel("Seuil")
    plt.ylabel("Accuracy")
    plt.show()
    return