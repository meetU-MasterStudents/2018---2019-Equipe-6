import numpy as np
import pandas as pd
import os
import sys



### FUNCTIONX ###

# Va chercher les .fasta du benchmark
# Dictionnaire de query + nom
#go into each folder to retrieve each query file.fasta

repertoire='test_dataset'

def RecupereBenchmarkQuery(repertoire):
    Listqueryfolder=[]
    Listfasta=[]
    for nom in os.listdir(repertoire):
        Listqueryfolder.append(nom)
        QueryFolder=repertoire+"/"+nom
        for queryfile in os.listdir(QueryFolder):
            if queryfile==nom+"."+"fasta": # engros je veux dire si mon extension de fichier est .fasta
                fastafile=open('test_dataset'+"/"+nom+"/"+"/"+nom+"."+"fasta", 'r')
                fastaread=pd.read_table(fastafile, engine="python")
                n,m=fastaread.shape
                Fil=fastaread.iloc[1:n-1,:m-1]
                fil=Fil.as_matrix()
                Listfasta.append(fil)
                fastafile.close()
                total = {}
            else:
                print("fasta file is missing")
    for x, y in zip(Listqueryfolder, fastafile):
            total[x] = y 
    return fastafile, Listqueryfolder, total

Data, Name, Dictionnaire_Fasta_Dataset =RecupereBenchmarkQuery(repertoire)


# Va dans Benchmark.list 
# output dictionnaire de dictionnaire nom + superfamille et fold associes    

def RecupereResultBenchmark():
    
    return



# creer dossier avec resultat + documentation
def main():
    #nomDossier=sys.argv[1]
    


if __name__ == "__main__":
    # execute only if run as a script
    main()

