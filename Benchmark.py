import numpy as np
import os
import sys



### FUNCTIONX ###


# Va chercher les .fasta du benchmark
# Dictionnaire de query + nom
def RecupereBenchmarkQuery(repertoire):
    
    return


# Va dans Benchmark.list 
# output dictionnaire de dictionnaire nom + superfamille et fold associes    
def RecupereResultBenchmark():
    path = "./Benchmark/Benchmark.list"
    file = open(path,"r")
    return file



# creer dossier avec resultat + documentation
def main():
    #nomDossier=sys.argv[1]
    RecupereResultBenchmark()
    

if __name__ == "__main__":
    # execute only if run as a script
    main()

