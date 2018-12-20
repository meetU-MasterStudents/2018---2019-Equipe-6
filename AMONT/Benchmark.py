import numpy as np
import pandas as pd
import os, shutil
import warnings
import sys, getopt
import matplotlib.pyplot as plt
from ThreadManager import *

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
def Affichage_Accuracy(List_Benchmark_Fold,  List_Benchmark_SF, Results, threshold, name):
    
    #stocke les accuracy
    acc=[0]*threshold
    
    
    for query in Results.keys():
        print("query is " +query)
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

    #plot    
    plt.figure()    
    plt.plot(range(threshold), acc, 'ro')
    y = np.linspace(0, 1, threshold, endpoint=False)
    plt.plot(range(threshold),y, linestyle='-.')
    plt.xlabel("Threshold")
    plt.ylabel("Semi-ROC")
    plt.title(name)
    plt.show()
    plt.savefig(name+".png")
    
    return acc

def usage(scriptFile):
    print("Usage: " + scriptFile + " \n" + 
          "Options:\n" + 
          "\t-h,--help\t\tShow this help message\n" +
          "\t-q,--qpath\t\t<Query path>\t\tPath to folder of queries\n" +
          "\t-m,--hpath\t\t<HOMSTRAD path>\t\tPath to the folder of HOMSTRAD dataset\n" +
          "\t-g,--confile\t\t<Config file>\t\tPath of the configuration file\n" +
          "\t-e,--evalue\t\t<e-Value>\t\te-Value for PSI-Blast\n" +
          "\t-d,--database\t\t<database>\t\tDatabase for PSI-Blast\n" +
          "\t-p,--prochoms\t\tCreate profiles for the HOMSTRAD dataset\n" +
          "\t-j,--qprof\t\tCreate profiles for the given queries\n" +
          "\t-w,--wInProf\t\tWeighing sequences during profile creation\n" +
          "\t-r,--mltproc\t\tRun in multiprocessing mode\n" +
          "\t-c,--recomp\t\tRecompile packages\n" +
          "\t-o,--output\t\tCreate output alignment\n" +
          "\t-x,--compare\t\tPerform profile-profile comparison\n" + 
          "\t-s,--secstru\t\tUse 2nd structure in profile-profile comparison\n" +
          "\t-l,--correl\t\tApply correlation in profile-profile comparison\n" +
          "\t-f,--rmtdb\t\tUse remote database for PSIBlast process")

def main(argv):
    benchmarkPath = ''
    homstradPath = ''
    justQueryProfiles = False
    processHomstrad = False
    multiProcess = False
    applyWeights = False
    configuration = False
    configPath = ''
    recompile = False
    printOutput = False
    performComparison = False
    applyCorrelation = False
    useSecStruct = False
    remoteDB = False
    homstradProfilesPath = "HomstradResults" #Constant
    queryProfilesPath = "QueryResults"       #Constant
    evalue = "1e-4"
    database = "swissprot"
    try:
        opts, args = getopt.getopt(argv[1:],"hq:m:jpe:d:wrcg:oxslf",["help", "qpath=", "hpath="
                                                                 "qprof", "prochoms", "evalue=",
                                                                 "database=", "wInProf", "mltproc"
                                                                 "recomp", "confile=", "output"
                                                                 "compare","secstru","correl","rmtdb"])
    except getopt.GetoptError:
        usage(argv[0])
        sys.exit(-1)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(argv[0])
            sys.exit()
        elif opt in ("-q", "--qpath"):
            benchmarkPath = arg
        elif opt in ("-m", "--hpath"):
            homstradPath = arg
        elif opt in ("-j","--qprof"):
            justQueryProfiles = True
        elif opt in ("-p","--prochoms"):
            processHomstrad = True
        elif opt in ("-e","--evalue"):
            evalue = arg
        elif opt in ("-d","--database"):
            database = arg
        elif opt in ("-w","--wInProf"):
            applyWeights = True
        elif opt in ("-r","--mltproc"):
            multiProcess = True
        elif opt in ("-c","--recomp"):
            recompile = True
        elif opt in ("-g","--confile"):
            configuration = True
            configPath = arg
        elif opt in ("-o","--output"):
            printOutput = True
        elif opt in ("-x","--compare"):
            performComparison = True
        elif opt in ("-s","--secstru"):
            useSecStruct = True   
        elif opt in ("-l","--correl"):
            applyCorrelation = True
        elif opt in ("-f","--rmtdb"):
            remoteDB = True                

    if(configuration):
        with open(configPath, 'r') as fHandler:
            lines = fHandler.readlines()
            for line in lines:
                option = line.replace('\n','').split(' ')
                if(option[0] == 'qpath'):
                    benchmarkPath = option[1]
                elif(option[0] == 'hpath'):
                    homstradPath = option[1]
                elif(option[0] == "qprof"):
                    justQueryProfiles = True
                elif(option[0] == "prochoms"):
                    processHomstrad = True
                elif(option[0] == "evalue"):
                    evalue = option[1]
                elif(option[0] == "database"):
                    database = option[1]
                elif(option[0] == "wInProf"):
                    applyWeights = True
                elif(option[0] == "mltproc"):
                    multiProcess = True
                elif(option[0] == "recomp"):
                    recompile = True
                elif(option[0] == "output"):
                    printOutput = True
                elif(option[0] == "compare"):
                    performComparison = True
                elif(option[0] == "secstru"):
                    useSecStruct = True 
                elif(option[0] == "correl"):
                    applyCorrelation = True
                elif(option[0] == "rmtdb"):
                    remoteDB = True

    # compare --> secstru --> correl!!
    if performComparison == False:
        useSecStruct = False
    if useSecStruct == True:
        applyCorrelation = True 

    warnings.filterwarnings("ignore")
    if(recompile):
        os.system('sudo g++ -o Profile main.cpp Profile.cpp -std=c++11')
    if(justQueryProfiles):
        if os.path.exists(queryProfilesPath):
            shutil.rmtree(queryProfilesPath)
        os.mkdir(queryProfilesPath)

    homstradDict = {}
    if(printOutput):
        homstrad = ReadDatabase(homstradPath)
        for i in range(len(homstrad)):
            homstradDict[homstrad[i][0]] = homstrad[i][1]
    
    
    queries = ReadDatabase(benchmarkPath)
    

    if(processHomstrad):
        if os.path.exists(homstradProfilesPath):
            shutil.rmtree(homstradProfilesPath)
        os.mkdir(homstradProfilesPath)
        os.system('./Profile -t '+homstradPath)

    dataProfileHomstrad = GetHomstradProfiles(homstradProfilesPath)

    if(multiProcess):
        Scores = MultiThreadQuery(queries,homstradDict,dataProfileHomstrad,evalue,database,
                                  justQueryProfiles,printOutput,performComparison,useSecStruct,
                                  applyCorrelation,applyWeights,remoteDB)
    else:
        Scores = MultiQuery(queries,homstradDict,dataProfileHomstrad,evalue,database,
                            justQueryProfiles,printOutput,performComparison,useSecStruct,
                            applyCorrelation,applyWeights,remoteDB)
    
    np.save(evalue+database,Scores)
    
    if performComparison==True:
        title=str(evalue)+"_"+database
        if applyCorrelation == True:
            title=title+"_correlation"
        else:
            title=title+"_dotProduct"
        if useSecStruct==True:
            title=title+"_withSecondaryStruct"
        Affichage_Accuracy(Benchmark_Fold,  Benchmark_SF, Scores, 400, title)

if __name__ == "__main__":
   main(sys.argv)
