from multiprocessing import Process,cpu_count,Manager
import pandas as pd
import os

from SW_vectors import *
from Printing import *

#Move them to the Benchmark.py
mu = 0.03316
sigma = 0.06867
gapPenalty = -2
misMatchPenalty = -1

queryProfilesPath = "QueryResults"

def GetQueryProfile(filePath):
    File=open(filePath,'r')
    file=pd.read_csv(File, sep="\t| " , header=None, engine="python") 
    n,m=file.shape
    Fil=file.iloc[:n-1,1:m-1]
    fil=Fil.as_matrix()
    data = fil
    File.close()
    return data

#Produit Matriciel qui renvoie une matrice de score
def calculateScore(profile1,profile2):
    n1,m1=np.shape(profile1)
    n2,m2=np.shape(profile2)
    Profile1=np.transpose(profile1)
    Profile2=np.transpose(profile2)
    MatrixScore=np.zeros((m1,m2))
    for i in range(m1):
        for j in range(m2):
            MatrixScore[i][j] = np.dot(Profile1[i],Profile2[j])
    return MatrixScore

def Profile2Comparison(query,seqHomstrad,profHomstrad,evalue,database,return_dict):
    print('Process id: {0}'.format(os.getpid()))
    os.system('./Profile -p -q ' + query[2] + ' -e ' + evalue + ' -d ' + database)
    queryPath = queryProfilesPath + '/' + query[0] + '/' + query[0] + '_PSSMProfile'
    profQuery = GetQueryProfile(queryPath)
    
    return_dict[query[0]] = {}
    list_results = []
    for i in range(len(profHomstrad)):
        print('Processing: ',query[0], ' VS ',profHomstrad[i][0])
        matrix = np.array(calculateScore(profQuery,profHomstrad[i][1]))
        matrix = normalise_SW(matrix,mu,sigma)
        score = forward_SW(matrix, gapPenalty, misMatchPenalty)
        traceback = traceback_SW(matrix)
        return_dict[query[0]][profHomstrad[i][0]] = score

        
        list_results.append(Result(query=query[0], name=profHomstrad[i][0],
                               score = traceback[0], qseq=query[1], tseq=seqHomstrad[profHomstrad[i][0]],
                               gaps = traceback[1], qbegin = traceback[2], qend = traceback[3], tbegin = traceback[4], 
                               tend = traceback[5], qal = traceback[6], tal = traceback[7]))
    #print_all(query[0], len(query[1]), list_results, 'firstOut'+query[0])
        

    #Save data!!


def MultiThreadQuery(queryList,homstradList,profilesHomstrad,evalue,database):
    nQueries = len(queryList)
    manager = Manager()
    return_dict = manager.dict()
    procs = []
    for i in range(nQueries):
        proc = Process(target=Profile2Comparison, args=(queryList[i],homstradList,profilesHomstrad,evalue,database,return_dict))
        procs.append(proc)
        proc.start()
 
    for proc in procs:
        proc.join()
    print('Benchmark score calculation for the given parameters is finished!')
    return return_dict


def MultiQuery(queryList,homstradList,profilesHomstrad,evalue,database):
    nQueries = len(queryList)
    manager = Manager()
    return_dict = manager.dict()
    for i in range(nQueries):
        Profile2Comparison(queryList[i],homstradList,profilesHomstrad,evalue,database,return_dict)
    print('Benchmark score calculation for the given parameters is finished!')
    return return_dict