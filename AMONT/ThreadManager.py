from multiprocessing import Process,cpu_count,Manager
import pandas as pd
import os

from SW_vectors import *
from Printing import *

#Move them to the Benchmark.py
mu = 0.03316
sigma = 0.06867
gapPenalty = 12
misMatchPenalty = 3
gapExtension = 1
queryProfilesPath = "QueryResults" #Constant

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

def dotProduct(profile1, profile2, mu, sigma):
    matrix = np.array(calculateScore(profile1,profile2))
    matrix = normalise_SW(matrix,mu,sigma)
    return matrix


from scipy.stats import pearsonr
def pearsonCorrelationCoefficient(profile1, profile2, mu, sigma):
    n1,m1=np.shape(profile1)
    n2,m2=np.shape(profile2)
    Profile1=np.transpose(profile1)
    Profile2=np.transpose(profile2)
    MatrixScore=np.zeros((m1,m2))
    for i in range(m1):
        for j in range(m2):
            A=Profile1[i]
            B=Profile2[j]
            MatrixScore[i][j] = pearsonr(A,B)[0]
    MatrixScore = normalise_SW(MatrixScore,mu,sigma)
    return MatrixScore

def spearmannCorrelationCoefficient(profile1, profile2, mu, sigma):
    n1,m1=np.shape(profile1)
    n2,m2=np.shape(profile2)
    Profile1=np.transpose(profile1)
    Profile2=np.transpose(profile2)
    MatrixScore=np.zeros((m1,m2))
    for i in range(m1):
        for j in range(m2):
            A=list(Profile1[i])
            Abar=np.mean(A)
            B=list(Profile2[j])
            Bbar=np.mean(B)
            num=0
            denom=0
            Aquadra=0
            Bquadra=0
            for k in range(n1):
                num=num+(k-Abar)*(k-Bbar)
                Aquadra=Aquadra+(k-Abar)**2
                Bquadra=Bquadra+(k-Bbar)**2
            denom=np.sqrt(Aquadra*Bquadra)
            MatrixScore[i][j] = abs(num/denom)
    MatrixScore = normalise_SW(MatrixScore,mu,sigma)
    return MatrixScore

from scipy.spatial import distance
def euclideanDistance(profile1, profile2, mu, sigma):
    n1,m1=np.shape(profile1)
    n2,m2=np.shape(profile2)
    Profile1=np.transpose(profile1)
    Profile2=np.transpose(profile2)
    MatrixScore=np.zeros((m1,m2))
    for i in range(m1):
        for j in range(m2):
            MatrixScore[i][j] = distance.euclidean(Profile1[i], Profile2[j])
    MatrixScore = normalise_SW(MatrixScore,mu,sigma)
    return MatrixScore

def SequenceAlignment(path,seq1Name,seq1Cont,seq2Name,seq2Cont):
    inputFile = path+seq1Name+'_'+seq2Name
    outputFile = inputFile + '_aligned.fasta'
    with open(inputFile, 'w') as fHandler:
        fHandler.write('>'+seq1Name+'\n')
        fHandler.write(seq1Cont+'\n')
        fHandler.write('>'+seq2Name+'\n')
        fHandler.write(seq2Cont)
    os.system('muscle -in ' + inputFile + ' -out ' + outputFile)
    with open(outputFile, 'r') as fHandler:
        lines = fHandler.readlines()
    alignedSeqs = []
    for line in lines:
        if '>' in line:
            alignedSeqs.append('') 
            continue
        alignedSeqs[-1] += line.replace('\n','')
    return alignedSeqs[0],alignedSeqs[1]

def ProfileProcessor(query,seqHomstrad,profHomstrad,evalue,database,qProf,printOut,comparison,useSS,applyCorrl,applyW,remoteDB,return_dict):
    print('Process id: {0}'.format(os.getpid()))
    if(qProf):
        command = './Profile -p -q ' + query[2] + ' -e ' + evalue + ' -d ' + database
        if(applyW):
            command += ' -w '
        if(remoteDB):
            command += ' -r '
        #command += ' -m'
        os.system(command)
    
    if(comparison):
        queryPath = queryProfilesPath + '/' + query[0] + '/'
        queryFilePath = queryPath + query[0] + '_PSSMProfile'
        profQuery = GetQueryProfile(queryFilePath)
    
        return_dict[query[0]] = {}
        list_results = []
        for i in range(len(profHomstrad)):
            print('Profile comparison: ',query[0], ' VS ',profHomstrad[i][0])
            if applyCorrl:
                matrix=pearsonCorrelationCoefficient(profQuery, profHomstrad[i][1],mu,sigma)
            else:
                matrix=dotProduct(profQuery, profHomstrad[i][1],mu,sigma)
                
            ### AFFINE
            #scoreMat = forward_SW(matrix, gapPenalty, misMatchPenalty)
            scoreMat = local_alignment_affine_gap_penalty(matrix, gapPenalty, gapExtension)
            traceback = traceback_SW(scoreMat)
            return_dict[query[0]][profHomstrad[i][0]] = traceback[0]

            if(printOut):
                seqAln1,seqAln2 = SequenceAlignment(queryPath,query[0],query[1],profHomstrad[i][0],seqHomstrad[profHomstrad[i][0]])
                list_results.append(Result(query=query[0], name=profHomstrad[i][0],
                                            score = traceback[0], qseq=seqAln1, tseq=seqAln2,
                                            gaps = traceback[1], qbegin = traceback[2]+1, qend = traceback[4]+1, tbegin = traceback[3]+1, 
                                            tend = traceback[5]+1, qal = traceback[6], tal = traceback[7],
                                            norm_score=0.5,qcov=0.8,identity=0.05,ssscore=6.30,allength=len(seqAln1),corr=50))
        if(printOut):
            list_results_sorted = sorted(list_results, key=lambda x:x.score, reverse=True)
            print_all(query[0], len(query[1]), list_results_sorted, 'Output_'+query[0]+'.foldrec')

def MultiThreadQuery(queryList,homstradList,profilesHomstrad,evalue,database,qProf,printOut,comparison,useSS,applyCorrl,applyW,remoteDB):
    nQueries = len(queryList)
    manager = Manager()
    return_dict = manager.dict()
    procs = []
    for i in range(nQueries):
        proc = Process(target=ProfileProcessor, args=(queryList[i],homstradList,profilesHomstrad,evalue,database,qProf,printOut,comparison,useSS,applyCorrl,applyW,remoteDB,return_dict))
        procs.append(proc)
        proc.start()
 
    for proc in procs:
        proc.join()
    print('Benchmark score calculation for the given parameters is finished!')
    return return_dict


def MultiQuery(queryList,homstradList,profilesHomstrad,evalue,database,qProf,printOut,comparison,useSS,applyCorrl,applyW,remoteDB):
    nQueries = len(queryList)
    return_dict = dict()
    for i in range(nQueries):
        ProfileProcessor(queryList[i],homstradList,profilesHomstrad,evalue,database,qProf,printOut,comparison,useSS,applyCorrl,applyW,remoteDB,return_dict)
        #np.save(queryList[i][0],return_dict)
    print('Benchmark score calculation for the given parameters is finished!')
    return return_dict
