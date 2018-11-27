import numpy as np
import pandas as pd
import os

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from SW_vectors import normalise_SW,forward_SW,traceback_SW,print_alignement
from Printing import *

repertoire = 'Results'
homstrad = '../HOMSTRAD/'

def GetResults(repertoire):
    homstradData = []
    for nom in os.listdir(repertoire):
        File=open(repertoire+"/"+nom,'r')
        file=pd.read_csv(File, sep="\t| " , header=None, engine="python") 
        n,m=file.shape
        Fil=file.iloc[:n-1,1:m-1]
        fil=Fil.as_matrix()
        if(nom == 'query.fasta_PSSMProfile'):
            queryData = ((nom.replace('_PSSMProfile',''),fil))
        else: 
            homstradData.append((nom.replace('_PSSMProfile',''),fil))
        File.close()
    return homstradData, queryData

def GetHomstrad(homstrad):
    recordsHomstrad = {}
    for famille in os.listdir(homstrad):
        for record in SeqIO.parse(homstrad+famille+'/'+famille+'.fasta',"fasta"):
            recordsHomstrad[record.id] = str(record.seq)
        for record in SeqIO.parse('query.fasta','fasta'):
            recordQuery = (record.id, str(record.seq))
    return recordsHomstrad,recordQuery

    
homstradData, queryData = GetResults(repertoire)
recordsHomstrad, recordQuery = GetHomstrad(homstrad)

#Produit Matriciel qui renvoir une matrice de score
def calculateScore(profile1,profile2):
    n1,m1=np.shape(profile1)
    n2,m2=np.shape(profile2)
    Profile1=np.transpose(profile1)
    Profile2=np.transpose(profile2)
    MatrixScore=np.zeros((m1,m2))
    for i in range(m1):
        for j in range(m2):
            MatrixScore[i][j] = np.dot(Profile1[i],Profile2[j]) #talk to irene
    return MatrixScore

mu = 0.03316
sigma = 0.06867
gapPenalty = -2
misMatchPenalty = -1

list_results = []

for i in range(2):
    matrix = np.array(calculateScore(queryData[1],homstradData[i][1]))
    matrix = normalise_SW(matrix,mu,sigma)
    score = forward_SW(matrix, gapPenalty, misMatchPenalty)
    traceback = traceback_SW(matrix)
    print('Processing: ',recordQuery[0], ' against ',homstradData[i][0])
    list_results.append(Result(query=recordQuery[0], name=homstradData[i][0],
                               score = traceback[0], qseq=recordQuery[1], tseq=recordsHomstrad[homstradData[i][0]],
                               gaps = traceback[1], qbegin = traceback[2], qend = traceback[3], tbegin = traceback[4], 
                               tend = traceback[5], qal = traceback[6], tal = traceback[7]))
print_all(recordQuery[0], len(recordQuery[1]), list_results, 'firstOut')
    