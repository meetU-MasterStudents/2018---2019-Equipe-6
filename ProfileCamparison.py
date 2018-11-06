import numpy as np
import os
from io import BytesIO
import pandas as pd
from random import seed
from random import choice
from datetime import datetime

repertoire = 'Results'

def RecupereFiles(repertoire):
    data=[]
    for nom in os.listdir(repertoire) :
        File=open(repertoire+"/"+nom,'r')
        file=pd.read_csv(File, sep="\t| " , header=None, engine="python") 
        n,m=file.shape
        Fil=file.iloc[:n-1,1:m-1]
        fil=Fil.as_matrix()
        data.append(fil)
        File.close()
    return data
    
Data=RecupereFiles(repertoire)
#print((Data[0]))

def calculateScore(profile1,profile2):
    n1,m1=np.shape(profile1)
    n2,m2=np.shape(profile2)
    Profile1=np.transpose(profile1)
    Profile2=np.transpose(profile2)
    MatrixScore=np.zeros((m1,m2))
    for i in range(m1):
        for j in range(m2):
            MatrixScore[i][j]=np.dot(Profile1[i],Profile2[j])
    return MatrixScore
            
def calculateMean_Deviation(ScoreMatrix):
    mat = np.matrix(ScoreMatrix)
    mean = mat.mean()   
    Std = mat.std()
    #print('Mu: ', mean,', Sigma: ',Std)
    return mean, Std

def Selection(n,data):
    N=len(data)
    seed(datetime.now())
    vec=np.arange(N)
    selection=[]
    for _ in range(n):
          selection.append(choice(vec))
    return selection

def RandomSampling(n,data):
    sampleData=[]
    RandomSelection=Selection(n,data)
    for i in RandomSelection:
        sampleData.append(data[i])
    return sampleData

def ScoringProfile(data):
    mu=[]
    sigma=[]
    for profile1 in data:
        for profile2 in data:
            Score=calculateScore(profile1,profile2)
            mu_i,sigma_j=calculateMean_Deviation(Score)
            mu.append(mu_i)
            sigma.append(sigma_j)
    return mu,sigma

for i in range(5):
    SampleData=RandomSampling(40,Data)
    Mu,Sigma=ScoringProfile(SampleData)
    print('Mu: ', np.mean(Mu),', Sigma: ',np.mean(Sigma))

