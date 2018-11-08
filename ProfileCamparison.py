import numpy as np
import os
from io import BytesIO
import pandas as pd
from random import seed
from random import choice
import random
from datetime import datetime


repertoire = 'Results'

#Fonction qui permet de lire les résultats des profils
def RecupereFiles(repertoire):
    data=[]
    reference=[]
    for nom in os.listdir(repertoire) :
        reference.append(nom)
        File=open(repertoire+"/"+nom,'r')
        file=pd.read_csv(File, sep="\t| " , header=None, engine="python") 
        n,m=file.shape
        Fil=file.iloc[:n-1,1:m-1]
        fil=Fil.as_matrix()
        data.append(fil)
        File.close()
    total = {} 
    for x, y in zip(reference, data):
            total[x] = y
    return data, reference, total
    
Data, Name, Dictionnaire_Profile_Homstrad =RecupereFiles(repertoire)
#print((Data[0]))

#Produit Matriciel qui renvoir une matrice de score
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

def RandomSampling(n,dico):
    RandomDico = random.sample(dico.items(), k=n)
    RandomDico=dict(RandomDico)
    return RandomDico

def ScoringProfile(data):
    mu=[]
    sigma=[]
    for cle1, valeur1 in data.items():
        for cle2, valeur2 in data.items():
            if cle2!=cle1:
                Score=calculateScore(valeur1,valeur2)
                mu_i,sigma_j=calculateMean_Deviation(Score)
                mu.append(mu_i)
                sigma.append(sigma_j)
            #else:
                #print(cle1, cle2)
    return mu,sigma


#Definir le mu et le sigma
def Mu_And_Sigma(Nb_iter, N_Matrix_Per_Iter, dico):
    print("Mu and Sigma determination for",Nb_iter,"iteration(s) with ",N_Matrix_Per_Iter," profiles per iteration.")
    Result_Mu=[]
    Result_Sigma=[]
    Sample_following=[]
    for i in range(Nb_iter):
        SampleData=RandomSampling(N_Matrix_Per_Iter,dico)
        Sample_following.append(SampleData)
        Mu,Sigma=ScoringProfile(SampleData)
        Result_Mu.append(Mu)
        Result_Sigma.append(Sigma)
        print('Mu: ', round(np.mean(Mu),5),', Sigma: ',round(np.mean(Sigma),5))

    return Sample_following, Result_Mu, Result_Sigma 
   

test,_,_=Mu_And_Sigma(5, 3, Dictionnaire_Profile_Homstrad)
#Effectuer le dot product entre tous les profils et la query, et stocker les scores.

