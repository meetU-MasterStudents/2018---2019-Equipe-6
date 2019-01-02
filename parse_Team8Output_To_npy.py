#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 11:11:01 2019

@author: roselyne
"""
import os
import numpy as np
def ReadDatabase(Chemin_Repertoire):
    DicoRes={}
    for nom in os.listdir(Chemin_Repertoire):
        queryname=nom.split("_")[1].split(".")[0]
        fastaFile=open(Chemin_Repertoire+"/"+nom,"r")
        dico_tempo={}
        for lines in fastaFile.readlines():
            to_list=lines.split()
            dico_tempo[to_list[1]]=to_list[2]
        DicoRes[queryname]=dico_tempo
    return DicoRes



chemin='/home/roselyne/Bureau/benchmark_equipe8/benchmark'
Results=ReadDatabase(chemin)
np.save('benchmark_aval.npy',Results)