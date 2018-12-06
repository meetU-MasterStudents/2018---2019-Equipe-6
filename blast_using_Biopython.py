#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 18:22:18 2018
"""

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import os



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


# Entree Adresse Query
# Sort fichier fasta avec la liste des hits sur seuil evalue demande
def FromQueryToBLAST(adress_query, adress_result, evalue):
    my_query = SeqIO.read(adress_query, format="fasta")
    result_handle = NCBIWWW.qblast("blastp","swissprot", my_query.seq)
    output = open(adress_result, "w")
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < evalue:
                output.write(">"+alignment.title.split(">")[0]+"\n")
                output.write(hsp.sbjct+"\n")
    result_handle.close()
    output.close()

#FromQueryToBLAST("Cohesin.fasta", "result_BLAST.fasta", 1e-5)


# Entree Adresse Query
# Sort fichier fasta avec la liste des hits sur seuil evalue demande
def FromQueryToPSIBLAST(adress_query, adress_result, evalue):
    my_query = SeqIO.read(adress_query, format="fasta")
    psi_blast = NCBIWWW.qblast("blastp", "swissprot", my_query.seq, service="psi")
    output = open(adress_result, "w")
    psiblast_record = NCBIXML.read(psi_blast)
    for alignment in psiblast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < evalue:
                output.write(">"+alignment.title.split(">")[0]+"\n")
                output.write(str(hsp.match)+"\n")
                #output.write(hsp.sbjct+"\n")
    psi_blast.close()
    output.close()
    
queries=ReadDatabase("/home/roselyne/Bureau/BIM-info/Test Biopython/test_dataset")

for i in range(len(queries)):
    FromQueryToPSIBLAST(queries[i][2], "./e5/"+queries[i][0]+".fasta", 1e-5)
    FromQueryToPSIBLAST(queries[i][2], "./e10/"+queries[i][0]+".fasta", 1e-10)
