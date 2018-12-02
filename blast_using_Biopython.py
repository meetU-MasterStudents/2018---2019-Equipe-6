#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 18:22:18 2018
"""

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


# Entree Adresse Query
# Sort fichier fasta avec la liste des hits sur seuil evalue demande
def FromQueryToBLAST(adress_query, adress_result, evalue):
    my_query = SeqIO.read(adress_query, format="fasta")
    result_handle = NCBIWWW.qblast("blastp", "nr", my_query.seq)
    output = open(adress_result, "w")
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < evalue:
                output.write(">"+alignment.title.split(">")[0]+"\n")
                output.write(hsp.sbjct+"\n")
    result_handle.close()
    output.close()

#FromQueryToBLAST("Agglutinin.fasta", "result_BLAST.fasta", 0.04)


# Entree Adresse Query
# Sort fichier fasta avec la liste des hits sur seuil evalue demande
def FromQueryToPSIBLAST(adress_query, adress_result, evalue):
    my_query = SeqIO.read(adress_query, format="fasta")
    psi_blast = NCBIWWW.qblast("blastp", "refseq_protein", my_query.seq, service="psi")
    output = open(adress_result, "w")
    psiblast_record = NCBIXML.read(psi_blast)
    for alignment in psiblast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < evalue:
                output.write(">"+alignment.title.split(">")[0]+"\n")
                output.write(hsp.sbjct+"\n")
    psi_blast.close()
    output.close()
    
FromQueryToPSIBLAST("Agglutinin.fasta", "result_PSIBLAST.fasta", 0.04)