class Result:
    
    def __init__(self, query = None, name = None, score = None, ungap = None, pvalueq = None, pscore = None, pqtscore = None, pvaluet = None, qbegin = None, qend = None, tbegin = None, tend = None, qseq = None, tseq = None, tal = None, qal = None, norm_score = None, qcov = None, identity = None, gaps = None, ssscore = None, allength = None, corr = None):
        
        self.query = query
        self.name = name
        self.score = score
        self.ungap = ungap
        self.pvalueq = pvalueq
        self.pscore = pscore
        self.pqtscore = pqtscore
        self.pvaluet = pvaluet
        self.qbegin = qbegin
        self.qend = qend
        self.tbegin = tbegin
        self.tend = tend
        self.qseq = qseq
        self.tseq = tseq
        self.tal = tal
        self.qal = qal
        self.norm_score = norm_score
        self.qcov = qcov
        self.identity = identity
        self.gaps = gaps
        self.ssscore = ssscore
        self.allength = allength
        self.corr = corr
        
    def print_result_1(self):
		if self.score != None:
			score = ' '*(9-len(str(round(self.score, 3))))+ str(round(self.score, 3))
		else:
			score = ' ' * 12
		if self.ungap != None:
			ungapped = ' '*(8 - len(str(round(self.ungap, 3)))) + str(round(self.ungap, 3))
		else:
			ungapped = ' ' * 12
		if self.pvalueq != None:
			pvalueq = ' '*(9 - len("{:.2E}".format(self.pvalueq))) + "{:.2E}".format(self.pscore)
		else:
			pvalueq = ' '*12
		if self.pscore != None:
			pscore = ' '*(7 - len(str(round(self.pscore, 3)))) + str(round(self.pscore, 3))
		else:
			pscore = ' '*12
		if self.pqtscore:
			pqtscore = ' ' * (7- len(str(round(self.pqtscore, 3)))) + str(round(self.pqtscore, 3))
		else:
			pqtscore = ' '*12
		if self.pvaluet != None:
			pvaluet = ' '*(9 - len("{:.2E}".format(self.pvalueq))) + "{:.2E}".format(self.pscore)
		else:
			pvaluet = ' '*12
		if self.qlen != None:
			qlength = ' '*(6-len(str(self.qend-self.qbegin))) + str(self.qend-self.qbegin)
		else:
			qlength = ' '*12
		if self.tlen != None:
			tlength = ' '*(6-len(str(self.tend-self.tbegin))) + str(self.tend-self.tbegin)
		else:
			tlength = ' '*12
		if self.qbegin != None and self.qend != None:
			qbeginend = ' '*(9 - len(str(self.qbegin)+'-'+str(self.qend))) + str(self.qbegin+1)+'-'+str(self.qend+1)
		else: 
			qbeginend = ' '*12
		if self.tbegin != None and self.tend != None:
			tbeginend = ' '*(9 - len(str(self.tbegin)+'-'+str(self.tend))) + str(self.tbegin+1)+'-'+str(self.tend+1)
		else:
			tbeginend = ' '*12
        hitname = '   ' + self.name 
        return(score + ungapped + pvalueq + pscore + pqtscore + pvaluet + qlength + tlength + qbeginend + tbeginend + hitname)
        
	def print_alignement(al1, al2, seq1, seq2):
    
        L = len(al1)
        print(L)
    
        seqal1 = ""
        seqal2 = ""
        if len(seq2)>len(seq1):
            Seq1=seq1
            Seq2=seq2
        else:
            Seq1=seq2
            Seq2=seq1
        for i in range(L):
            j=al1[i]
            k=al2[i]
       
           if Seq1[j]==Seq2[k]:
               seqal1 += Seq1[j]
               seqal2 += Seq2[k]
        
           else:
               if Seq1[j]!=Seq2[k]:
                   if j>k:
                       seqal1 += Seq1[j]
                       seqal2 += "-"
                   else:
                       seqal1 += "-" 
                       seqal2 += Seq2[k]
        return(seqal1, seqal2)

TestAlign1, TestAlign2=print_alignement(align1, align2, seq1, seq2)
    
    def print_result_2(self):
        names = 'Alignment : ' + self.query + ', ' + str(self.qseq)+' aa. vs  ' + self.name + '\n'
        info = "Score :  " + str(round(self.score, 3))+ " | Normalized score :    " + str(round(self.norm_score, 3))+ " | Query coverage : " + str(round(self.norm_score, 2)) 
        info += "% | Identity :   " + str(round(self.identity, 2)) + "% | Gaps :     "
        info += str(round(self.gaps, 2)) + "% | SS Score :    " + str(round(self.ssscore, 2))
        info += " | Alignment length :    " + str(allength) + " | Corr score :   " + str(round(self.corr, 2)) + "\n\n"
        seqs = self.get_alseq(self, self.qal, self.tal, self.qseq, self.tseq)
        align_query = "Query       " + str(self.qbegin) + seqs[0] + "    " + str(self.qend) + '\n'
        template_query = "Template    "+ str(self.tbegin) + seqs[1] + "    " + str(self.tend) + '\n'
        
        return(names+info+align_query+template_query)



def print_all(name_query, len_query, list_results, output_file):
    
    to_write = "*** HITS RANKED *** \n"
    to_write += "\n"
    to_write += "SEQUENCE QUERY FILE : " + name_query + ", " + len_query + " aa.\n"
    to_write += "\n"
    to_write += "#| Score | Ungaped_score | Pvalue_Q | Pscore | PQTscore | P-Value_T | Q. Length | T. Length | Q. begin-end |  T. begin-end | HITS"
    to_write += '\n---------------------------------------------------------'
    
    for i in range(len(list_results)):
        
        first_part = ' '*(5-len(str(i+1))) + str(i+1) + list_results[i].print_result_1()
        to_write += first_part
        
    to_write += '\n\n\n'
    to_write += '*** ALIGNMENTS DETAILS ***\n'
    
    for i in range(len(list_results)):
        to_write += '\n'
        to_write += 'No '+str(i+1)+'\n'
        to_write += list_results[i].print_result_2()
        
    with open(output_file, 'w') as f:
        for line in to_write:
            f.write(line)
