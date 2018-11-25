class Result:
    
    def __init__(self, query, name, score, upgap, pvalueq, pscore, pqtscore, pvaluet, qbegin, qend, tbegin, tend, qseq, tseq, tal, qal, norm_score, qcov, identity, gaps, ssscore, allength, corr):
        
        self.query = query
        self.name = name
        self.score = score
        self.upgap = ungap
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
        score = ' '*(9-len(str(round(self.score, 3))))+ str(round(self.score, 3))
        ungapped = ' '*(8 - len(str(round(self.ungap, 3)))) + str(round(self.ungap, 3))
        pvalueq = ' '*(9 - len("{:.2E}".format(self.pvalueq))) + "{:.2E}".format(self.pscore)
        pscore = ' '*(7 - len(str(round(self.pscore, 3)))) + str(round(self.pscore, 3))
        pqtscore = ' ' * (7- len(str(round(self.pqtscore, 3)))) + str(round(self.pqtscore, 3))
        pvaluet = ' '*(9 - len("{:.2E}".format(self.pvalueq))) + "{:.2E}".format(self.pscore)
    
        qlength = ' '*(6-len(str(self.qend-self.qbegin))) + str(self.qend-self.qbegin)
        tlength = ' '*(6-len(str(self.tend-self.tbegin))) + str(self.tend-self.tbegin)
        qbeginend = ' '*(9 - len(str(self.qbegin)+'-'+str(self.qend))) + str(self.qbegin+1)+'-'+str(self.qend+1)
        tbeginend = ' '*(9 - len(str(self.tbegin)+'-'+str(self.tend))) + str(self.tbegin+1)+'-'+str(self.tend+1)

        hitname = '   ' + self.name 
        return(score + ungapped + pvalueq + pscore + pqtscore + pvaluet + qlength + tlength + qbeginend + tbeginend + hitname)
        
    def get_alseq(self, al1, al2, seq1, seq2):
    
        L = np.max([len(seq1), len(seq2)])
    
        seqal1 = ""
        seqal2 = ""
    
        for i in range(L):
        
            # si match
            if al1[i] == al2[i]:
                seqal1 += seq1[i]
                seqal2 += seq2[i]
            
            # si gap dans seq1
            elif al1[i] > al2[i]:
                seqal1 += seq1[i]
                seqal2 += "-"
            
            # si gap dans seq2
            elif al1[i] < al2[i]:
                seqal1 += "-"
                seqal2 += seq2[i]
    
        return(seqal1, seqal2)   
        
    
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