import numpy as np

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
        if self.qbegin != None and self.qend != None:
            qlength = ' '*(6-len(str(self.qend-self.qbegin))) + str(self.qend-self.qbegin)
        else:
            qlength = ' '*12
        if self.tbegin != None and self.tend != None:
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
        
    
    def print_result_2(self):
        if self.query != None and self.qseq != None and self.name != None:
            names = 'Alignment : ' + self.query + ', ' +' aa. vs  ' + self.name + ' :  \n'
        else:
            names = 'Alignment : ' + '          ' + ', ' + '          '+' aa. vs  ' + '          ' + '\n'
        if self.score != None and self.norm_score != None and self.qcov != None:
            #info = "Score :  " + str(round(self.score, 3))+ " | Normalized score :    " + str(round(self.norm_score, 3))+ " | Query coverage : " + str(round(self.qcov, 2)) 
            info = "Score :  " + str(round(self.score, 3))+ " | Normalized score :    " + str(round(self.norm_score, 3))+ " | Query coverage : " + str(round(self.qcov, 2)*100) + "%"
        else:
            info = "Score :  " + '       '+ " | Normalized score :    " + '       ' + " | Query coverage : " + '       '
        if self.identity != None:
            #info += "% | Identity :   " + str(round(self.identity, 2)) + "% | Gaps :     "
            info += "% | Identity :   " + str(round(self.identity, 2)*100) + "% | Gaps :     " 
        else:
            info += "% | Identity :   " + '        ' + "% | Gaps :     "
        if self.gaps != None and self.ssscore != None:
            #info += str(round(self.gaps, 2)) + "% | SS Score :    " + str(round(self.ssscore, 2))
            info += str(round(self.gaps, 2)*100) + "% | SS Score :    " + str(round(self.ssscore, 2))
        else:
            info += '        ' + "% | SS Score :    " + '        '
        if self.allength != None and self.corr != None:
            info += " | Alignment length :    " + str(self.allength) + " | Corr score :   " + str(round(self.corr, 2)) + "\n\n"
        else: 
            info += " | Alignment length :    " + '        ' + " | Corr score :   " + '        ' + "\n\n"
        if self.qseq != None and self.tseq != None:
            seqs = [self.qseq, self.tseq]
        else:
            seqs = 'N/A'
        if self.qbegin != None and self.qend != None: 
            align_query = "Query       " + str(self.qbegin) +' '*(5 - len(str(self.qbegin)))+ seqs[0] + "    " + str(self.qend) + '\n'
        else:
            align_query = "Query       " + '           ' + "    " + '         ' + '\n'
        if self.tbegin != None and self.tend != None:
            template_query = "Template    "+ str(self.tbegin) +' '*(5 - len(str(self.tbegin)))+ seqs[1] + "    " + str(self.tend) + '\n'
        else:
            template_query = "Template    "+ '           ' + "    " + '           ' + '\n'
        
        return(names+info+align_query+template_query)



def print_all(name_query, len_query, list_results, output_file):

    len_query = str(len_query)
    
    to_write = "*** HITS RANKED *** \n"
    to_write += "\n"
    to_write += "SEQUENCE QUERY FILE : " + name_query + ", " + len_query + " aa.\n"
    to_write += "\n"
    to_write += "#| Score | Ungaped_score | Pvalue_Q | Pscore | PQTscore | P-Value_T | Q. Length | T. Length | Q. begin-end |  T. begin-end | HITS"
    to_write += '\n---------------------------------------------------------\n'
    
    for i in range(len(list_results)):
        
        first_part = ' '*(5-len(str(i+1))) + str(i+1) + list_results[i].print_result_1()
        to_write += first_part
        to_write += '\n'
        
    to_write += '\n\n\n'
    to_write += '*** ALIGNMENTS DETAILS ***\n'
    
    for i in range(len(list_results)):
        to_write += '\n'
        to_write += 'No '+str(i+1)+'\n'
        to_write += list_results[i].print_result_2()
        
    with open(output_file, 'w') as f:
        for line in to_write:
            f.write(line)