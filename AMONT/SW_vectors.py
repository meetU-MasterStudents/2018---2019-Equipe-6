import numpy as np

def normalise_SW(matrix, mu, sigma):
    matrix -= mu
    matrix /= sigma
    
    return(matrix)


def forward_SW(MatrixN, gp, match):
    x,y = 1,1
    start, end = MatrixN.shape[0], MatrixN.shape[1]
    scoresSW = np.zeros((start, end))
    
    # initialisation 1ère ligne et 1ère colonne
    for i in range(start):
        scoresSW[i][0] = 0 
    for i in range(end):
        scoresSW[0][i] = 0 
    
    
    for i in range(1, start):
        for j in range(1, end):
            possible = [0,MatrixN[i-1][j-1] + match, MatrixN[i][j-1] + gp, MatrixN[i-1][j] + gp]
            scoresSW[i][j] += np.max(possible)
    
    return(scoresSW)
        



def traceback_SW(ScoresSW, give_alignment = True):
    
    # trouver la valeur maximale
    x, y = np.unravel_index(np.argmax(ScoresSW, axis=None), ScoresSW.shape)[0], np.unravel_index(np.argmax(ScoresSW, axis=None), ScoresSW.shape)[1]
    traceback = [(x,y)]
    #print(traceback)
    align_seq1 = [x]
    align_seq2 = [y]
    ngap = 0
    
    #print(matrix[x][y])
    #print(x,y)    
    
    # tant que l'on est pas arrivé au bout de la matrice
    while x > 0 and y > 0:
        possible = [ScoresSW[x-1][y-1], ScoresSW[x-1][y], ScoresSW[x][y-1]]  #Gabriela
        max_possible = np.argmax(possible)
        #print(max_possible)
        if max_possible == 0:
            x -= 1
            y -= 1
        elif max_possible == 1:
            x -= 1
            ngap += 1
        elif max_possible == 2:
            y -= 1
            ngap += 1
        
        traceback.append((x,y))
        align_seq1.append(x)
        align_seq2.append(y)
        
    # 1ere et derniere pos alignée
    start_query = np.min(align_seq1)
    start_template = np.min(align_seq2)
    
    end_query = np.max(align_seq1)
    end_template = np.max(align_seq2)
    
        
    if give_alignment:
        return(np.max(ScoresSW),ngap, start_query, start_template, end_query, end_template, sorted(align_seq1), sorted(align_seq2))
    else:
        return(traceback)