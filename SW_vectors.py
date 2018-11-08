import numpy as np

def normalise_SW(matrix, mu, sigma):
    matrix -= mu
    matrix /= sigma
    
    return(matrix)

def forward_SW(matrix, gp, mm):
    x,y = 1,1
    start, end = matrix.shape[0], matrix.shape[1]
    scores = np.zeros((start, end))
    
    # initialisation 1ère ligne et 1ère colonne
    for i in range(start):
        scores[i][0] = matrix[i][0] - gp
    for i in range(end):
        scores[0][i] = matrix[0][i] - gp
    
    
    for i in range(1, start):
        for j in range(1, end):
            possible = [matrix[i-1][j-1] + mm, matrix[i][j-1] + gp, matrix[i-1][j] + gp]
            scores[i][j] += np.max(possible)
    
    return(scores)
        

def traceback_SW(matrix, give_alignment = True):
    
    # trouver la valeur maximale
    x, y = np.unravel_index(np.argmax(matrix, axis=None), matrix.shape)[0], np.unravel_index(np.argmax(matrix, axis=None), matrix.shape)[1]
    traceback = [(x,y)]
    align_seq1 = [x]
    align_seq2 = [y]
    
    # tant que l'on est pas arrivé au bout de la matrice
    while x > 0 and y > 0:
        
        possible = [matrix[x-1][y-1], matrix[x-1][y], matrix[y-1][x]]
        max_possible = np.argmax(possible)
         
        if max_possible == 0:
            x -= 1
            y -= 1
        elif max_possible == 1:
            x -= 1
        elif max_possible == 2:
            y -= 1
        
        traceback.append((x,y))
        align_seq1.append(x)
        align_seq2.append(y)
        
    if give_alignment:
        return(align_seq1, align_seq2)
    else:
        return(traceback)