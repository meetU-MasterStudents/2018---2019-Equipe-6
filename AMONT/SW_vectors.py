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



def local_alignment_affine_gap_penalty(scoring_matrix, gap_opening_penalty, gap_extension_penalty):
    '''
    Returns the score and local alignment substrings for strings v, w with the
    given scoring matrix, gap opening penalty, and gap extension penalty.
    '''
    v,w=scoring_matrix.shape
    # Initialize the matrices.
    S_lower = [[0 for j in range(w+1)] for i in range(v+1)]
    scoresSW = [[0 for j in range(w+1)] for i in range(v+1)]
    S_upper = [[0 for j in range(w+1)] for i in range(v+1)]

    # Initialize the maximum score below the lowest possible score.
    max_score = -1
    max_i, max_j = 0, 0

    # Fill in the Score and Backtrack matrices.
    for i in range(1, v+1):
        for j in range(1, w+1):
            S_lower[i][j] = max([S_lower[i-1][j] - gap_extension_penalty, scoresSW[i-1][j] - gap_opening_penalty])
            S_upper[i][j] = max([S_upper[i][j-1] - gap_extension_penalty, scoresSW[i][j-1] - gap_opening_penalty])
            middle_scores = [S_lower[i][j], scoresSW[i-1][j-1] + scoring_matrix[i-1,j-1], S_upper[i][j], 0]
            scoresSW[i][j] = max(middle_scores)

            if scoresSW[i][j] > max_score:
                max_score = scoresSW[i][j]
                max_i, max_j = i, j

#    # Initialize the indices to start at the position of the high score.
#    i, j = max_i, max_j
#
#    # Initialize the aligned strings as the input strings up to the position of the high score.
#    v_aligned, w_aligned = v[:i], w[:j]
#
#    # Backtrack to start of the local alignment starting at the highest scoring cell.
#    # Note: the solution format specifically asks for substrings, so no indel insertion necessary.
#    while backtrack[i][j] != 3 and i*j != 0:
#        if backtrack[i][j] == 0:
#            i -= 1
#        elif backtrack[i][j] == 1:
#            i -= 1
#            j -= 1
#        elif backtrack[i][j] == 2:
#            j -= 1
#
#    # Cut the strings at the ending point of the backtrack.
#    v_aligned = v_aligned[i:]
#    w_aligned = w_aligned[j:]

    return(np.array(scoresSW))  



def traceback_SW(ScoresSW, give_alignment = True):
    
    # trouver la valeur maximale
    x, y = np.unravel_index(np.argmax(ScoresSW, axis=None), ScoresSW.shape)[0], np.unravel_index(np.argmax(ScoresSW, axis=None), ScoresSW.shape)[1]
    traceback = [(x,y)]
    #print(traceback)
    align_seq1 = [x]
    align_seq2 = [y]
    ngap = 0
       
    
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





