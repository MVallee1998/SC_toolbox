def low(mat, j):
    """calculates low(j), low(j) is the highest index i with a non-zero element"""
    if mat[j] == []:
        return 0
    else:
        return mat[j][-1]


def create_D(Ksimplicies):
    '''creates the matrices of the boundary operator'''
    D = list(Ksimplicies)
    return (D)


def fusion(T1, T2):
    '''this is method to sum 2 Z2 vectors represented with their indexes of non-0'''
    if T1 == []: return T2
    if T2 == []: return T1
    if T1[0] < T2[0]:
        return [T1[0]] + fusion(T1[1:], T2)
    elif T1[0] == T2[0]:
        return fusion(T1[1:], T2[1:])
    else:
        return [T2[0]] + fusion(T1, T2[1:])


def reduction_block(block):
    '''This function reduces a given matrix (finding linear relations between its columns'''
    listLow = [low(block, j) for j in range(len(block))]
    for j in range(len(block)):
        exists = True
        while listLow[j] != 0 and exists:
            exists = False
            j0 = 0
            while j0 < j and not exists:
                if listLow[j0] == listLow[j]:
                    exists = True
                    block[j] = fusion(block[j], block[j0])
                    listLow[j] = low(block, j)
                else:
                    j0 = j0 + 1


def reduction(mat):
    """reduces the mat matrix, which should be  a list of matrices to be reduced"""
    dimension = len(mat)
    for d in range(1, dimension):
        reduction_block(mat[d])


def computeBettiNumbers(D):
    '''Computes the Z2 betti numbers of a simplicial complex represented by its boundary operator matrix'''
    dimension = len(D)
    if dimension == 0:
        return []
    reduction(D)
    kernelDim = []
    imageDim = []
    for d in range(dimension + 1):
        if d == 0:
            kernelDim.append(len(D[d]))
        elif d == dimension:
            imageDim.append(0)
        else:
            kernelDim.append(0)
            imageDim.append(0)
            for i in range(len(D[d])):
                if low(D[d], i) == 0:
                    kernelDim[-1] += 1
                else:
                    imageDim[-1] += 1
    bettiNbrs = []
    for k in range(len(kernelDim)):
        bettiNbrs.append(kernelDim[k] - imageDim[k])
    return bettiNbrs
