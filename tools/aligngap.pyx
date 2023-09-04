# distutils: language=c++
#cython: language_level=3
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef int[2][2000][2000] dp

def align(str a, str b):
    cdef int na, nb
    cdef int gap_open=-1, gap_cont=0, gap_malign=-20
    na, nb = len(a), len(b)
    cdef int i, j
    for i in range(2):
        for j in range(0, na+1):
            dp[i][j][0] = i * gap_open
        for j in range(0, nb+1):
            dp[i][0][j] = i * gap_open
    dp[1][0][0] = 0
    
    for j in range(1, na+1):
        for k in range(1, nb+1):
            for i in range(1, -1, -1):
                dp[i][j][k] = -1000
                if a[j-1] == b[k-1]:
                    dp[i][j][k] = max(dp[1][j-1][k-1]+1, dp[1][j-1][k-1]+gap_malign)
                    dp[i][j][k] = max(dp[i][j][k], dp[0][j-1][k]+(i*gap_open))
                    dp[i][j][k] = max(dp[i][j][k], dp[0][j][k-1]+(i*gap_open))
                else:
                    dp[i][j][k] = max(dp[0][j-1][k]+(i*gap_open), dp[0][j][k-1]+(i*gap_open))
                    dp[i][j][k] = max(dp[i][j][k], dp[1][j-1][k-1]+gap_malign)

    cdef int c, d
    c, d = len(a), len(b)
    cdef vector[int] resulta, resultb
    i = 1
    while c > 0 or d > 0:
        if c > 0 and dp[i][c][d] == dp[0][c-1][d]+i*gap_open:
            c -= 1
            i = 0
        elif d > 0 and dp[i][c][d] == dp[0][c][d-1]+i*gap_open:
            d -= 1
            i = 0
        elif c > 0 and d > 0 and (dp[i][c][d] == dp[1][c-1][d-1]+gap_malign or dp[i][c][d] == dp[1][c-1][d-1]+1):
            resulta.push_back(c-1)
            resultb.push_back(d-1)
            c -= 1
            d -= 1
            i = 1
        else:
            c -= 1
            d -= 1
    return resulta[::-1], resultb[::-1]