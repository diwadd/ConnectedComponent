import math

class ConnectedComponent:
    def permute(self, matrix):
        ret = []
        S = int(math.sqrt(len(matrix)))
        for i in range(S):
            ret.append(S - 1 - i)
        return ret

# -------8<------- end of solution submitted to the website -------8<-------

import sys
M = int(raw_input())
matrix = []
for i in range(M):
    matrix.append(int(raw_input()))

cc = ConnectedComponent()
ret = cc.permute(matrix)
print len(ret)
for num in ret:
    print num
sys.stdout.flush()
