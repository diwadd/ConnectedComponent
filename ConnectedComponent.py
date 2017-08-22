import math


def eprint(s="", end="\n", verbose=True):
    """
    This is a custom print function since
    only logging into the standard error is allowed.
    """
    if verbose == True:
        sys.stderr.write(s + end)
    else:
        pass


def resize_string(s, l=2):

    while (len(s) < l):
        s = " " + s
    return s


def print_matrix(mat):

    s = len(mat)
    for i in range(s):
        for j in range(s):
            if type(mat[i][j]) == str:
                number = mat[i][j]
            else:
                number = resize_string(str(mat[i][j]))
            eprint(number, end=" ")
        eprint()


def convert_row_matrix_to_2d_array(matrix):

    s = int(math.sqrt(len(matrix)))
    mat = [[matrix[j*s + i] for i in range(s)] for j in range(s)]

    return mat


def permute_matrix(mat, per):

    s = len(mat)
    mat_p = [[mat[per[i]][per[j]] for j in range(s)] for i in range(s)]

    return mat_p    



class ConnectedComponent:
    def permute(self, matrix):

        mat = convert_row_matrix_to_2d_array(matrix)

        """
        mat2 = ["00", "01", "02", "03", 
                "10", "11", "12", "13",
                "20", "21", "22", "23",
                "30", "31", "32", "33"]
        """

        mat2 = [" 1", "-1", " 0", 
                "-1", " 1", " 0",
                " 0", " 0", " 0"]
        mat2 = convert_row_matrix_to_2d_array(mat2)


        eprint("Original matrix: ")
        print_matrix(mat2)

        eprint("Permuted matrix: ")
        #mat_p = permute_matrix(mat2, [1, 0, 2, 3])
        mat_p = permute_matrix(mat2, [0, 1, 2])
        print_matrix(mat_p)

        eprint("Permuted matrix: ")
        mat_p = permute_matrix(mat_p, [0, 2, 1])
        #mat_p = permute_matrix(mat_p, [0, 1, 3, 2])
        print_matrix(mat_p)

        eprint("Permuted matrix: ")
        mat_p = permute_matrix(mat_p, [1, 0, 2])
        #mat_p = permute_matrix(mat_p, [0, 1, 3, 2])
        print_matrix(mat_p)

        eprint("Permuted matrix: ")
        mat_p = permute_matrix(mat_p, [1, 2, 0])
        #mat_p = permute_matrix(mat_p, [0, 1, 3, 2])
        print_matrix(mat_p)

        eprint("Permuted matrix: ")
        mat_p = permute_matrix(mat_p, [2, 0, 1])
        #mat_p = permute_matrix(mat_p, [0, 1, 3, 2])
        print_matrix(mat_p)

        eprint("Permuted matrix: ")
        mat_p = permute_matrix(mat_p, [2, 1, 0])
        #mat_p = permute_matrix(mat_p, [0, 1, 3, 2])
        print_matrix(mat_p)



        mat3 = ["00", "01", "02", "03", 
                "10", "11", "12", "13",
                "20", "21", "22", "23",
                "30", "31", "32", "33"]
        mat3 = convert_row_matrix_to_2d_array(mat3)

        eprint("Permuted matrix: ")
        mat_p = permute_matrix(mat3, [2, 3, 1, 0])
        print_matrix(mat_p)


        ret = []
        S = int(math.sqrt(len(matrix)))
        for i in range(S):
            ret.append(S - 1 - i)
        #return [4, 3, 8, 7, 2, 0, 1, 9, 6, 5]
        return [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        #return ret

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
