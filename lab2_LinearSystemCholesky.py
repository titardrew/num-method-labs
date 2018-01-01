import numpy as np
import cmath as cm


# SQRT method for both sym/not sym matrix (Cholesky's method)
# The program has an file-type input
# Data is structured like 'A b1 b2 ... bn', where A - is system's matrix, b1...bn - set of right vectors
# Works for hermit matrix + sym

PATH = "lab2_input"


def read_matrix(path):
    try:
        with open(path) as textFile:
            lines = [line.split() for line in textFile]
        return lines
    except FileNotFoundError:
        print("Error 0. File is not found.")


def get_system_data(adv_matrix):
    b = []
    matrix = []

    for j in range(0, len(adv_matrix)):
        matrix.append([])
        for i in range(0, len(adv_matrix)):
            matrix[j].append(float(adv_matrix[i][j]))

    for j in range(len(adv_matrix), len(adv_matrix[0])):
        temp = []
        for i in range(0, len(adv_matrix)):
            temp.append(float(adv_matrix[i][j]))
        b.append(temp)

    return {
            'matrix': matrix,
            'vectors': b
            }


def get_left_matrix(data):

    l = np.zeros((len(data['matrix']), len(data['matrix'])), dtype=complex)
    for i in range(0, len(data['matrix'])):
        for j in range(0, i + 1):
            if i == j:
                l[i][i] = cm.sqrt(data['matrix'][i][i] - sum(l[i][x]**2 for x in range(i)))
            elif i > j:
                if l[j][j] == 0:
                    raise ZeroDivisionError
                l[i][j] = (data['matrix'][i][j] - sum(l[i][x] * l[j][x] for x in range(j))) / l[j][j]
    return l


def get_first_vectors(vectors, left_matrix):
    first_vectors = []
    for vector in vectors:
        temp_vector = np.zeros(len(left_matrix), dtype=complex)
        for i in range(len(left_matrix)):
            if left_matrix[i][i] == 0:
                raise ZeroDivisionError
            temp_vector[i] = ((1 / left_matrix[i][i]) *
                                (vector[i] - sum(left_matrix[i][k] * temp_vector[k] for k in range(i))))
        first_vectors.append(temp_vector)
    return first_vectors


def get_second_vectors(first_vectors, left_matrix):
    second_vectors = []
    for vector in first_vectors:
        temp_vector = np.zeros(len(left_matrix), dtype=complex)
        for i in reversed(range(0, len(left_matrix))):
            if left_matrix[i][i] == 0:
                raise ZeroDivisionError
            temp_vector[i] = ((1 / left_matrix[i][i]) *
                                (vector[i] - (sum(np.conj(left_matrix[k][i]) * temp_vector[k] for k in range(i + 1, len(vector))))))
        second_vectors.append(temp_vector)
    return second_vectors


def det_by_l_matrix(l):
    print("-----------------")
    print("Det(A) = %f" % np.prod(list(l[x][x] for x in range(len(l)))).real ** 2)
    print("-----------------")


def cholesky_(data):
    try:
        l = get_left_matrix(data)
        fu = get_second_vectors(get_first_vectors(data['vectors'], l), l)
        return fu
    except ZeroDivisionError:
        print("Wrong matrix..")
        return None


def print_inverse_matrix(data):
    print("-----------------")
    print("\t\t\t\t\t\t\tA^-1: \n")
    ident = []
    for i in range(len(data['matrix'])):
        ident.append([])
        for j in range(len(data['matrix'])):
            ident[i].append(1 if i==j else 0)
    data['vectors'] = ident
    m = cholesky_(data)
    for i in range(len(m)):
            print(np.transpose(m)[i].real)
    print("-----------------")
    print("-----------------")
    print("\t\t\t\t\t\t\tA•A^-1: \n")
    m = np.array(data['matrix']).dot(np.array(np.transpose(m)))
    for i in range(len(m)):
        print(m[i].real)
    print("-----------------")


def print_roots(m, data):
    print("-----------------")
    print("Roots are: ")
    j = 0
    for str in m:
        i = 0
        print()
        for element in str:
            i += 1
            print("x[%d] = %f" % (i, element.real))
        print("Residual vector :")
        print((data['vectors'][j-1] - np.array(data['matrix']).dot(str)).real)
        print()
        j += 1
    print("-----------------")


def print_left_matrix(data):
    print("-----------------")
    print("\t\t\t\t\t\t\tL:\n")
    for l in get_left_matrix(data).real:
        print(l)
    print("-----------------")
    print("-----------------")
    print("\t\t\t\t\t\t\tL(T*):\n")
    for l in np.conj(np.transpose(get_left_matrix(data).real)):
        print(l)
    print("-----------------")
    print("-----------------")
    print("\t\t\t\t\t\t\tL•L(T*):\n")
    for l in get_left_matrix(data).dot(np.conj(np.transpose(get_left_matrix(data)))).real:
        print(l)
    print("-----------------")

# Read and build data dict
data = get_system_data(read_matrix(PATH))

# Print det(L)
det_by_l_matrix(get_left_matrix(data))

# Print L, L(T*), L•L(T*)
print_left_matrix(data)

# Print roots with residual vectors
print_roots(cholesky_(data), data)

# Print A^-1, A•A^-1
print_inverse_matrix(data)

'''[[ 0.1556313  -0.04963902 -0.00347847 -0.02841869  0.00248072]
 [-0.04963902  0.35024567 -0.11289702  0.11865354  0.09875967]
 [-0.00347847 -0.11289702  0.33632411 -0.37827949 -0.33299015]
 [-0.02841869  0.11865354 -0.37827949  0.90287604  0.75116793]
 [ 0.00248072  0.09875967 -0.33299015  0.75116793  0.83076572]]'''