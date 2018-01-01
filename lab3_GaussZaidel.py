import numpy as np
PATH = "lab3_input"


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


def solve_zaidel(matrix, right_vector, eps):
    value = np.zeros(len(matrix))
    iter = 0
    finished = False
    while not finished:

        next_value = value.copy()
        for i in range(len(matrix)):
            sum1 = sum(matrix[i][j] * next_value[j] for j in range(i))
            sum2 = sum(matrix[i][j] * value[j] for j in range(i + 1, len(matrix)))
            next_value[i] = (right_vector[i] - sum1 - sum2) / matrix[i][i]
        iter += 1

        print("iter = ", iter, "\t\tValue = ", next_value, "")
        print('\t\t\t\tResidual vector = [',end='')
        for i in range(len(matrix)):
            print(next_value[i] - value[i], end=' ')
        print(']', end='\n')
        error = np.sqrt(sum((next_value[i] - value[i]) ** 2 for i in range(len(matrix))))
        print("\t\t\t\tError = ", error, end='\n\n')
        finished = error <= eps
        value = next_value
    return value


data = get_system_data(read_matrix(PATH))
print()
try:
    print("Value = ", solve_zaidel(np.transpose(data['matrix']), data['vectors'][0], 0.000001))
except Exception:
    print("Exception matrix...")
