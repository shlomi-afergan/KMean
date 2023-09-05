import numpy as np
import os
import sys
import enum
import mykmeansspf


args = sys.argv
MAX_ITER = 300
EPS = 0


class GOALS(enum.Enum):
    spk = 'spk'
    wam = 'wam'
    ddg = 'ddg'
    lnorm = 'lnorm'
    jacobi = 'jacobi'


def invalid_input():
    print("Invalid Input!")
    exit(1)


def validateArg(arg_str, mode):
    if mode == "k":
        try:
            arg = int(arg_str)
        except:
            invalid_input()
        arg = int(arg_str)
        if arg < 0:
            invalid_input()

    elif mode == "file":
        length = len(arg_str)
        if length <= 4 or (arg_str[length-4:length] != ".txt" and arg_str[length-4:length] != ".csv"):
            invalid_input()
    elif mode == "goal":
        try:
            if GOALS(arg_str) not in GOALS: invalid_input()
        except:
            invalid_input()


def print_mat(matrix):
    for row in matrix:
        print(",".join(addZero(str(d)) for d in row))


# Calculate the euclidean distance between two vectors.
def euclideanDistance(vec1, vec2):
    return sum([(vec1[i] - vec2[i]) ** 2 for i in range(len(vec1))])


def addZero(num):
    numLst = num.split(".")
    if len(numLst[1]) == 3:
        numLst[1] += '0'
    elif len(numLst[1]) == 2:
        numLst[1] += '00'
    elif len(numLst[1]) == 1:
        numLst[1] += '000'
    return ".".join(numLst)


def writeFile(vec_array):
    str_lst = [[addZero(str(x)) for x in lst] for lst in vec_array]
    str_lst = [",".join(lst) for lst in str_lst]
    with open("file_in.txt",'w') as f:
        for line in str_lst:
            f.write(line)
            f.write('\n')
    f.close()


def mat_to_lst(centroids):
    res = []
    for lst in centroids:
        for num in lst:
            res.append(num)
    return res


def kmeans_pp():
    np.random.seed(0)
    idxs = [i for i in range(len(matrix))]
    vec_array = np.asarray(matrix)
    dimension = vec_array.shape[1]
    number_of_rows = vec_array.shape[0]
    if k > number_of_rows:
        invalid_input()
    first = np.random.choice(number_of_rows)
    dl_dict = dict()
    centroids = [vec_array[first]]
    indices = [first]
    while len(centroids) < k:
        for i in range(number_of_rows):
            tmp = float('inf')
            for center in centroids:
                dl = euclideanDistance(vec_array[i], center)
                if dl < tmp:
                    tmp = dl
            dl_dict[i] = tmp
        denominator = sum(dl_dict.values())
        for key in dl_dict:
            if denominator != 0:
                dl_dict[key] = dl_dict[key] / denominator
        sel = np.random.choice(number_of_rows, p=list(dl_dict.values()))
        indices.append(sel)
        centroids.append(vec_array[sel])
    writeFile(vec_array)
    result = []
    for i in range(k):
        result.append(int(idxs[indices[i]]))
    last_centroids = mykmeansspf.fit_kmeanspp(k, MAX_ITER, EPS, "file_in.txt", indices)[:k*dimension]
    last_centroids = np.asarray(last_centroids).reshape(k, dimension).round(4)
    indices_str = ",".join(str(j) for j in result)
    print(indices_str)
    for cent in last_centroids:
        print(",".join(addZero(str(j)) for j in cent))
    os.remove("file_in.txt")


if len(args) != 4:
    invalid_input()
validateArg(args[1], "k")
validateArg(args[2], "goal")
k = int(args[1])
goal = GOALS(args[2])
file_name = args[3]
validateArg(file_name, "file")


try:
    matrix = mykmeansspf.fit_by_goal(k, goal.name, file_name)
    if goal == GOALS.spk:
        matrix = np.asarray(matrix)
        k = matrix.shape[1]
        kmeans_pp()
    else:
        matrix = np.asarray(matrix).round(4)
        print_mat(matrix)

except:
    print("An Error Has Occurred")
    exit(1)
