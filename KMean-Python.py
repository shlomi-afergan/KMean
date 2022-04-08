import sys

args = sys.argv
k = int(args[1])
if len(args) == 5:
    iterations = int(args[2])
    inputFileName = args[3]
    outputFileName = args[4]
else:
    iterations = 200
    inputFileName = args[2]
    outputFileName = args[3]


# Calculate the euclidean distance between two vectors.
def euclideanDistance(vec1: list[float], vec2: list[float]) -> float:
    return (sum([(vec1[i] - vec2[i]) ** 2 for i in range(len(vec1))])) ** 0.5


# Calculate cluster centroid.
def clusterCentroid(vectors: list[list[float]]) -> list[float]:
    vec_len = len(vectors[0])
    centroid = []
    for i in range(vec_len):
        summ = 0
        for lst in vectors:
            summ += lst[i]
        center = summ / len(vectors)
        center = round(center, 4)
        centroid.append(center)
    return centroid


def clustersFar(old: list[list[float]], new: list[list[float]]) -> bool:
    for i in range(len(old)):
        if euclideanDistance(old[i], new[i]) > 0.001:
            return True
    return False


def toStrLst(centroids: list[int]) -> list[list[str]]:
    result = [[addZero(str(x)) for x in lst] for lst in centroids]
    result = [",".join(lst) for lst in result]
    return result


def addZero(num: str):
    numLst = num.split(".")
    if len(numLst[1]) == 3:
        numLst[1] += '0'
    return ".".join(numLst)


def kMeans(k: int, max_iter: int, myInput: str, myOutput: str):
    """
    :param k: The number of desired clusters.
    :param myInput: A text file contains all the points in R**d.
    :param max_iter: max iterations allowed.
    :param myOutput: txt file to output.
    :return: K centers of points.
    """
    # process input
    inputLst = []
    with open(myInput) as f:
        for line in f:
            lst = [float(x) for x in line[:-1].split(",")]
            inputLst.append(lst)
    f.close()

    # Initialize clusters centers.
    centroids = [inputLst[i] for i in range(k)]

    counter = 0

    # Repeat until convergence.
    new_centroids = [[float("inf") for x in lst] for lst in centroids]

    while counter < max_iter and clustersFar(centroids, new_centroids):
        if counter != 0:
            centroids = new_centroids.copy()
        clusters = [[] for i in range(len(centroids))]

        # Assign points to centroid by Euclidean distance.
        for vec in inputLst:
            min_euc = float("inf")
            cluster_idx = len(clusters) + 1
            for i in range(len(centroids)):
                tmp = euclideanDistance(vec, centroids[i])
                if tmp < min_euc:
                    cluster_idx = i
                    min_euc = tmp
            clusters[cluster_idx].append(vec)

        # Update centroids.
        for i in range(len(clusters)):
            new_centroids[i] = clusterCentroid(clusters[i])

        counter += 1  # update counter
    # output result
    strLst = toStrLst(new_centroids)
    with open(outputFileName, 'w') as f:
        for line in strLst:
            f.write(line)
            f.write('\n')
    return new_centroids


kMeans(k, iterations, inputFileName, outputFileName)





