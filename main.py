# Calculate the euclidean distance between two vectors.
def euclideanDistance(vec1: list, vec2: list) -> int:
    return (sum([(vec1[i] - vec2[i]) ** 2 for i in range(len(vec1))])) ** 0.5


# Calculate cluster centroid.
def clusterCentroid(vectors: list[list]) -> list[float]:
    vec_len = len(vectors[0])
    centroid = []
    for i in range(vec_len):
        summ = 0
        for lst in vectors:
            summ += lst[i]
        center = summ / len(vectors)
        center = "%.4f" % center
        center = float(center)
        centroid.append(center)
        # append(summ / len(vectors))
    return centroid


def clustersFar(old: list[list], new: list[list]) -> bool:
    for i in range(len(old)):
        if euclideanDistance(old[i], new[i]) > 0.01:
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


print(clusterCentroid([[0, 0], [4, 4], [0, 4], [4, 0]]))  # -----------------Code Review---------------------------
print(clustersFar([[0, 0], [4, 4], [0, 4], [4, 0]], [[0,3], [3,4], [2,2], [1,4]]))

print(euclideanDistance([7, 11], [11, 8]))  # -----------------------------Code Review-----------------------------


def kMeans(k: int, myInput: str, max_iter=200):
    """
    :param k: The number of desired clusters.
    :param myInput: A text file contains all the points in R**d.
    :param max_iter: max iterations allowed.
    :return: K centers of points.
    """
    # process input
    inputLst = []
    with open(myInput) as f:
        for line in f:
            lst = [float(x) for x in line[:-1].split(",")]
            inputLst.append(lst)
    f.close()

    #  print(inputLst)  # -----------------------------Code Review---------------------------------------

    # Initialize cluster centers.
    centroids = [inputLst[i] for i in range(k)]
    # print("centroids are: ", centroids)  # -----------------------------Code Review---------------------------------

    counter = 0

    # Repeat until convergence.
    # new_centroids = [[float("inf") for x in vec] for vec in centroids]
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
        # print("cent; ", centroids)
        # print("new cent: ", new_centroids)

        # update counter
        counter += 1
    # output result
    strLst = toStrLst(new_centroids)
    with open('myFile.txt', 'w') as f:
        for line in strLst:
            f.write(line)
            f.write('\n')
    print("result is: ", new_centroids)
    return new_centroids


m = kMeans(3, "input_1.txt")  # -----------------------------Code Review---------------------------------------


print(toStrLst(m))