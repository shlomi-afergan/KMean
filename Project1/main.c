/*Preprocessor statements*/
#include <stdio.h>
#include "math.h"
#include "stdlib.h"

/* functions declaration*/
double ** readFile(int k, char *fileName);
int kMean(int k, int max_iter,char *fileName, char *outputName);
double euclideanDistance(const double* vec1, const double* vec2);
int centroidsFar(double **old, double** new, int k);
void copyCentroids(double **target, double **source, int k);
double** createCentroid(int k);
void updateSumCentroid(double* vec1, const double* vec2);
void resetCentroids(double** centroids, int k);
void updateCentroids(double **centroids, const int arr[], int k);
void resetSizes(int* arr, int k);
void writeFile(double **centroids, int k, char* outputName);

/* Global variables*/
int dimensions = 1;
int lineCounter = 0;


/** The function reads all the vectors from the input file.
 *  Counts the number of lines of text as the number of vectors.
 *  Counts the number of values in each vector as the number of dimension **/
double ** readFile(int k, char *fileName) {
    FILE *fp = NULL;
    char ch;
    double **inputData, *buff;
    int size_of_buff, i=0, buff_idx=0, vec_idx, line_idx;

    /* Get Number of vectors and dimensions in the file.*/
    fp = fopen(fileName, "r");
    if(!fp) {
        printf("Invalid Input!\n");
        exit(1);
    }
    while((ch=fgetc(fp))!= EOF) {
        if (lineCounter == 0 && ch == ',')
            dimensions++;
        if(ch=='\n')
            lineCounter++;
    }
    fclose(fp);

    /* Create list of pointers in the right size.*/
    fp = fopen(fileName, "r");
    inputData = (double **)calloc(lineCounter, sizeof(double *));
    if(!inputData) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    size_of_buff = dimensions*lineCounter;
    buff = calloc(size_of_buff, sizeof(double));
    if(!buff) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    while(fscanf(fp, "%lf,", &buff[i++])!=EOF){
    }
    fclose(fp);

    /* copy buff to my inputData matrix. */
    for(line_idx=0; line_idx<lineCounter; line_idx++) {
        double *vector = (double *)calloc(dimensions, sizeof(double));
        if(!vector) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
        for(vec_idx=0;vec_idx<dimensions;vec_idx++)  {
            vector[vec_idx] = buff[buff_idx++];
        }
        inputData[line_idx] = vector;
    }

    free(buff);
    if (k > lineCounter){             /* The num of vectors should be >= k */
        printf("Invalid Input!\n");
        exit(1);
    }
    return inputData;
}


/** Calculates the Euclidean distance between two vectors **/
double euclideanDistance(const double* vec1, const double* vec2) {
    int i;
    double sum = 0, distance;
    for(i=0;i<dimensions;i++) {
        double tmp = vec1[i] - vec2[i];
        sum += pow(tmp,2);
    }
    distance = sqrt(sum);
    return distance;
}

/** Checks whether the distance between the old and new
 *  centroid is smaller than epsilon **/
int centroidsFar(double **old, double** new, int k) {
    int i;
    for(i=0; i<k; i++) {
        if(euclideanDistance(*old, *new) > 0.001) {
            return 1;
        }
        old += 1;
        new += 1;
    }
    return 0;
}

/** Copies the current list of centroids,
 *  before calculating the new next list. **/
void copyCentroids(double **target, double **source, int k) {
    int i, j;
    for(i=0; i<k; i++) {
        for(j=0; j<dimensions; j++) {
            target[i][j] = source[i][j];
        }
    }
}

/** Creates the list of clusters for the first time.
 *  The centroids will be the first k vectors  **/
double** createCentroid(int k) {
    int i, j;
    double inf = INFINITY;
    double **new_centroids = (double **) calloc(k, sizeof(double *));
    if(!new_centroids) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < k; i++) {
        double *centroid = (double *)calloc(dimensions, sizeof(double));
        if(!centroid) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
        for (j = 0; j < dimensions; j++) {
            centroid[j] = inf;
        }
        new_centroids[i] = centroid;
    }
    return new_centroids;
}

/** Summarize all vectors in each cluster into specific vector **/
void updateSumCentroid(double* vec1, const double* vec2) {
    int i;
    for (i=0; i<dimensions;i++) {
        vec1[i] += vec2[i];
    }
}

/** Resets all centroids to values 0.0 **/
void resetCentroids(double** centroids, int k){
    int i,j;
    for (i=0;i<k;i++) {
        for(j=0;j<dimensions;j++) {
            centroids[i][j] = 0.0;
        }
    }
}

/** Calculates the new centroids values by dividing
 * the sum of vectors in each cluster divided by the
 * number of vectors in it. **/
void updateCentroids(double **centroids, const int arr[], int k) {
    int i, j;
    for(i=0;i<k;i++) {
        for(j=0;j<dimensions;j++) {
            centroids[i][j] = centroids[i][j] / arr[i];
        }
    }
}

/** Resets the number of vectors in each cluster after
 * each iteration, while optimizing the centroid **/
void resetSizes(int* arr, int k) {
    int i;
    for(i=0;i<k;i++) {
        arr[i] = 0;
    }
}

/** Creates the output file with the list of clusters that
 *  each contains the relevant vectors. **/
void writeFile(double **centroids, int k, char* outputName) {
    FILE *fp = NULL;
    int i,j;
    fp = fopen(outputName, "w");
    if (fp == NULL) {
        printf("An Error Has Occurred!");
    }
    for(i=0;i<k;i++) {
        for(j=0;j<dimensions-1;j++) {
            fprintf(fp, "%0.4lf,", centroids[i][j]);
        }
        fprintf(fp,"%0.4lf\n",centroids[i][j]);
    }
    fclose(fp);
}

/** The main function
 * Creates K clusters, which are a division of all
 * vectors into K groups so that each vector belongs
 * to the cluster with its nearest centroid. **/
int kMean(int k, int max_iter,char *fileName, char *outputName) {
    /* convert data in File to matrix of points.*/
    double **inputData, **centroids, **new_centroids;
    int i, j, counter=0, *sizes;

    inputData = readFile(k, fileName);

    /* Initialize clusters centers.*/
    centroids = (double **)calloc(k, sizeof(double *));
    if(!centroids) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    for(i=0; i < k; i++) {
        double *centroid = (double*) calloc(dimensions,sizeof(double));
        if(!centroid) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
        for(j=0; j<dimensions; j++) {
            centroid[j] = inputData[i][j];
        }
        centroids[i] = centroid;
    }
    new_centroids = createCentroid(k);
    sizes = (int *) calloc(k, sizeof(int));
    if(!sizes) {
        printf("An Error Has Occurred\n");
        exit(1);
    }

    while(counter<max_iter && centroidsFar(centroids, new_centroids, k)) {
        int cluster_idx;
        if(counter != 0) {
            copyCentroids(centroids, new_centroids, k);
        }
        resetCentroids(new_centroids,k);
        resetSizes(sizes,k);

        for(i=0; i<lineCounter;i++) {
            double min_euc = INFINITY;
            cluster_idx = k+1;
            for(j=0;j<k;j++) {
                double tmp = euclideanDistance(inputData[i],centroids[j]);
                if(tmp<min_euc) {
                    cluster_idx = j;
                    min_euc = tmp;
                }
            }

            updateSumCentroid(new_centroids[cluster_idx], inputData[i]);
            sizes[cluster_idx] += 1;
        }
        updateCentroids(new_centroids, sizes, k);
        counter++;
    }
    free(inputData);
    free(new_centroids);
    writeFile(centroids, k, outputName);
    free(centroids);
    return 0;
}


int main(int argc, char* argv[]) {
    int k, max_iter = 200;
    char *fileName, *outputName;

    if (argc == 5){                   /*received max_iter:*/
        k = atoi(argv[1]);
        max_iter = atoi(argv[2]);
        fileName = argv[3];
        outputName = argv[4];
    }
    else {                            /*max_iter = default*/
        k = atoi(argv[1]);
        fileName = argv[2];
        outputName = argv[3];
    }

    if (k <= 0) {
        printf("Invalid Input!\n");
        exit(1);
    }
    kMean(k, max_iter, fileName, outputName);
}
