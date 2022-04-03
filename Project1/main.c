// Preprocessor statements
#include <stdio.h>
#include "math.h"
#include "stdlib.h"

// functions declaration
double ** readFile(char *fileName);
int kMean(int k, char *fileName, int max_iter);
double euclideanDistance(const double* vec1, const double* vec2);
int centroidsFar(double **old, double** new, int k);
void copyCentroids(double **target, double **source, int k);
double** createCentroid(double num, int k);
void updateSumCentroid(double* vec1, double* vec2);
void resetCentroids(double** centroids, int k);
void updateCentroids(double **centroids, const int arr[], int k);
void resetSizes(int* arr, int k);
void writeFile(double **centroids, int k);


// Global variables
int dimensions = 1;
int lineCounter = 0;

double ** readFile(char *fileName) {
    FILE *fp = NULL;

    // Get Number of vectors and dimensions in the file.
    fp = fopen(fileName, "r");
    if(fp==NULL) {
        printf("Invalid Input!\n");
        // Needs to terminate the program.
    }
    char ch;
    // int dimensions = 1;
    while((ch=fgetc(fp))!= EOF) {
        if (lineCounter == 0 && ch == ',')
            dimensions++;
        if(ch=='\n')
            lineCounter++;
    }
    fclose(fp);

    // Create list of pointers in the right size.
    fp = fopen(fileName, "r");
    double **inputData = (double **)calloc(lineCounter, sizeof(double *));
    // double *vector = (double *)calloc(dimensions, sizeof(double));
    int size_of_buff = dimensions*lineCounter;
    double *buff = calloc(size_of_buff, sizeof(double));
    int i=0;
    while(fscanf(fp, "%lf,", &buff[i++])!=EOF){
//        printf("%0.4lf ", buff[i-1] );
    }
    fclose(fp);

    // copy buff to my inputData matrix.
    int buff_idx=0;
    int vec_idx, line_idx;
    for(line_idx=0; line_idx<lineCounter; line_idx++) {
        double *vector = (double *)calloc(dimensions, sizeof(double));
        for(vec_idx=0;vec_idx<dimensions;vec_idx++)  {
            vector[vec_idx] = buff[buff_idx++];
        }
        inputData[line_idx] = vector;
    }
    free(buff);

    return inputData;
}


double euclideanDistance(const double* vec1, const double* vec2) {
    int i;
    double sum = 0;
    for(i=0;i<dimensions;i++) {
        double tmp = vec1[i] - vec2[i];
        sum += pow(tmp,2);
    }
    double distance = sqrt(sum);
    return distance;
}
int centroidsFar(double **old, double** new, int k) {
    int i;
    for(i=0; i<k; i++) {
        if(euclideanDistance(*old, *new) > 0.01) {
            return 1;
        }
        old += 1;
        new += 1;
    }
    return 0;
}

void copyCentroids(double **target, double **source, int k) {
    int i, j;
    for(i=0; i<k; i++) {
        for(j=0; j<dimensions; j++) {
            target[i][j] = source[i][j];
        }
    }
}

double** createCentroid(double num, int k) {
    int i, j;
    double inf = INFINITY;
    double **new_centroids = (double **) calloc(k, sizeof(double *));
    for (i = 0; i < k; i++) {
        double *centroid = (double *)calloc(dimensions, sizeof(double));
        for (j = 0; j < dimensions; j++) {
            centroid[j] = inf;
        }
        new_centroids[i] = centroid;
    }
    return new_centroids;
}

void updateSumCentroid(double* vec1, double* vec2) {
    int i;
    for (i=0; i<dimensions;i++) {
        vec1[i] += vec2[i];
    }
}

void resetCentroids(double** centroids, int k){
    int i,j;
    for (i=0;i<k;i++) {
        for(j=0;j<dimensions;j++) {
            centroids[i][j] = 0.0;
        }
    }
}
void updateCentroids(double **centroids, const int arr[], int k) {
    int i, j;
    for(i=0;i<k;i++) {
        for(j=0;j<dimensions;j++) {
            centroids[i][j] = centroids[i][j] / arr[i];
        }
    }
}

void resetSizes(int* arr, int k) {
    int i;
    for(i=0;i<k;i++) {
        arr[i] = 0;
    }
}

void writeFile(double **centroids, int k) {
    FILE *fp = NULL;
    fp = fopen("myCFile.txt", "w");
    if (fp == NULL) {
        printf("An Error Has Occurred!");
    }
    int i,j;
    for(i=0;i<k;i++) {
        for(j=0;j<dimensions-1;j++) {
            fprintf(fp, "%0.4lf,", centroids[i][j]);
        }
        fprintf(fp,"%0.4lf\n",centroids[i][j]);
    }
    fclose(fp);
    // Need to add EOF for the File.
}

int kMean(int k, char *fileName, int max_iter) {
    // convert data in File to matrix of points.
    double **inputData = readFile(fileName);

    // Initialize clusters centers.
    double **centroids = (double **)calloc(k, sizeof(double *));
    int i, j;
    for(i=0; i < k; i++) {
        double *centroid = (double*) calloc(dimensions,sizeof(double));
        for(j=0; j<dimensions; j++) {
            centroid[j] = inputData[i][j];
        }
        centroids[i] = centroid;
    }
    int counter=0;
    double inf = INFINITY;
    double **new_centroids = createCentroid(inf, k);
    int *sizes = (int *) calloc(k, sizeof(int));

    while(counter<max_iter && centroidsFar(centroids, new_centroids, k)) {
        if(counter != 0) {
            copyCentroids(centroids, new_centroids, k);
        }
        resetCentroids(new_centroids,k);
        int cluster_idx;
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
    writeFile(centroids, k);
    return 0;
    // printf("counter is: %d\n", counter);
//    for(i=0;i<k;i++) {
//        for(j=0; j<dimensions;j++) {
//            printf("%0.4lf,", centroids[i][j]);
//        }
//        printf("\n");
//    }

}




int main() {
    kMean(7,"input_2.txt", 200);
}