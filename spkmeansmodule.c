/*Preprocessor statements*/
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "string.h"
#include "ctype.h"
#include "spkmeans.h"
#include <stdio.h>
#include "math.h"
#include "stdlib.h"


/* functions declaration*/
GOAL set_GOAL_Python_vers(char *goal);
static void mat_to_py(PyObject *res,double ** mat, int rows, int cols);
int centroidsFar(double **old, double** new, int k, double eps);
void copyCentroids(double **target, double **source, int k);
double** createCentroid(int k);
void updateSumCentroid(double* vec1, const double* vec2);
void resetCentroids(double** centroids, int k);
void updateCentroids(double **centroids, const int arr[], int k);
void resetSizes(int* arr, int k);
int* lst_to_array(PyObject *indices, int k);
void free_GENERAL(double **vectors);
static PyObject* kMean(int k, int max_iter,char *fileName, PyObject *initCentroids, double eps);
static PyObject *fit_by_goal(PyObject *self, PyObject *args);
static PyObject* fit_kmeanspp(PyObject *self, PyObject *args);


GOAL set_GOAL_Python_vers(char *goal) {
    int i;
    char *goals[] = {"wam", "ddg", "lnorm", "jacobi","spk"};
    for (i = 0; i < 5; i++) {
        if (strcmp(goals[i], goal) == 0) return i;
    }
    invalid_input();
    return wam;
}

static void mat_to_py(PyObject *res,double ** mat, int rows, int cols){
    int i,j;
    PyObject *vec;
    for (i=0; i<rows; i++){
        vec = PyList_New(cols);
        for (j=0; j<cols; j++){
            PyList_SetItem(vec, j, PyFloat_FromDouble(mat[i][j]));
        }
        PyList_SetItem(res, i, vec);
    }
}


/** Checks whether the distance between the old and new
 *  centroid is smaller than epsilon **/
int centroidsFar(double **old, double** new, int k, double eps) {
    int i;
    for(i=0; i<k; i++) {
        if(euclideanDistance(*old, *new) > eps) {
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
        for(j=0; j<DIMS; j++) {
            target[i][j] = source[i][j];
        }
    }
}

/** Creates the list of clusters for the first time.
 *  The centroids will be the first k vectors  **/
double** createCentroid(int k) {
    int i, j;
    double inf = HUGE_VAL;
    double **new_centroids = (double **) calloc(k, sizeof(double *));
    if(!new_centroids) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < k; i++) {
        double *centroid = (double *)calloc(DIMS, sizeof(double));
        if(!centroid) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
        for (j = 0; j < DIMS; j++) {
            centroid[j] = inf;
        }
        new_centroids[i] = centroid;
    }
    return new_centroids;
}

/** Summarize all vectors in each cluster into specific vector **/
void updateSumCentroid(double* vec1, const double* vec2) {
    int i;
    for (i=0; i<DIMS;i++) {
        vec1[i] += vec2[i];
    }
}

/** Resets all centroids to values 0.0 **/
void resetCentroids(double** centroids, int k){
    int i,j;
    for (i=0;i<k;i++) {
        for(j=0;j<DIMS;j++) {
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
        for(j=0;j<DIMS;j++) {
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

int* lst_to_array(PyObject *indices, int k){
    int* arr = (int*) calloc(k, sizeof(int));
    if(!arr) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    int i;
    for(i=0; i<k; i++){
        PyObject *item = PyList_GetItem(indices, i);
        arr[i] = PyLong_AsLong(item);
    }
    return arr;
}

/** The main function
 * Creates K clusters, which are a division of all
 * vectors into K groups so that each vector belongs
 * to the cluster with its nearest centroid. **/
static PyObject* kMean(int k, int max_iter,char *fileName, PyObject *initIndices, double eps) {
    /* convert data in File to matrix of points.*/
    double **inputData, **centroids, **new_centroids;
    int i, j,d, counter=0, *sizes, index;
    int *indices;
    PyObject *res;
    inputData = readFile(fileName);
    /* Initialize clusters centers.*/
    centroids = (double **)calloc(k, sizeof(double *));
    if(!centroids) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    indices = lst_to_array(initIndices, k);
    for(i=0; i < k; i++) {
        double *centroid = (double*) calloc(DIMS,sizeof(double));
        if(!centroid) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
        index = indices[i];
        for(j=0; j<DIMS; j++) {
            centroid[j] = inputData[index][j];
        }
        centroids[i] = centroid;
    }
    free(indices);
    new_centroids = createCentroid(k);
    sizes = (int *) calloc(k, sizeof(int));
    if(!sizes) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    while(counter<max_iter && centroidsFar(centroids, new_centroids, k, eps)) {
        int cluster_idx;
        if(counter != 0) {
            copyCentroids(centroids, new_centroids, k);
        }
        resetCentroids(new_centroids,k);
        resetSizes(sizes,k);

        for(i=0; i<NUM_OF_VECS;i++) {
            double min_euc = HUGE_VAL;
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
    res = PyList_New(DIMS*k);
    d = 0;
    for (i = 0; i< k; i++){
        for (j = 0; j < DIMS; j++){
            PyList_Insert(res, d, PyFloat_FromDouble(centroids[i][j]));
            d++;
        }
    }
    for (i=0; i<k;i++){
        free(new_centroids[i]);
    }
    free(sizes);
    for (i=0; i<NUM_OF_VECS;i++){
        free(inputData[i]);
    }
    free(inputData);
    free(new_centroids);
    for (i=0; i<k;i++){
        free(centroids[i]);
    }
    free(centroids);
    return res;
}

static PyObject *fit_by_goal(PyObject *self, PyObject *args){
    GOAL goal;
    int k;
    char *fileName, *str_goal;
    double **vectors, **W=NULL, **D=NULL,**V=NULL,**A=NULL,**T=NULL, **sqrtD=NULL, **Lnorm=NULL, **AandV=NULL;
    PyObject *res_matrix;
    if (!PyArg_ParseTuple(args, "iss:fit_by_goal", &k,&str_goal, &fileName)) {
        return NULL;
    }
    goal = set_GOAL_Python_vers(str_goal);
    vectors = readFile(fileName);
    if (k > NUM_OF_VECS) invalid_input();
    GENERAL_MAT = create_mat(NUM_OF_VECS, NUM_OF_VECS);
    GENERAL_I_MAT = create_I_mat(NUM_OF_VECS);
    switch(goal){
        case spk:
            W = create_W(vectors);
            D = wTOd(W);
            sqrtD = Dsqrt(D);
            Lnorm = create_Lnorm(sqrtD,W);
            Jacobi_algorithm(&Lnorm,&V);
            T = create_T(Lnorm,V,NUM_OF_VECS,&k);
            res_matrix = PyList_New(NUM_OF_VECS);
            mat_to_py(res_matrix ,T,NUM_OF_VECS,k);
            free_mat(W,NUM_OF_VECS);
            free_mat(D,NUM_OF_VECS);
            free_mat(Lnorm,NUM_OF_VECS);
            free_mat(T,k);
            free_GENERAL(vectors);
            return res_matrix;
        case wam:
            W = create_W(vectors);
            res_matrix = PyList_New(NUM_OF_VECS);
            mat_to_py(res_matrix, W,NUM_OF_VECS,NUM_OF_VECS);
            free_mat(W,NUM_OF_VECS);
            free_GENERAL(vectors);
            return res_matrix;
        case ddg:
            W = create_W(vectors);
            D = wTOd(W);
            res_matrix = PyList_New(NUM_OF_VECS);
            mat_to_py(res_matrix ,D,NUM_OF_VECS,NUM_OF_VECS);
            free_mat(W,NUM_OF_VECS);
            free_mat(D,NUM_OF_VECS);
            free_GENERAL(vectors);
            return res_matrix;
        case lnorm:
            W = create_W(vectors);
            D = wTOd(W);
            sqrtD = Dsqrt(D);
            Lnorm = create_Lnorm(sqrtD, W);
            res_matrix = PyList_New(NUM_OF_VECS);
            mat_to_py(res_matrix,Lnorm,NUM_OF_VECS,NUM_OF_VECS);
            free_mat(W,NUM_OF_VECS);
            free_mat(D,NUM_OF_VECS);
            free_mat(Lnorm,NUM_OF_VECS);
            free_GENERAL(vectors);
            return res_matrix;
        case jacobi:
            A = vectors;
            if (NUM_OF_VECS != DIMS){
                free_GENERAL(vectors);
                invalid_input();
            }
            Jacobi_algorithm(&A,&V);
            AandV = jacobi_val_for_py(A,V, NUM_OF_VECS);
            res_matrix = PyList_New(NUM_OF_VECS+1);
            mat_to_py(res_matrix,AandV,NUM_OF_VECS+1,NUM_OF_VECS);
            free_mat(AandV,NUM_OF_VECS);
            free_mat(V,NUM_OF_VECS);
            free_mat(A,NUM_OF_VECS);
            free_GENERAL(vectors);
            return res_matrix;
        default:
            error_occurred();
    }
    exit(1);
}


static PyObject* fit_kmeanspp(PyObject *self, PyObject *args) {
    int max_iter, k;
    double eps;
    char *filename;
    PyObject *init_centroids;
    if (!PyArg_ParseTuple(args, "iidsO:fit", &k, &max_iter, &eps, &filename,
                          &init_centroids)) {
        return NULL;
    }
    return Py_BuildValue("O", kMean(k, max_iter, filename, init_centroids, eps));
}


static PyMethodDef kmeans_methods[] ={
        {"fit_by_goal", (PyCFunction)fit_by_goal, METH_VARARGS,PyDoc_STR("Kmeans_by_goal c python API")},
        {"fit_kmeanspp", (PyCFunction)fit_kmeanspp, METH_VARARGS,PyDoc_STR("Kmeans_pp c python API")},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeansspf",
        NULL,
        -1,
        kmeans_methods};

PyMODINIT_FUNC
PyInit_mykmeansspf(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m){
        return NULL;
    }
    return m;
}
