
#ifndef FINAL_CLEAN_2_SPKMEANS_H
#define FINAL_CLEAN_2_SPKMEANS_H

/* Global variables*/
extern int DIMS;
extern int NUM_OF_VECS;
extern double** GENERAL_MAT;
extern double** GENERAL_I_MAT;

typedef enum {wam, ddg, lnorm, jacobi,spk}GOAL;

typedef struct eigenTuple{
    double eigenValue;
    double* eigenVector;
}eigenTuple;

/* functions declaration */

void invalid_input(void);
void error_occurred(void);
void check_allocate(double** p);
void copy_mat(double **source, double **dest, int rows, int cols);
void reset_mat(double** mat, int rows, int cols);
void free_mat(double **mat, int rows);
void free_GENERAL(double **vectors);
GOAL set_GOAL(char *goal);
void validateArg(char* arg_str);
double ** readFile(char *fileName);
void print_mat(double** mat, int rows,int cols);
void print_eigenValues(double **mat);
double euclideanDistance(const double* vec1, const double* vec2);
double Weighted_distance(const double* vec1, const double* vec2);
double** create_mat(int rows, int cols);
double** create_W(double **vectors);
double calc_deg(double* row);
double** wTOd(double **W);
double** Dsqrt(double **D);
double** create_I_mat(int MAT_SIZE);
void mult_mat(double** Mat1, double **Mat2, double** dest_res);
void subtract_mat(double **mat1, double **mat2, int MAT_SIZE,double** dest);
double **create_Lnorm(double** sqrtD, double** W);
double** create_transpose(double **mat, int rows, int cols,double** P_transpose);
double off_mat(double **mat, int MAT_SIZE);
int check_if_diag(double **A,double **A_tag, double eps, int MAT_SIZE);
void max_abs_value(double** mat,int *row, int *col, int MAX_SIZE);
double calc_theta(double **L_norm, int max_val_ROW, int max_val_COL);
double calc_t(double theta);
double calc_c(double t);
double** create_P(double **A, int MAT_SIZE,double** P);
void fill_eigenTuples(eigenTuple* eigenTuples, double **A_tag, double **V, int MAT_SIZE);
int compare_by_eigenValue(const void *p1, const void *p2);
int find_k(eigenTuple* eigenTuples,int len);
void mat_renormalization(double** U, int dims, int k);
double **create_U(eigenTuple* eigenTuples, int k, int dims);
void Jacobi_algorithm(double*** A, double*** V);
double** create_T(double **A, double **V, int MAT_SIZE, int *k);
double** jacobi_val_for_py(double **A, double **V, int MAT_SIZE);

#endif
