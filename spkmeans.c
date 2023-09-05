#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spkmeans.h"

int DIMS = 1;
int NUM_OF_VECS = 0;
double** GENERAL_MAT;
double** GENERAL_I_MAT;

void invalid_input(void){
    printf("Invalid Input!\n");
    exit(1);
}

void error_occurred(void){
    printf("An Error Has Occurred\n");
    exit(1);
}

void reset_mat(double** mat, int rows, int cols){
    int i,j;
    for(i=0; i<rows; i++){
        for(j=0; j<cols; j++){
            mat[i][j] = 0;
        }
    }
}

void reset_I_mat(double** mat, int dims){
    int i;
    reset_mat(mat,dims, dims);
    for(i=0; i<dims; i++){
        mat[i][i] = 1;
    }
}

void copy_mat(double **source, double **dest, int rows, int cols){
    int i,j;
    for(i=0; i<rows; i++){
        for(j=0; j<cols; j++){
            dest[i][j] = source[i][j];
        }
    }
}

void free_mat(double **mat, int rows){
    int i;
    for (i=0; i<rows;i++){
        free(mat[i]);
    }
    free(mat);
}

void free_GENERAL(double **vectors){
    free_mat(vectors,NUM_OF_VECS);
    free_mat(GENERAL_MAT,NUM_OF_VECS);
    free_mat(GENERAL_I_MAT,NUM_OF_VECS);
}

GOAL set_GOAL(char *goal) {
    int i;
    char *goals[] = {"wam", "ddg", "lnorm", "jacobi"};
    for (i = 0; i < 4; i++) {
        if (strcmp(goals[i], goal) == 0) return i;
    }
    invalid_input();
    return wam;
}

void validateArg(char* arg_str) {
    int i = 0;
    char *filename, *extension, *arg_str_cpy;
    while (arg_str[i] != '\0') {
        i++;
    }
    if (strcmp(&arg_str[i - 1], ".") == 0) {
        invalid_input();
    }
    arg_str_cpy = (char *) calloc(strlen(arg_str)+1, sizeof(char));
    if(!arg_str_cpy) invalid_input();
    strcpy(arg_str_cpy, arg_str);
    filename = strtok(arg_str_cpy, ".");
    while (filename != NULL) {
        extension = filename;
        filename = strtok(NULL, ".");
    }
    if ((strcmp(extension, "txt") != 0) && (strcmp(extension, "csv") != 0)) {
        invalid_input();
    }
    free(arg_str_cpy);
}


/** The function reads all the vectors from the input file.
 *  Counts the number of lines of text as the number of vectors.
 *  Counts the number of values in each vector as the number of dimension **/
double ** readFile(char *fileName) {
    FILE *fp = NULL;
    char ch;
    double **inputData, *buff;
    int size_of_buff, i=0, buff_idx=0, vec_idx, line_idx;
    DIMS = 1;
    NUM_OF_VECS = 0;
    /* Get Number of vectors and DIMS in the file.*/
    fp = fopen(fileName, "r");
    if(!fp) invalid_input();
    while((ch=fgetc(fp))!= EOF) {
        if (NUM_OF_VECS == 0 && ch == ',')
            DIMS++;
        if(ch=='\n')
            NUM_OF_VECS++;
    }
    fclose(fp);

    /* Create list of pointers in the right size.*/
    fp = fopen(fileName, "r");
    inputData = (double **)calloc(NUM_OF_VECS, sizeof(double *));
    if(!inputData) error_occurred();
    size_of_buff = DIMS * NUM_OF_VECS;
    buff = calloc(size_of_buff, sizeof(double));
    if(!buff) error_occurred();
    while(fscanf(fp, "%lf,", &buff[i++])!=EOF){
    }
    fclose(fp);

    /* copy buff to my inputData matrix. */
    for(line_idx=0; line_idx < NUM_OF_VECS; line_idx++) {
        double *vector = (double *)calloc(DIMS, sizeof(double));
        if(!vector) error_occurred();
        for(vec_idx=0; vec_idx < DIMS; vec_idx++)  {
            vector[vec_idx] = buff[buff_idx++];
        }
        inputData[line_idx] = vector;
    }

    free(buff);
    return inputData;
}


void print_mat(double** mat, int rows,int cols){
    int i,j;
    for(i = 0; i<rows; i++){
        for(j=0; j<cols-1; j++){
            printf("%.4f,", mat[i][j]);
        }
        printf("%.4f", mat[i][j]);
        printf("\n");
    }
}

void print_eigenValues(double **mat){
    int i;
    for(i=0; i<(NUM_OF_VECS-1); i++){
        printf("%.4f,", mat[i][i]);
    }
    printf("%.4f\n", mat[i][i]);
}

/** Calculates the Euclidean distance between two vectors **/
double euclideanDistance(const double* vec1, const double* vec2) {
    int i;
    double sum = 0, distance;
    for(i=0; i < DIMS; i++) {
        double tmp = vec1[i] - vec2[i];
        sum += pow(tmp,2);
    }
    distance = sqrt(sum);
    return distance;
}

double Weighted_distance(const double* vec1, const double* vec2){
    double euc_dist = euclideanDistance(vec1, vec2);
    double weighted_dist = exp(-euc_dist/2);
    return weighted_dist;
}

double** create_mat(int rows, int cols){
    int i;
    double** mat = (double **) calloc(rows, sizeof (double *));
    for(i=0; i < rows; i++){
        mat[i] = (double *) calloc(cols, sizeof (double*));
    }
    return mat;
}

double** create_W(double **vectors){
    int i,j;
    double dist;
    double** W = create_mat(NUM_OF_VECS, NUM_OF_VECS);
    for(i=0; i < NUM_OF_VECS; i++){
        for(j=i+1; j < NUM_OF_VECS; j++){
            dist = Weighted_distance(vectors[i], vectors[j]);
            W[i][j] = W[j][i] = dist;
        }
    }
    return W;
}

double calc_deg(double* row){
    int i;
    double sum = 0;
    for(i=0; i<NUM_OF_VECS; i++){
        sum += row[i];
    }
    return sum;
}

double** wTOd(double **W){
    int i;
    double **D = create_mat(NUM_OF_VECS, NUM_OF_VECS);
    for(i=0; i < NUM_OF_VECS; i++) {
        D[i][i] = calc_deg(W[i]);
    }
    return D;
}

double** Dsqrt(double **D){
    int i;
    for(i=0; i < NUM_OF_VECS; i++){
        if(D[i][i] != 0){
            D[i][i] = (1/(sqrt(D[i][i])));
        }
    }
    return D;
}

double** create_I_mat(int MAT_SIZE){
    int i;
    double** I = create_mat(MAT_SIZE, MAT_SIZE);
    for(i=0; i<MAT_SIZE; i++){
        I[i][i] = 1;
    }
    return I;
}

void mult_mat(double** Mat1, double **Mat2, double** dest_res){
    int i,j,k;
    reset_mat(GENERAL_MAT,NUM_OF_VECS, NUM_OF_VECS);
    for(i=0; i < NUM_OF_VECS; i++){
        for(j=0; j < NUM_OF_VECS; j++){
            for(k=0; k<NUM_OF_VECS; k++){
                GENERAL_MAT[i][j] += Mat1[i][k] * Mat2[k][j];
            }
        }
    }
    copy_mat(GENERAL_MAT,dest_res,NUM_OF_VECS,NUM_OF_VECS);
}

void subtract_mat(double **mat1, double **mat2, int MAT_SIZE,double** dest){
    int i,j;
    reset_mat(GENERAL_MAT,MAT_SIZE,MAT_SIZE);
    for(i=0; i<MAT_SIZE; i++){
        for(j=0; j<MAT_SIZE; j++){
            GENERAL_MAT[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    copy_mat(GENERAL_MAT, dest,MAT_SIZE,MAT_SIZE);
}

double **create_Lnorm(double** sqrtD, double** W){
    double **Lnorm = create_mat(NUM_OF_VECS,NUM_OF_VECS);
    reset_I_mat(GENERAL_I_MAT,NUM_OF_VECS);
    mult_mat(sqrtD, W,Lnorm);
    mult_mat(Lnorm, sqrtD,Lnorm);
    subtract_mat(GENERAL_I_MAT, Lnorm, NUM_OF_VECS,Lnorm);
    return Lnorm;
}

double** create_transpose(double **mat, int rows, int cols,double** P_transpose){
    int i,j;
    for(i=0; i < cols; i++){
        for(j=0; j < rows; j++){
            P_transpose[i][j] = mat[j][i];
        }
    }
    return P_transpose;
}

double off_mat(double **mat, int MAT_SIZE){
    int i,j;
    double sum = 0;
    for(i=0; i<MAT_SIZE; i++){
        for(j=0; j<MAT_SIZE; j++){
            if (i != j){
                sum += pow(mat[i][j],2);
            }
        }
    }
    return sum;
}

/* I independently added abs for the expression! */
int check_if_diag(double **A,double **A_tag, double eps, int MAT_SIZE){
    if (fabs(off_mat(A,MAT_SIZE) - off_mat(A_tag,MAT_SIZE)) <= eps) return 1;
    else return 0;
}

/* check ONLY the top triangle of the matrix */
void max_abs_value(double** mat,int *row, int *col, int MAT_SIZE){
    int i,j;
    double maxi = -1;
    for(i=0; i<MAT_SIZE; i++){
        for(j=i+1; j<MAT_SIZE; j++){
            if (fabs(mat[i][j]) > maxi){
                *row = i, *col = j;
                maxi = fabs(mat[i][j]);
            }
        }
    }
}

double calc_theta(double **L_norm, int max_val_ROW, int max_val_COL){
    double theta;
    theta = ((L_norm[max_val_COL][max_val_COL]- L_norm[max_val_ROW][max_val_ROW])/
             (2*L_norm[max_val_ROW][max_val_COL]));
    return theta;
}

double calc_t(double theta){
    int sign = (theta < 0 ? -1: 1);
    double t = (sign / (fabs(theta) + sqrt(pow(theta,2)+1)));
    return t;
}

double calc_c(double t){
    double c = (1 / sqrt(pow(t,2)+1));
    return c;
}

/* It is required to validate before that MAT_SIZE >= 2 for the symmetric attribute */
double** create_P(double **A, int MAT_SIZE,double** P){
    int i,j;
    double c,t,s,theta;
    reset_I_mat(P,MAT_SIZE);
    max_abs_value(A, &i, &j, MAT_SIZE);
    theta = calc_theta(A, i, j);
    t = calc_t(theta);
    c = calc_c(t);
    s = t*c;
    P[i][i] = P[j][j] = c;
    P[i][j] = s;
    P[j][i] = -s;
    return P;
}


void fill_eigenTuples(eigenTuple* eigenTuples, double **A_tag, double **V, int MAT_SIZE){
    int i;
    double **V_transpose = create_mat(MAT_SIZE,MAT_SIZE);
    V_transpose = create_transpose(V, MAT_SIZE, MAT_SIZE,V_transpose);
    for(i=0; i<MAT_SIZE; i++){
        eigenTuples[i].eigenValue = A_tag[i][i];
        eigenTuples[i].eigenVector = V_transpose[i];
    }
}

int compare_by_eigenValue(const void *p1, const void *p2){
    const eigenTuple *q1 = p1; const eigenTuple *q2 = p2;
    if (fabs(q1->eigenValue - q2->eigenValue) <= (1.0 * pow(10, -5))) return 0;
    else if (q1->eigenValue < q2->eigenValue) return 1;
    else return -1;
}

/* what about [8,1,1] ? k == 0 ?? so I wrote k = i+1 */
int find_k(eigenTuple* eigenTuples,int len){
    int i,k=0;
    double maxGap = -1;
    for(i=0; i < (len/2); i++){
        if (maxGap < fabs(eigenTuples[i].eigenValue- eigenTuples[i+1].eigenValue)){
            maxGap = fabs(eigenTuples[i].eigenValue- eigenTuples[i+1].eigenValue);
            k = i+1;
        }
    }
    return k;
}


void mat_renormalization(double** U, int dims, int k){
    int i,j;
    double denominator;
    for(i=0; i<dims; i++){
        denominator = 0;
        for(j=0; j<k; j++){
            denominator += pow(U[i][j],2);
        }
        denominator = sqrt(denominator);
        for(j=0; j<k; j++){
            if (denominator != 0) {
                U[i][j] = (U[i][j] / denominator);
            }
        }
    }
}

double **create_U(eigenTuple* eigenTuples, int k, int dims){
    int i;
    double** U = (double **) calloc(k, sizeof (double *));
    double** U_transpose = create_mat(dims,k);
    for(i=0; i<k; i++){
        U[i] = eigenTuples[i].eigenVector;
    }
    U_transpose = create_transpose(U, k, dims,U_transpose);
    free(U);
    return U_transpose;
}


void Jacobi_algorithm(double*** A, double*** V) {
    int iterations = 0;
    double eps = 1.0 * pow(10, -5);
    double **P=NULL, **P_transpose=NULL, **A_tag=NULL;
    A_tag = create_mat(NUM_OF_VECS, NUM_OF_VECS);
    *V = create_I_mat(NUM_OF_VECS);
    P = create_I_mat(NUM_OF_VECS);
    P_transpose = create_mat(NUM_OF_VECS,NUM_OF_VECS);
    while (iterations < 100) {
        P = create_P(*A, NUM_OF_VECS,P);
        P_transpose = create_transpose(P, NUM_OF_VECS, NUM_OF_VECS,P_transpose);
        mult_mat(P_transpose, *A,A_tag);
        mult_mat(A_tag, P,A_tag);
        mult_mat(*V, P,*V);
        if (check_if_diag(*A, A_tag, eps, NUM_OF_VECS) == 1) break;
        else {
            copy_mat(A_tag,*A,NUM_OF_VECS,NUM_OF_VECS);
            iterations++;
        }
    }
    *A = A_tag;
    free_mat(P_transpose,NUM_OF_VECS);
    free_mat(P,NUM_OF_VECS);
}

double** create_T(double **A, double **V, int MAT_SIZE, int *k) {
    double **U;
    eigenTuple *eigenTuples = (eigenTuple*)calloc(MAT_SIZE,sizeof(eigenTuples));
    fill_eigenTuples(eigenTuples, A, V, MAT_SIZE);
    qsort(eigenTuples, MAT_SIZE, sizeof(eigenTuple), compare_by_eigenValue);
    if (*k == 0) *k = find_k(eigenTuples, MAT_SIZE);
    U = create_U(eigenTuples, *k, MAT_SIZE);
    mat_renormalization(U,MAT_SIZE,*k);
    return U;  /* T = re-normal U */
}

double** jacobi_val_for_py(double **A, double **V, int MAT_SIZE){
    int i,j;
    double **res = create_mat(MAT_SIZE+1, MAT_SIZE);
    for(j=0; j<MAT_SIZE; j++){
        res[0][j] = A[j][j];
    }
    for(i=1; i<(MAT_SIZE+1); i++){
        for(j=0; j<MAT_SIZE; j++){
            res[i][j] = V[i-1][j];
        }
    }
    return res;
}

int main(int argc, char* argv[]){
    GOAL goal;
    char *fileName;
    double **vectors, **W=NULL, **D=NULL,**V=NULL,**A=NULL, **sqrtD=NULL, **Lnorm=NULL;
    if (argc != 3) invalid_input();
    goal = set_GOAL(argv[1]);
    fileName = argv[2];
    validateArg(fileName);
    vectors = readFile(fileName);
    GENERAL_MAT = create_mat(NUM_OF_VECS, NUM_OF_VECS);
    GENERAL_I_MAT = create_I_mat(NUM_OF_VECS);
    switch(goal){
        case wam:
            W = create_W(vectors);
            print_mat(W, NUM_OF_VECS, NUM_OF_VECS);
            free_mat(W,NUM_OF_VECS);
            break;
        case ddg:
            W = create_W(vectors);
            D = wTOd(W);
            print_mat(D, NUM_OF_VECS, NUM_OF_VECS);
            free_mat(W,NUM_OF_VECS);
            free_mat(D,NUM_OF_VECS);
            break;
        case lnorm:
            W = create_W(vectors);
            D = wTOd(W);
            sqrtD = Dsqrt(D);
            Lnorm = create_Lnorm(sqrtD, W);
            print_mat(Lnorm, NUM_OF_VECS, NUM_OF_VECS);
            free_mat(W,NUM_OF_VECS);
            free_mat(D,NUM_OF_VECS);
            free_mat(Lnorm,NUM_OF_VECS);
            break;
        case jacobi:
            A = vectors;
            if (NUM_OF_VECS != DIMS){
                free_GENERAL(vectors);
                invalid_input();
            }
            Jacobi_algorithm(&A,&V);
            print_eigenValues(A);
            print_mat(V,NUM_OF_VECS,NUM_OF_VECS);
            free_mat(V,NUM_OF_VECS);
            free_mat(A,NUM_OF_VECS);
            break;
        default:
            error_occurred();
    }
    free_GENERAL(vectors);
    exit(0);
}
