## This is my implementation of the KMeans algorithm, completed during my Computer Science degree at Tel Aviv University.

(https://serokell.io/files/q4/q49pm3tx.K-Means_Clustering_Algorithm_pic1_(1).png)

### This repo includes the following files:
1. spkmeans.py: Python interface of my code.
2. spkmeans.h: C header file.
3. spkmeans.c: C interface of my code.
4. spkmeansmodule.c: Python C API wrapper.
5. setup.py: The setup file.
6. comp.sh: my compilation script.


### How to run this code:
1. first run this command: 
* $python setup.py build_ext --inplace (try python3)
2. run:
* $bash comp.sh
3. After successful compilation, the program can be executed, for example:
* ./spkmeans lnorm input_for_example.txt

### Input and Output of the program: 
1. Reading user CMD arguments:
    - (a) **goal** (enum): Can get the following values:
        *   i. **wam:** Calculate and output the Weighted Adjacency Matrix.
        *   ii. **ddg:** Calculate and output the Diagonal Degree Matrix.
        *   iii. **lnorm:** Calculate and output the Normalized Graph Laplacian
        *   iv. **jacobi:** Calculate and output the eigenvalues and eigenvectors.
    - (b) **file name** (.txt or .csv): The path to the Input file, it will contain N data points for all above goals except Jacobi, in case the goal is Jacobi the input is a symmetric matrix, the file extension is .txt or .csv.

2. Outputting the following:
    * Case of ’Jacobi’: The first line will be the eigenvalues, second line onward will be the corresponding eigenvectors (printed as columns).
    * Else: output the required matrix separated by a comma, such that each row is in a line of its own.
