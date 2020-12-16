#include <omp.h> 
#include <iostream>
#include <string>
#include <climits>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

class Matrix {
    private:
        double** matrix;
        int dimension;
    
    public:
        Matrix(int theDimension, double** theMatrix) {
            dimension = theDimension;
            matrix = theMatrix;
        }

        Matrix(int inputDimension, int inputChoice) {
            dimension = inputDimension;            
            matrix = (double**) malloc(dimension * sizeof(double*));
            for (int i = 0; i < dimension; i++) {
                matrix[i] = (double*) malloc(dimension * sizeof(double));
            }

            switch (inputChoice) {
                case 1:
                    
                    printf("Enter %d numbers for the matrix:\n", (dimension * dimension));
                    for (int i = 0; i < dimension; i++) {
                        for (int j = 0; j < dimension; j++) {
                            int input = 0;
                            scanf("%d", &input);
                            matrix[i][j] = (double)input;
                        }
                    }
                    printf("Done.\n");

                    printf("The matrix is: \n");
                    for (int i = 0; i < dimension; i++) {
                        for (int j = 0; j < dimension; j++) {
                            printf("%f ", matrix[i][j]);
                        }
                        printf("\n");
                    }
                    printf("\n");

                    break;
                case 2:
                
                    printf("Generating numbers... ");
                    for (int i = 0; i < dimension; i++) {
                        for (int j = 0; j < dimension; j++) {
                            int randInt = rand() % 100 + 1; 
                            // double f = (double)rand() / RAND_MAX;
                            // matrix[i][j] = 1.0 + f * (100.0 - 1.0);
                            matrix[i][j] = (double) randInt;
                            // printf("Random int: %f\n", (double) randInt);
                        }
                    }
                    printf("Done.\n\n");
                    break;
                default:
                    printf("Invalid response. Aborting program...\n");
                    exit(0);
                    break;
            }
        }

        //row parallelization?
        //parallelize the loop that starts on line 78? maybe only outer loop
        Matrix multiply(Matrix otherArray) {
            int i, j;
            double** matrixProduct;
            matrixProduct = (double**)malloc(dimension * sizeof(double*));
            for (int i = 0; i < dimension; i++) {
                matrixProduct[i] = (double*)malloc(dimension * sizeof(double));
            }

            #pragma omp parallel for private(j)
            for (i = 0; i < dimension; i++) {
                for (j = 0; j < dimension; j++) {
                    double tempSum = 0;
                    for (int k = 0; k < dimension; k++) {
                        tempSum += matrix[i][k] * otherArray.getMatrix()[k][j];
                        // printf("tempSum = %f\n", tempSum);
                    }                        
                    // printf("final tempSum = %f\n", tempSum);
                    matrixProduct[i][j] = tempSum;
                }
            }

            Matrix m(dimension, matrixProduct);
            return m;
        }

        //row parallelization? and how to implement blocked algorithm?
        void luDecomposition() {

            //print matrix
            // printf("Matrix before operation: \n");
            // for (int i = 0; i < dimension; i++) {
            //     for (int j = 0; j < dimension; j++) {
            //         printf("%f ", matrix[i][j]);
            //     }
            //     printf("\n");
            // }
            // printf("\n");

            // printf("test1\n");
            int* p = (int*) malloc(dimension * sizeof(int));
            double** u = (double**) malloc(dimension * sizeof(double*));
            for (int i = 0; i < dimension; i++) {
                u[i] = (double*) malloc(dimension * sizeof(double));
            }            
            double** l = (double**) malloc(dimension * sizeof(double*));
            for (int i = 0; i < dimension; i++) {
                l[i] = (double*) malloc(dimension * sizeof(double));
            }
            
            // printf("test2\n");
            
            int uJ;
            #pragma omp parallel for private(uJ)
            for (int i = 0; i < dimension; i++) {
                for (uJ = 0; uJ < dimension; uJ++) {
                    if (uJ < i) {
                        //should be a 0
                        u[i][uJ] = 0.0;
                    }
                    else {
                        u[i][uJ] = matrix[i][uJ];
                    }
                }
            }
            // printf("test3\n");

            int lJ;
            #pragma omp parallel for private(lJ)
            for (int i = 0; i < dimension; i++) {
                for (lJ = 0; lJ < dimension; lJ++) {
                    if (lJ == i) {
                        l[i][lJ] = 1.0;
                    }
                    else if (lJ > i) {
                        l[i][lJ] = 0.0;
                    }
                    else {
                        l[i][lJ] = matrix[i][lJ];
                    }
                }
            }
            
            // printf("test4\n");

            for (int i = 0; i < dimension; i++) {
                p[i] = i;
            }
            
            // printf("test5\n");

            int k_prime = 0;
            for (int k = 0; k < dimension; k++) {
                double max = 0.0;
                for (int i = k; i < dimension; i++) {
                    if (max < std::abs(matrix[i][k])) {
                        max = std::abs(matrix[i][k]);
                        k_prime = i;
                    }
                }
                if (max = 0) {
                    printf("Error: singular matrix.\n");
                    exit(0);
                }
                //swap pi[k] and pi[k_prime]
                int temp_p = p[k];
                p[k] = p[k_prime];
                p[k_prime] = temp_p;

                //swap matrix[k] and matrix[k_prime]
                double* temp_matrix = matrix[k];
                matrix[k] = matrix[k_prime];
                matrix[k_prime] = temp_matrix;

                //swap l[k][0 to k-1] and l[k_prime][0 to k-1]
                for (int i = 0; i <= k-1; i++) {
                    double temp_l = l[k][i];
                    l[k][i] = l[k_prime][i];
                    l[k_prime][i] = temp_l;
                }


                u[k][k] = matrix[k][k];
                for (int i = k+1; i < dimension; i++) {
                    l[i][k] = matrix[i][k]/u[k][k];
                    u[k][i] = matrix[k][i];
                }

                int lastJ;
                #pragma omp parallel for private(lastJ)
                for (int i = k+1; i < dimension; i++) {
                    for (lastJ = k+1; lastJ < dimension; lastJ++) {
                        matrix[i][lastJ] = matrix[i][lastJ] - (l[i][k]*u[k][lastJ]);
                    }
                }
            }

            //print matrix
            // printf("Elements of the original matrix:\n");
            // for (int i = 0; i < dimension; i++) {
            //     for (int j = 0; j < dimension; j++) {
            //         printf("%f ", matrix[i][j]);
            //     }
            //     printf("\n");
            // }
            // printf("\n");

            // //print p
            // printf("Elements of p:\n");
            // for (int i = 0; i < dimension; i++) {
            //     printf("%d ", p[i]);
            // }
            // printf("\n");

            // //print l
            // printf("Elements of l:\n");
            // for (int i = 0; i < dimension; i++) {
            //     for (int j = 0; j < dimension; j++) {
            //         printf("%f ", l[i][j]);
            //     }
            //     printf("\n");
            // }
            // printf("\n");

            // //print u
            // printf("Elements of u:\n");
            // for (int i = 0; i < dimension; i++) {
            //     for (int j = 0; j < dimension; j++) {
            //         printf("%f ", u[i][j]);
            //     }
            //     printf("\n");
            // }
            // printf("\n");
            
        }
        
        double** getMatrix() {
            return matrix;
        }

        void printMatrix() {
            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    printf("Element matrix[%d][%d] is %f.\n", i, j, matrix[i][j]);
                }
            }
        }
};

int main() {
     
    int opChoice = 0;
    printf("Choose: Enter 1 for matrix-multiply, or 2 for LU decomposition\n");
    scanf("%d", &opChoice);

    switch (opChoice) {
        case 1: {
            int multDimens = 0;
            int multChoice = 0;
            printf("Enter the dimension of matrix:\n");
            scanf("%d", &multDimens);
            printf("Choose option:\n Enter 1 for manual input.\n Enter 2 for automatic random input.\n");
            scanf("%d", &multChoice);
            // printf("First matrix is:\n");
            Matrix sm(multDimens, multChoice);
            // sm.printMatrix();
            // printf("Second matrix is:\n");
            Matrix sm2(multDimens, multChoice);
            // sm2.printMatrix();
            printf("Product is:\n");
            // sm.multiply(sm2).printMatrix();
            sm.multiply(sm2);
            break;
        }
            

        case 2: {
            int luDimens = 0;
            int luChoice = 0;
            printf("Enter the dimension of matrix:\n");
            scanf("%d", &luDimens);
            printf("Choose option:\n Enter 1 for manual input.\n Enter 2 for automatic random input.\n");
            scanf("%d", &luChoice);
            Matrix m(luDimens, luChoice);

            m.luDecomposition();
            break;
        }
            

        default: {
            printf("Invalid input. Aborting....\n");
            exit(0);
            break;
        }
            
    }

    return 0;
}