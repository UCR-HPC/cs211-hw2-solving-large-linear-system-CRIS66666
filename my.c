#ifndef __MY_C__
#define __MY_C__

#include "include.h"

int mydgetrf(double *A,int *ipiv,int n)
{
    //TODO
    //The return value (an integer) can be 0 or 1
    //If 0, the matrix is irreducible and the result will be ignored
    //If 1, the result is valid
    int i, j, k, t, maxind;
    double  max, temps, tempv;

    for(i=0;i<n-1;i++){
        maxind=i;
        max=fabs(A[i*n+i]);
        for(t=i+1;t<n;t++){
            if(fabs(A[t*n+i])>max){
                maxind=t;
                max=fabs(A[t*n+i]);
            }
        }
        if(fabs(max)==0.0){
            printf("error here.\n");
            return 0;
        }
        else{
            if(maxind!=i){
                // swap pivot indices
                temps=ipiv[i];
                ipiv[i]=ipiv[maxind];
                ipiv[maxind]=temps;
                // swap rows
                for(k=0;k<n;k++){
                    tempv=A[i*n+k];
                    A[i*n+k]=A[maxind*n+k];
                    A[maxind*n+k]=tempv;
                }
            }
        }
        for(j=i+1;j<n;j++){
            A[j*n+i] /=A[i*n+i];
            for(k=i+1;k<n;k++){
                A[j*n+k]-=A[j*n+i]*A[i*n+k];
            }
        }
    }
    return 1;
}

void mydtrsv(char UPLO,double *A,double *B,int n,int *ipiv)
{
    double *y = (double *)malloc(n * sizeof(double)); // for forward substitution
    double *x = (double *)malloc(n * sizeof(double)); // for backward substitution

    if (UPLO == 'L') {
        // Forward substitution for lower triangular matrix
        y[0] = B[ipiv[0]];
        for (int i = 1; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += y[j] * A[i * n + j];  // A(i, 1:i-1) in row-major
            }
            y[i] = B[ipiv[i]] - sum;
        }

        // Copy y to B
        for (int i = 0; i < n; i++) {
            B[i] = y[i];
        }
    }
    else if (UPLO == 'U') {
        // Backward substitution for upper triangular matrix
        x[n - 1] = B[n - 1] / A[(n - 1) * n + (n - 1)];  // B(n) / A(n,n)
        for (int i = n - 2; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += x[j] * A[i * n + j];  // A(i, i+1:n)
            }
            x[i] = (B[i] - sum) / A[i * n + i];
        }

        // Copy x to B
        for (int i = 0; i < n; i++) {
            B[i] = x[i];
        }
    }
}

void my_f(double *A,double *B,int n)
{
    int *ipiv=(int*)malloc(n*sizeof(int));
    for (int i=0;i<n;i++)
        ipiv[i]=i;
    if (mydgetrf(A,ipiv,n)==0) 
    {
        printf("LU factoration failed: coefficient matrix is singular.\n");
        return;
    }
    mydtrsv('L',A,B,n,ipiv);
    mydtrsv('U',A,B,n,ipiv);
}

#endif