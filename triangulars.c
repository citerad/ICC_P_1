/*COGNOM: CITERA NOM: DOMENICO DNI: AK1183068*/
/* 
 * File:   triangulars.c
 * Author: domenicocitera
 *
 * Created on 22 de octubre de 2013, 15:27
 * 
 * 
 * El fitxer d'entrada 'triagulars.dad' conté com el primer nombre,
 * la mida de la matriu, seguit pels elements de la matriu,
 * i, finalment, els termes coneguts.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void exit(int);
void printMatrix(int, double **, double *);
int resoltrisup(int, double **, double *, double *, double);
double residui(double **, double *, double *, int);
double *instVec(int);
double **instMat(int, int);

int main(int argc, char** argv) {

    int n = 0, i = 0, j = 0;
    double tol = 0, base = 0, exp = 0;
    double **a = NULL, *x = NULL, *b = NULL;

    //Obro streaming en llegir el fitxer
    FILE * entrada;
    entrada = fopen("triangulars.dad", "r");
    if (entrada == NULL) {
        printf("Error en obrir el fitxer %s \n", "triangulars.dad");
        exit(1);
    }

    fscanf(entrada, "%d", &n);

    //crec matriu
    a = instMat(n, n);
    x = instVec(n);
    b = instVec(n);

    //carregar matriu
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            fscanf(entrada, "%lf", &a[i][j]);
            if ((j < i && a[i][j] != 0) || (i == j && a[i][j] == 0)) {
                printf("No es un sistema triangular superior");
                exit(3);
            }
        }
    }

    //carregar b
    for (i = 0; i < n; i++) {
        fscanf(entrada, "%lf", &b[i]);
    }

    fclose(entrada); //chiudo il file

    printMatrix(n, a, b);

    printf("\n# Doneu base i exponente de la tolerancia = ");
    scanf("%lf %lf", &base, &exp);
    tol = pow(base, exp);
    printf("\n# Tollerancia = %16.7e\n", tol);

    //caluc
    if (resoltrisup(n, a, b, x, tol) != 1) {
        //stampo le x
        for (i = 0; i < n; i++)
            printf("\nx%d = %16.7e", i + 1, x[i]);

        // residuo
        double err = residui(a, x, b, n);
        printf("\n\n# Residuo: = %16.7e\n", err);

    } else {
        printf("# El calcul no es pot fer\n");
        exit(4);
    }

    free(a);
    free(b);
    free(x);
    return 0;
}

/**
 * Risolve un sistema di equazioni sapendo che A è triangolare superiore
 * @param n dimensione matrice
 * @param A matrice
 * @param b termini noti
 * @param x puntatore per risultati
 * @param tol tollerancia
 * @return 0 se ha successo 1 altrimenti
 */
int resoltrisup(int n, double **A, double *b, double *x, double tol) {
    //resoldre
    int j = 0, k = 0;
    double sum = 0;

    for (j = n - 1; j >= 0; j--) {
        if (fabs(A[j][j]) < tol)
            return 1;

        for (k = j + 1; k < n; k++) {
            sum = sum + (A[j][k] * x[k]);
        }
        x[j] = (b[j] - sum) / A[j][j];
        sum = 0;
    }
    return 0;
}

/**
 * imprimeix la matriu completa
 * @param n dimensione
 * @param a matrice coefficienti
 * @param vettori termini noti
 */
void printMatrix(int n, double **a, double *b) {
    int i, j;

    printf("# La matrice completa è\n\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%16.7e ", a[i][j]);
        }
        printf("| %16.7e\n", b[i]);
    }
}

/**
 * instanazia memoria ad una matrice
 * 
 * @param n righe
 * @param m colonne
 * @return puntatore alla matrice
 */
double **instMat(int n, int m) {
    int i = 0, x = 0;
    double **M = NULL;

    M = (double **) malloc(n * sizeof (double *));
    if (M != NULL) {
        for (i = 0; i < n; i++) {
            M[i] = (double *) malloc(m * sizeof (double));
            if (M[i] == NULL) { //Se non ho memoria libero
                for (x = 0; x < i; x++) {
                    free(M[x]);
                }
                free(M);
                M = NULL;
                printf("No hi ha prou memoria");
                exit(2);
            }
        }
    }
    return M;
}

/**
 * Istanzia memoria ad un vettore
 * @param n dimensione
 * @return puntatore all'array
 */
double *instVec(int n) {
    double *v = (double *) malloc(n * sizeof (double));
    if (v == NULL) {
        printf("No hi ha prou memoria");
        exit(2);
    }
    return v;
}

/**
 * calcul el residu
 * @param a matrice coefficienti
 * @param x valori delle incognite
 * @param b termini noti
 * @param n dimensione matrice
 * @return residuo
 */
double residui(double **a, double *x, double *b, int n) {

    int i = 0, j = 0;
    double ris = 0, temp = 0;

    //A*x in e
    for (i = 0; i < n; i++) {
        temp = 0;
        for (j = 0; j < n; j++) {
            temp = temp + a[i][j] * x[j];
        }
        temp = fabs(temp - b[i]);
        ris += temp*temp;
    }
    return sqrt(ris);
}

/**
 * calcul el residu i el posa en e
 */
/*void calcoloResiduoError(double **a, double *x, double *b, double *e, int n) {

    int i, j;
    //A*x in e
    for (i = 0; i < n; i++) {
        e[i] = 0.;
        for (j = 0; j < n; j++) {
            e[i] += a[i][j] * x[j];
        }
    }
    for (i = 0; i < n; i++){
        e[i] = fabs(e[i] - b[i]);
        printf("\n   x%d = %16.7e;  e%d = %16.7e", i+1, x[i], i+1, e[i]);
    }
    
    double ea, er, sume = 0, sumb = 0;
    for (i = 0; i < n; i++) {
        sume = sume + pow(e[i], 2);
        sumb = sumb + pow(b[i], 2);
    }
    ea = sqrt(sume);
    er = ea/sqrt(sumb);
    printf("\n\n   Ea = %16.7e \n   Er = %16.7e \n", ea, er);
}*/

