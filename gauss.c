/*COGNOM: CITERA NOM: DOMENICO DNI: AK1183068*/
/* 
 * File:   gauss.c
 * Author: domenicocitera
 *
 * Created on 22 de octubre de 2013, 15:27
 * 
 * L'arxiu 'gaus.dad' d'entrada conté:
 * primer número = mida de la matriu de coeficients,
 * segon número = mida de la matriu de termes coneguts
 * seguir els nxn elements de la matriu de coeficients 
 * i els nxm elements de matriu dels termes coneguts (cada equació a la columna)
 * següent el nombre que indica la mida de la matriu per a calcular la inversa,
 * i els nxn elements de matriu.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void exit(int);
void printMatrix(int, int, double **);
int gauss(double **, double **, int, int, double);
double **instMat(int, int);
double *instVec(int);
double **copiaMatrice(int, int, double **);
int resoltrisup(int, double **, double *, double *, double);
double residui(double **, double *, double *, int);
int calcolaInversa(double **, int, double);
int triangolarizza(double **, double*, int, int, int);

int main(int argc, char** argv) {

    int n = 0, m = 0, i = 0, j = 0;
    double tol = 0, base = 0, esponente = 0;
    double **a, **b;

    FILE * entrada; //apro stream
    entrada = fopen("gauss.dad", "r");
    if (entrada == NULL) {
        printf("Error en obrir el fitxer %s \n", "gauss.dad");
        exit(1);
    }

    fscanf(entrada, "%d %d", &n, &m); //indica dim matrice

    a = instMat(n, n);
    b = instMat(n, m);

    for (i = 0; i < n; i++) { //inserisco nella matrice
        for (j = 0; j < n; j++) {
            fscanf(entrada, "%lf", &a[i][j]);
        }
    }

    for (i = 0; i < n; i++) { //inserisco termini noti
        for (j = 0; j < m; j++) {
            fscanf(entrada, "%lf", &b[i][j]);
        }

    }

    printf("# A:\n");
    printMatrix(n, n, a);
    printf("\n# B:\n");
    printMatrix(n, m, b);

    printf("\n# Doneu base i exponente de la tolerancia = ");
    scanf("%lf %lf", &base, &esponente);
    tol = pow(base, esponente);
    printf("\n# Tollerancia %16.7e\n", tol);


    if (gauss(a, b, n, m, tol) == 1) {
        printf("# El calcul no es pot fer\n");
        exit(4);
    }
    free(b);

    fscanf(entrada, "%d", &n);
    a = instMat(n, n);
    for (i = 0; i < n; i++) { //inserisco nella matrice
        for (j = 0; j < n; j++) {
            fscanf(entrada, "%lf", &a[i][j]);
        }
    }
    printf("\n# C: \n");
    printMatrix(n, n, a);
    if (calcolaInversa(a, n, tol) == 1) {
        printf("# El calcul no es pot fer\n");
        exit(4);
    }
    fclose(entrada); //chiudo il file

    return 0;
}

/**
 * Stampa la matrice
 * @param n righe
 * @param m colonne
 * @param a matrice
 */
void printMatrix(int n, int m, double **a) {
    int i, j;
    printf("\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            printf("%16.7e ", a[i][j]);
        }
        printf("\n");
    }
}

/**
 * Copia una matrice
 * @param n righe
 * @param m colonne
 * @param M matrice
 * @return puntatore a nuova matrice copiata
 */
double **copiaMatrice(int n, int m, double **M) {
    int i = 0, j = 0;
    double **cm = NULL;

    cm = instMat(n, m);
    if (cm != NULL && M != NULL) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                cm[i][j] = M[i][j];
            }
        }
    }
    return cm;
}

/**
 * Risolve m sistemi lineari con la stessa matrice A e termini indipendenti in B
 * @param A matrice coefficienti
 * @param B matrice termini noti
 * @param n numero incognite
 * @param m numero vettori termini noti
 * @param tol tolleranza
 * @return 0 se ha successo
 */
int gauss(double **A, double **B, int n, int m, double tol) {
    int i = 0, k = 0, l = 0, j = 0;
    double **mat = instMat(n, n);
    double *e = instVec(m);

    for (i = 0; i < m; i++) { //per ogni colonna di B

        double **cda = copiaMatrice(n, n, A);
        double *b = instVec(n);
        double *x = instVec(n);


        for (j = 0; j < n; j++) //copio colonna B in b
            b[j] = B[j][i];

        for (k = 0; k < n - 1; k++) {
            if (fabs(cda[k][k]) < tol) {
                printf("\n# El calcul no se pot fer\n");
                return 1;
            }

            for (l = k + 1; l < n; l++) {
                mat[l][k] = cda[l][k] / cda[k][k];
                cda[l][k] = 0;

                for (j = k + 1; j < n; j++) {
                    cda[l][j] = cda[l][j] - (mat[l][k] * cda[k][j]);
                }
                b[l] = b[l] - (mat[l][k] * b[k]);
            }
        }

        if (resoltrisup(n, cda, b, x, tol) == 1) { //risolvo
            printf("\n# El calcul no se pot fer\n");
            return 1;
        }

        for (j = 0; j < n; j++) { //sostituisco in B
            B[j][i] = x[j];
        }
        e[i] = residui(cda, x, b, n); //calcolo residuo
    }

    for (i = 0; i < m; i++) {
        printf("\n#   Soluzioni sistema %d:\n", i + 1);

        for (j = 0; j < n; j++) //stampo x
            printf("\n    x%d = %16.7e", j + 1, B[j][i]);
        //stampo residui
        printf("\n\n#    Residuo sistema %d: = %16.7e\n\n", i + 1, e[i]);
    }
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
    int j, k;
    double sum = 0;
    for (j = n - 1; j >= 0; j--) {
        if (fabs(A[j][j]) < tol)
            return 1;

        for (k = j + 1; k < n; k++) {
            sum = sum + A[j][k] * x[k];
        }
        x[j] = (1 / A[j][j])*(b[j] - sum);
        sum = 0;
    }
    return 0;
}

/**
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
 * calcola l'inversa della matrice A se possibile
 * 
 * @param A matrice
 * @param n righe e colonne
 * @param tol tolleranza
 */
int calcolaInversa(double **A, int n, double tol) {

    int i = 0, j = 0;

    //controllo determinante
    double **cda = copiaMatrice(n, n, A);
    if (triangolarizza(cda, NULL, n, n, tol) == 1) {
        printf("# El calcul no es pot fer\n");
        return 1;
    }

    for (i = 0; i < n; i++) {
        if (cda[i][i] == 0) {
            printf("# Det=0, non ha inversa");
            return 2;
        }
    }
    free(cda);
    double **inversa = instMat(n, n);
    double **completa = instMat(n, n * 2);

    //affianco matrice identità
    for (i = 0; i < n; i++) {
        for (j = 0; j < n * 2; j++) {
            if (j < n)
                completa[i][j] = A[i][j];
            else if (i + n != j)
                completa[i][j] = 0;
            else completa[i][j] = 1;
        }
    }

    //triangolarizzo
    if (triangolarizza(completa, NULL, n, n * 2, tol) == 1) {
        printf("# El calcul no es pot fer\n");
        return 1;
    }

    //salvo coefficienti
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i][j] = completa[i][j];
        }
    }
    for (i = 0; i < n; i++) {
        for (j = n; j < n * 2; j++) {
            inversa[i][j - n] = completa[i][j];
        }
    }
    free(completa);
    if (gauss(A, inversa, n, n, tol) == 1) {
        printf("# El calcul no es pot fer\n");
        return 1;
    }
    printf("# L'inversa è: \n");
    printMatrix(n, n, inversa);
    return 0;
}

/**
 * triangolarizza mediante mosse di gauss la matrice A e i termini noti in b se
 * presenti
 * @param A matrice coefficienti
 * @param B termini noti
 * @param n dimensione righe
 * @param m dimensioni colonne
 * @param tol tolleranza
 * @return 0 se ha successo 1 altrimenti
 */
int triangolarizza(double **A, double *B, int n, int m, int tol) {

    int j = 0, k = 0, l = 0;
    double **mat = instMat(n, m);

    //applico gauss
    for (k = 0; k < n - 1; k++) {
        if (fabs(A[k][k]) < tol) {
            return 1;
        }
        for (l = k + 1; l < n; l++) {
            mat[l][k] = A[l][k] / A[k][k];
            A[l][k] = 0;
            for (j = k + 1; j < m; j++) {
                A[l][j] = A[l][j]-(mat[l][k] * A[k][j]);
            }
            if (B != NULL)
                B[l] = B[l] - (mat[l][k] * B[k]);
        }
    }
    return 0;
}
