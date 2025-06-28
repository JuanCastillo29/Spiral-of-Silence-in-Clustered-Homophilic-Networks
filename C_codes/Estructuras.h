#ifndef ESTRUCTURAS_H_INCLUDED
#define ESTRUCTURAS_H_INCLUDED

unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;
#define NormRANu (2.3283063671E-10F)
#define Pi 3.14159265358979323846
#define NSim 1
typedef struct{
    int *vecinos;
    double Q;
    int NVecinos;
    int opinion;
    double prob;
}Nodo;

typedef struct{
    int NNodos, NNodos0;
    Nodo *nodos;
    double beta;
}Lattice;

#endif // ESTRUCTURAS_H_INCLUDED
