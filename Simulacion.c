#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "Estructuras.h"

const double alfa = 0.01;

//Generador de números de Parisi-Rapuano.
//Reproducción del código de los apuntes de la asignatura "Física Computacional"
//impartida por Alfonso Tarancon Lafita (Universidad de Zaragoza)
 void ini_ran(int SEMILLA)
 {
    int INI,FACTOR,SUM,i;
    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;
    for(i=0;i<256;i++)
    {
    INI=(INI*FACTOR+SUM);
    irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
 }

 float Random(void)
 {
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;
    return r;
 }

 //Generador de los números aleatorios según distribución gaussiana por el método de Box-Muller
double DistrGauss(double Av, double Desv){
    double d1, d2;
    //Se generan los numeros planos necesarios
    d1 = Random();
    d2 = Random();
    if(d1 == 0)
        d1 = d2;
    //Se usa la fórmula para crear dos números independientes según una gaussiana
    return Av -Desv*sqrt(-2*log(d1))*cos(2*Pi*d2);
}

double Reward(Nodo *vecino, int opinion){

    //Comprobamos si el vecino habla
    if(vecino->prob > Random()){

        //Si ambos nodos tienen la misma opinion, reward 1.
        if(vecino->opinion == opinion)
            return +1;

        //Si los nodos tienen opiniones diferentes, reward -1.
        else
            return -1;
    }

    //El vecino no habla.
    else
        return 0;
}

void PasoMonteCarlo(Lattice *red, double coste, double beta){
    Nodo *nodo = NULL, *vecino = NULL;
    double r =0;
    int n = 0;
    for(int i=0; i<red->NNodos; i++){
            r=0;
        n=red->NNodos*Random();
        nodo = &red->nodos[n];
        //Comprobamos si el nodo escogido habla en este paso temporal.
       // if(nodo->prob>Random()){
        for(int k=0; k<20; k++){
            r += -coste;
            n = nodo->NVecinos*Random();
            vecino = &red->nodos[nodo->vecinos[n]];
            r += Reward(vecino, nodo->opinion);
        }
            nodo->Q = (1-alfa)*nodo->Q + alfa*r/20;
            nodo->prob = 1.0/(1+exp(-beta*nodo->Q));
        //}
    }
    nodo = NULL;
    vecino = NULL;
}

double AverageQ(Lattice *red, int opinion){
    double average = 0.0;
    int N;
    if(opinion == 1){
        N = red->NNodos - red->NNodos0;
    }
    else
        N = red->NNodos0;

    for(int i=0; i<red->NNodos; i++){
        if(red->nodos[i].opinion==opinion){
            average += red->nodos[i].Q;
        }
    }
    return average/N;
}

void Histograma(char *filename, Lattice *lattice){
    FILE *histograma = fopen(filename, "w");
    double QMax = +1, QMin = -1;
    int NBins = 100, hist0[NBins], index, hist1[NBins];
    double delta = (QMax-QMin)/NBins ;
    for(int i=0; i<NBins; i++){
        hist0[i] =0;
        hist1[i] = 0;
    }
    for(int i=0; i<lattice->NNodos; i++){
        if(lattice->nodos[i].opinion == 0){
            index = (lattice->nodos[i].Q - QMin)/delta;
            if(index >=0 && index <NBins)
                hist0[index]++;
        }
        else{
            index = (lattice->nodos[i].Q - QMin)/delta;
            if(index >=0 && index <NBins)
                hist1[index]++;
        }
    }
    for(int i=0; i<NBins; i++)
        fprintf(histograma, "%lf\t%lf\t%lf\n", QMin + (i+0.5)*delta, hist0[i]/delta/lattice->NNodos0, hist1[i]/delta/(lattice->NNodos - lattice->NNodos0));
    fclose(histograma);

}

void Simulacion(int TTotal, Lattice *lattice, double coste, double beta, double *Q){
    double QAv[2];
    QAv[0] = Q[0];
    QAv[1] = Q[1];
    for(int i=0; i<lattice->NNodos; i++){
        lattice->nodos[i].Q = DistrGauss(QAv[lattice->nodos[i].opinion],0.01);
        lattice->nodos[i].prob = 1.0/(1+exp(-beta*lattice->nodos[i].Q));
    }
    for(int t=0; t<=TTotal; t++){
        PasoMonteCarlo(lattice, coste, beta);
    }
    Q[0] = AverageQ(lattice, 0);
    Q[1] = AverageQ(lattice, 1);
}

Lattice* CreacionRedErdos(int nnodos, double prob) {
    Lattice *red = malloc(sizeof(Lattice));
    red->nodos = malloc(nnodos * sizeof(Nodo));
    red->NNodos = nnodos;

    // Inicializar nodos
    for (int i = 0; i < nnodos; i++) {
        red->nodos[i].vecinos = NULL;
        red->nodos[i].NVecinos = 0;
    }

    // Construir conexiones
    for (int i = 0; i < nnodos; i++) {
        for (int j = i + 1; j < nnodos; j++) {
            if (Random() < prob) {
                // Añadir j como vecino de i
                Nodo *ni = &red->nodos[i];
                ni->vecinos = realloc(ni->vecinos, (ni->NVecinos + 1) * sizeof(int));
                ni->vecinos[ni->NVecinos++] = j;

                // Añadir i como vecino de j
                Nodo *nj = &red->nodos[j];
                nj->vecinos = realloc(nj->vecinos, (nj->NVecinos + 1) * sizeof(int));
                nj->vecinos[nj->NVecinos++] = i;
            }
        }
    }

    return red;
}

void DFS_iterativo(Lattice* red, int inicio, int* visitado, int* componente, int* tam) {
    int N = red->NNodos;
    int* stack = malloc(N * sizeof(int));
    int top = 0;
    stack[top++] = inicio;
    visitado[inicio] = 1;

    while (top > 0) {
        int nodo = stack[--top];
        componente[(*tam)++] = nodo;

        Nodo actual = red->nodos[nodo];
        for (int i = 0; i < actual.NVecinos; i++) {
            int vecino = actual.vecinos[i];
            if (!visitado[vecino]) {
                visitado[vecino] = 1;
                stack[top++] = vecino;
            }
        }
    }

    free(stack);
}


Lattice* ExtraerComponenteMayor(Lattice* red) {
    int N = red->NNodos;
    int* visitado = calloc(N, sizeof(int));
    int* mejor_componente = malloc(N * sizeof(int));
    int mejor_tam = 0;

    int* temp_componente = malloc(N * sizeof(int));
    int tam = 0;

    for (int i = 0; i < N; i++) {
        if (!visitado[i]) {
            tam = 0;
            DFS_iterativo(red, i, visitado, temp_componente, &tam);
            if (tam > mejor_tam) {
                mejor_tam = tam;
                memcpy(mejor_componente, temp_componente, tam * sizeof(int));
            }
        }
    }

    free(temp_componente);
    free(visitado);

    // Crear nueva red con la mayor componente
    Lattice* nueva = malloc(sizeof(Lattice));
    nueva->NNodos = mejor_tam;
    nueva->nodos = malloc(mejor_tam * sizeof(Nodo));
    nueva->NNodos0 = 0;
    // Mapeo de nodos originales a nuevos índices
    int* mapa = malloc(N * sizeof(int));
    for (int i = 0; i < N; i++)
        mapa[i] = -1;
    for (int i = 0; i < mejor_tam; i++)
        mapa[mejor_componente[i]] = i;

    for (int i = 0; i < mejor_tam; i++) {
        int old = mejor_componente[i]; //Indice del nodo en la red general
        Nodo* nodo = &red->nodos[old]; //Dirección del nodo original
        nueva->nodos[i].NVecinos = 0; //Inicializacion de la variable en la nueva red
        nueva->nodos[i].vecinos = malloc(nodo->NVecinos * sizeof(int)); // Máximo posible
        nueva->nodos[i].opinion = nodo->opinion;
        if(nodo->opinion == 0){
            nueva->NNodos0++;
        }
        for (int j = 0; j < nodo->NVecinos; j++) {
            int v = nodo->vecinos[j];
            if (mapa[v] != -1) {
                nueva->nodos[i].vecinos[nueva->nodos[i].NVecinos++] = mapa[v];
            }
        }

        // Redimensionar memoria
        nueva->nodos[i].vecinos = realloc(nueva->nodos[i].vecinos,
                                          nueva->nodos[i].NVecinos * sizeof(int));
    }

    free(mejor_componente);
    free(mapa);
    return nueva;
}


Lattice* Upgrading(Lattice *red, double s0, double s1, double s01, double s10, double c){
    double tabla[4] = {s0, s01, s10, s1};  // valores de s para las combinaciones
    double s;
    Nodo *nodo = NULL, *candidato = NULL;
    int n, n2, n3;
    int conectado;

    for(int i=0; i<red->NNodos; i++){
        n=red->NNodos*Random();
        nodo = &red->nodos[n];
        while(nodo->NVecinos == 0){
            n=red->NNodos*Random();
            nodo = &red->nodos[n];
        }
        if(Random()< c){
            n2 = n;
            candidato = nodo;

            n2 = nodo->NVecinos*Random();
            candidato = &red->nodos[red->nodos[n].vecinos[n2]];
            n2 = candidato->NVecinos*Random();
            n2 = candidato->vecinos[n2];
            if(n2==n)
                continue;
            candidato = &red->nodos[n2];
        }
        else{
            n2 = red->NNodos*Random();
            if(n2==n)
                continue;
            candidato = &red->nodos[n2];
        }
        s = tabla[2*nodo->opinion + candidato->opinion];
        if(s>Random()){
            // Verificar si ya existe una conexión
            conectado = 0;
            for (int k = 0; k < candidato->NVecinos; k++) {
                if (candidato->vecinos[k] == n) {
                    conectado = 1;
                    break;
                }
            }

            if (!conectado) {
                //Actualizo el nodo conectado con la nueva conexión.
                candidato->vecinos = (int *)realloc(candidato->vecinos, sizeof(int) * (candidato->NVecinos + 1));
                if (candidato->vecinos == NULL) {
                    fprintf(stderr, "Error al hacer realloc\n");
                    exit(EXIT_FAILURE);  // O maneja el error de forma apropiada
                }
                candidato->vecinos[candidato->NVecinos] = n;
                candidato->NVecinos++;

                //Actualizo el nodo que se conecta.

                    n3 = nodo->NVecinos*Random();
                    int vecino_original = nodo->vecinos[n3];
                    Nodo *vecino_a_desconectar = &red->nodos[vecino_original];
                    nodo->vecinos[n3] = n2;

                    //Actualizo el nodo cuya conexión se rompe.
                    int o, found = 0;
                    // Buscar el índice donde está el valor n
                    for (o = 0; o < vecino_a_desconectar->NVecinos; o++) {
                        if (vecino_a_desconectar->vecinos[o] == n) {
                            found = 1;
                            break;
                        }
                    }

                    if (found) {
                        // Desplazar los elementos a la izquierda para eliminar el elemento n
                        for (int j = o; j < vecino_a_desconectar->NVecinos - 1; j++) {
                            vecino_a_desconectar->vecinos[j] = vecino_a_desconectar->vecinos[j + 1];
                        }
                        vecino_a_desconectar->NVecinos--;
                        // Redimensionar el array
                        vecino_a_desconectar->vecinos = (int *)realloc(vecino_a_desconectar->vecinos, sizeof(int) * vecino_a_desconectar->NVecinos);
                        if (vecino_a_desconectar->NVecinos > 0 && vecino_a_desconectar->vecinos == NULL) {
                            fprintf(stderr, "Error al hacer realloc tras eliminar vecino\n");
                            exit(EXIT_FAILURE);
                        }
                    }

            }
        }
    }
    nodo = NULL;
    candidato = NULL;
    return red;
}

double Homophily(Lattice *red, int Op){
    double media = 0;
    int count =0;
    for(int i=0; i<red->NNodos; i++){
        if(red->nodos[i].opinion == Op && red->nodos[i].NVecinos>0){
            count++;
            for(int j=0; j<red->nodos[i].NVecinos; j++){
                if(red->nodos[red->nodos[i].vecinos[j]].opinion==Op)
                    media += 1.0/red->nodos[i].NVecinos;
            }
        }

    }
    if(count != 0)
        return media/count;
    return 0;
}

Lattice* TriadicClosure(int TTotal, int nnodos, double n, double k, double c, double s0, double s1, double s01, double s10){
    Lattice *red = CreacionRedErdos(nnodos, k/(nnodos-1));
    red->NNodos0 = 0;
    for(int i=0; i<nnodos; i++){
        if(Random()<n){
            red->nodos[i].opinion = 0;
            red->NNodos0++;
        }
        else
            red->nodos[i].opinion = 1;
    }
    for(int t=0; t<TTotal; t++)
            Upgrading(red, s0, s1, s01, s10, c);
    return red;
}

Lattice *LeerRedFichero(void){
 FILE *f=fopen("estructura_red.txt", "r");
 if(f!=NULL){
    Lattice *red = (Lattice *)malloc(sizeof(Lattice));
    red->nodos = (Nodo *)malloc(1000*sizeof(Nodo));
    red->NNodos = 1000;
    red->NNodos0=0;
    for(int i=0; i<1000; i++){
        fscanf(f, "%d\t%d", &red->nodos[i].opinion, &red->nodos[i].NVecinos);
        red->nodos[i].vecinos = (int *)malloc(red->nodos[i].NVecinos * sizeof(int));
        for(int j=0; j< red->nodos[i].NVecinos; j++){
            fscanf(f, "\t %d", &red->nodos[i].vecinos[j]);
        }
        fscanf(f, "\n");
    }
    fclose(f);
    return red;
 }
 return NULL;
}

int EsVecino(Nodo *nodos, int a, int b) {
    for (int i = 0; i < nodos[a].NVecinos; i++) {
        if (nodos[a].vecinos[i] == b) return 1;
    }
    return 0;
}

void ClusteringPorOpinion(Lattice *red, double *clust_op0, double *clust_op1) {
    double suma_op0 = 0.0, suma_op1 = 0.0;
    int cuenta_op0 = 0, cuenta_op1 = 0;

    for (int i = 0; i < red->NNodos; i++) {
        Nodo nodo = red->nodos[i];
        int k = nodo.NVecinos;

        if (k < 2) continue; // No se puede formar triángulo

        int conexiones = 0;
        // Comparar todos los pares de vecinos
        for (int v1 = 0; v1 < k; v1++) {
            for (int v2 = v1 + 1; v2 < k; v2++) {
                int n1 = nodo.vecinos[v1];
                int n2 = nodo.vecinos[v2];
                if (EsVecino(red->nodos, n1, n2)) {
                    conexiones++;
                }
            }
        }

        double max_conexiones = (k * (k - 1)) / 2.0;
        double C = conexiones / max_conexiones;

        if (nodo.opinion == 0) {
            suma_op0 += C;
            cuenta_op0++;
        } else {
            suma_op1 += C;
            cuenta_op1++;
        }
    }

    *clust_op0 = (cuenta_op0 > 0) ? (suma_op0 / cuenta_op0) : 0.0;
    *clust_op1 = (cuenta_op1 > 0) ? (suma_op1 / cuenta_op1) : 0.0;
}

void GuardarRedFichero(Lattice *red) {
    FILE *f = fopen("estructura_red.txt", "w");
    if (f == NULL) {
        perror("Error al abrir el archivo para escribir");
        return;
    }

    for (int i = 0; i < red->NNodos; i++) {
        fprintf(f, "%d\t%d", red->nodos[i].opinion, red->nodos[i].NVecinos);
        for (int j = 0; j < red->nodos[i].NVecinos; j++) {
            fprintf(f, "\t%d", red->nodos[i].vecinos[j]);
        }
        fprintf(f, "\n");
    }

    fclose(f);
}

void FreeLatice(Lattice *red){
    for(int i=0; i<red->NNodos; i++){
        free(red->nodos[i].vecinos);
        red->nodos[i].vecinos=NULL;
    }
    free(red->nodos);
    free(red);
}


