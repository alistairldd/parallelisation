#include <chrono>
#include <iostream>
#include "mpi.h"


/*
mpic++ -O2 -std=c++11 tp2.cpp -o tp2
mpirun -np 2 ./tp2

*/


/* part 1
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int procid; 
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    int numproc;
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    
    printf("Hello World! rang proc : %d/%d \n", procid, numproc);
    
    MPI_Finalize();
    return 0;


}
*/

/* part 2
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int procid; 
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    int numproc;
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    
   

    if (procid%2 == 0) 
    {
        MPI_Send(&procid, 1, MPI_INT, procid+1, 0, MPI_COMM_WORLD);
        int recuPair;
        MPI_Recv(&recuPair, 1, MPI_INT, procid+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("proc : %d, res : %d \n", procid, recuPair);
    } 

    else 
    {
        int recuImpair;
        MPI_Recv(&recuImpair, 1, MPI_INT, procid-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        recuImpair += 10*procid;
        MPI_Send(&recuImpair, 1, MPI_INT, procid-1, 0, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    return 0;
}
*/

inline double f(double x)
{
  return (4.0 / (1 + x * x));
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int i;
  const int N = 100000000;
  double pi = 0.0;
  double s = 1.0 / (double)N;

  int procid, numproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);

  // Calculer le pi tel que chaque processus parcours N/size iterations de la boucle. A savoir:
  // P0: (0...N/size - 1)
  // P1: (N/size ... 2N/size - 1), etc.
  // A FAIRE ...


    
    int nb_iter = N/numproc;
    int deb = procid*nb_iter;
    int fin = deb + nb_iter;
    if (procid == numproc-1) {fin = N;}

    for (int i = deb; i < fin; i++) {
        pi = pi + s * (f(i * s) + f((i + 1) * s)) / 2;

    }  

    
    double temp = 0.;
   

    if (procid == 0){
        
        for (int i = 1; i < numproc; i++){
            
            MPI_Recv(&temp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            pi = pi + temp;

        }

        for (int i = 1; i < numproc; i++){
            MPI_Send(&pi, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            
            
        }
        // Afficher la valeur pi finale depuis le processus 0
        printf("Pi: %lf\n", pi);

    }else 
    {
        MPI_Send(&pi, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(&pi, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


  MPI_Finalize();
  return 0;
}