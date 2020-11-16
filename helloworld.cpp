#include <mpi.h>
#include <iostream>
#include <stdlib.h>

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	int numero_de_processus, nombre_de_processus;

	MPI_Comm_rank(MPI_COMM_WORLD, &numero_de_processus);

	MPI_Comm_size(MPI_COMM_WORLD, &nombre_de_processus);

	MPI_Finalize();

	std::cout << "Bonjour, je suis la tache numero " << numero_de_processus << " sur " << nombre_de_processus << " taches." << std::endl;

	int jeton;
	int tag = 100;
	MPI_Status status;

	if ()

	return 0;
}