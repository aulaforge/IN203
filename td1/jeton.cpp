#include <mpi.h>
#include <iostream>
#include <stdlib.h>

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	int rang, nombre_de_processus;

	MPI_Comm_rank(MPI_COMM_WORLD, &rang);

	MPI_Comm_size(MPI_COMM_WORLD, &nombre_de_processus);
	int jeton;
	int tag = 100;
	MPI_Status status;

	if (rang == 0) {
		jeton = 123;
		std::cout << "Je suis le rang " << rang << " et je renvoie le jeton " << jeton << "." << std::endl;
		MPI_Send(&jeton, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
		MPI_Recv(&jeton, 1, MPI_INT, nombre_de_processus - 1, tag, MPI_COMM_WORLD, &status);
		std::cout << "[ " << rang << " ]" << "J'ai reçu le jeton " << jeton <<"."<< std::endl;
	}

	else {
		std::cout << "[" << rang << "] " << "je suis en attente" << std::endl;
		MPI_Recv(&jeton, 1, MPI_INT, rang - 1, tag, MPI_COMM_WORLD, &status);
		std::cout << "[" << rang << "] "<< "j'ai recu le jeton " << jeton<< " et je renvoie " << jeton + 1 << std::endl;

		jeton += 1;
		MPI_Send(&jeton, 1, MPI_INT, (rang + 1) % nombre_de_processus, tag, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
	

