#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <string>
#include<math.h>


int deter(int test) {
	int d = 0;
	while (test - std::pow(2,d+1)>= 0) {
		d++;
	}
	return d;
}


int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	int rang, nombre_de_processus;

	MPI_Comm_rank(MPI_COMM_WORLD, &rang);

	MPI_Comm_size(MPI_COMM_WORLD, &nombre_de_processus);
	int jeton;
	int i;
	int tag = 100;
	MPI_Status status;
	int d = std::stoi(argv[1]);


	if (rang == 0) {
		jeton = 123;
		for (i = 0; i < d; i++) {

			MPI_Send(&jeton, 1, MPI_INT, std::pow(2, i), tag, MPI_COMM_WORLD);
		}
		std::cout << "Je suis le processus "<<rang<<" et le jeton est " << jeton << std::endl;
	}

	else {
		int j = deter(rang);

		MPI_Recv(&jeton, 1, MPI_INT, rang - std::pow(2 ,j), tag, MPI_COMM_WORLD, &status);
		for (i = j + 1; i < d; i++) {

			MPI_Send(&jeton, 1, MPI_INT, std::pow(2, i) + rang, tag, MPI_COMM_WORLD);

		}
		std::cout << "Je suis le processus " << rang << " et le jeton est " << jeton << std::endl;
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
