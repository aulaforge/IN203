# include <chrono>
# include <random>
#include <string> 
#include <iostream>
#include <mpi.h>


// Attention , ne marche qu'en C++ 11 ou supérieur :
double approximate_pi(unsigned long nbSamples) {
	typedef std::chrono::high_resolution_clock myclock;
	myclock::time_point beginning = myclock::now();
	myclock::duration d = myclock::now() - beginning;
	unsigned seed = d.count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution <double > distribution(-1.0, 1.0);
	unsigned long nbDarts = 0;
	// Throw nbSamples darts in the unit square [ -1:1] x [ -1:1]
	for (unsigned sample = 0; sample < nbSamples; ++sample) {
		double x = distribution(generator);
		double y = distribution(generator);
		// Test if the dart is in the unit disk
		if (x * x + y * y <= 1) nbDarts++;
	}
	// Number of nbDarts throwed in the unit disk
	double ratio = double(nbDarts) / double(nbSamples);
	return ratio;
}


int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	int rang, nombre_de_processus;

	MPI_Comm_rank(MPI_COMM_WORLD, &rang);

	MPI_Comm_size(MPI_COMM_WORLD, &nombre_de_processus);

	MPI_Status status;


	int tag = 100;

	double nbPoint = std::stoul(argv[1]);


	if (rang == 0) {

		double recep;

		double ratio;

		ratio = approximate_pi(nbPoint - (nombre_de_processus - 1) * floor(nbPoint / nombre_de_processus));

		int i;
		for (i = 1; i < nombre_de_processus; i++) {

			MPI_Recv(&recep, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);

			ratio += recep;
		}

		ratio = ratio / nombre_de_processus;
		std::cout << "Pi vaut environ " << 4*ratio << " ." << std::endl;

	}

	else {

		double send = approximate_pi(floor(nbPoint / nombre_de_processus));

		MPI_Send(&send, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}