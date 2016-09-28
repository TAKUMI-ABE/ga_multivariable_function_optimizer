#include <iostream>
#include <ctime>
#include "ga_util.h"


int main() {

	// init rand()
	srand((unsigned int)time(NULL));

	// configure the parameter
	GaConfiguration conf;
	conf.gtypeLength = 20;
	conf.gtypeMax = 1;
	conf.numOfVariable = 2;
	conf.population = 50;
	conf.ptypeMax = { 5.0, 4.0 };
	conf.ptypeMin = { -5.0, -4.0 };

	std::cout << "Start GA Multivariable Function Optimiser" << std::endl;

	// GA Initialization
	GaGenaration generation;
	generation.initGeneration(conf);
	generation.evaluation();

	// GA Alternation Loop of generations 
	for (int i = 0; i < NUM_OF_GENERATION; i++) {
		generation.selection();
		std::cout << i << generation << std::endl;
		generation.crossover();
		generation.mutation();
		generation.evaluation();
	}

	return 0;
}