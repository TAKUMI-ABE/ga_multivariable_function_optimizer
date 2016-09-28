#ifndef GA_UTIL_H
#define GA_UTIL_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "ga_target.h"

// params of GA
const int NUM_OF_GENERATION = 50;
const double ELITE_RATE = 0.1;
const double MUTATE_RATE = 0.02;



/*
 * Geno Type
 */
class Gtype {
public:
	Gtype() {};

	void setGtypeRamdom(int gtypeCodeLength, int gtypeCodeMax, int numOfVariable);
	void crossoverGtype(Gtype& g, int intersection);
	int mutateGtype();
	int getGCodeLength() { return gCodeLength; }
	int getGCodeMax() { return gCodeMax; }
	int getNumOfVariable() { return static_cast<int>(gVector.size()); }

	std::vector<std::vector<int>> gVector;

private:
	int gCodeLength;
	int gCodeMax;
};

void Gtype::setGtypeRamdom(int gtypeCodeLength, int gtypeCodeMax, int numOfVariable) {
	gCodeLength = gtypeCodeLength;
	gCodeMax = gtypeCodeMax;
	gVector.resize(numOfVariable);

	for (int i = 0; i < numOfVariable; i++) {
		gVector[i].resize(gtypeCodeLength);
		for (int j = 0; j < gtypeCodeLength; j++)
			gVector[i][j] = rand() % (gtypeCodeMax + 1);
	}
}

void Gtype::crossoverGtype(Gtype & g, int intersection){
	for (int i = 0; i < getNumOfVariable(); i++)
		for (int j = 0; j < getGCodeLength(); j++)
			if (j < intersection)
				gVector[i][j] = g.gVector[i][j];
}

int Gtype::mutateGtype() {
	int mutateCount = 0;
	for (int i = 0; i < getNumOfVariable(); i++) {
		for (int j = 0; j < gCodeLength; j++) {
			if (((double)rand() / RAND_MAX) <= MUTATE_RATE) {
				gVector[i][j] = (gVector[i][j] == 1) ? 0 : 1;
				mutateCount++;
			}
		}
	}

	return mutateCount;
}

std::ostream &operator<<(std::ostream &out, Gtype& g) {
	out << "[";
	for (int j = 0; j < g.getGCodeLength(); j++)
		std::cout << g.gVector[0][j];
	for (int i = 1; i < g.getNumOfVariable(); i++) {
		std::cout << "|";
		for (int j = 0; j < g.getGCodeLength(); j++)
			std::cout << g.gVector[i][j];
	}
	out << "]";

	return out;
}

/*
 * Pheno Type
 */
class Ptype {
public:
	Ptype() {};

	void decodeGtype(Gtype& g);
	void setRange(std::vector<double> min, std::vector<double> max, int numVariable);
	std::vector<double> getPtype() { return p; }
	int getNumOfVariable() { return static_cast<int>(p.size()); }

private:
	std::vector<double> p;
	std::vector<double> pMin;
	std::vector<double> pMax;
};

void Ptype::decodeGtype(Gtype& g) {
	for (int i = 0; i < g.getNumOfVariable(); i++) {
		double gap = pMax[i] - pMin[i];
		double decoded_value = pMin[i];


		for (int j=0;j<g.getGCodeLength();j++)
			if (g.gVector[i][j])
				decoded_value += gap / pow(2, j + 1);

		p[i] = decoded_value;
	}
}

void Ptype::setRange(std::vector<double> min, std::vector<double> max, int numVariable) {
	pMin = min;
	pMax = max;
	p.resize(numVariable);
}

std::ostream &operator<<(std::ostream &out, Ptype& p) {
	out << p.getPtype()[0];
	for (int i = 1; i < p.getNumOfVariable(); i++)
		out << ", " << p.getPtype()[i];

	return out;
}

/*
 * individual
 */
class GaIndividual {
public:
	GaIndividual() {};

	void setGaIndividual(int gtypeCodeLength, int gtypeCodeMax, std::vector<double> ptyeMin, std::vector<double> ptypeMax, int numVariable);
	void calcFitness(Ptype& p);

	static bool compareGaIndividualPredicate(GaIndividual a, GaIndividual b) { return (a.fitness > b.fitness); }

	double getFitness() { return fitness; }
	int getRank() { return rank; }
	int getParent1() { return parent1; }
	int getParent2() { return parent2; }
	int getCrossPoint() { return crossPoint; }
	void setRank(int gRank) { rank = gRank; }
	void setParent1(int p1) { parent1 = p1; }
	void setParent2(int p2) { parent2 = p2; }
	void setCrossPoint(int cp) { crossPoint = cp; }

	Gtype gtype; // Geno Type
	Ptype ptype; // Pheno Type

private:
	double fitness;
	int rank;  // rank after sorting
	int parent1; // index of parent 1
	int parent2; // index of parent 2
	int crossPoint; // crossover point
};

void GaIndividual::setGaIndividual(int gtypeCodeLength, int gtypeCodeMax, std::vector<double> ptyeMin, std::vector<double> ptypeMax, int numVariable){
	rank = 0;
	parent1 = 0;
	parent2 = 0;
	crossPoint = 0;

	gtype.setGtypeRamdom(gtypeCodeLength, gtypeCodeMax, numVariable);
	ptype.setRange(ptyeMin, ptypeMax, numVariable);
}

void GaIndividual::calcFitness(Ptype& p) {
	fitness = gy(fx(p.getPtype()));
}

/*
 * a generation
 */
class GaGenaration {
public:
	GaGenaration() {};

	void initGeneration(GaConfiguration& conf);
	void coutResult(int num);
	int getPopulation() { return static_cast<int>(genes.size()); }

	void evaluation();
	void selection();
	void crossover();
	void mutation();

private:

	GaIndividual decideParentRoulette(std::vector<GaIndividual>& candidate);
	GaIndividual getChildCrossOver(GaIndividual& parent1, GaIndividual& parent2);

	std::vector<GaIndividual> genes;
	int mutateCount;  // total number of mutation
	double maxFitness;
	double minFitness;
	double avgFitness;
};

void GaGenaration::initGeneration(GaConfiguration& conf) {
	mutateCount = 0;
	maxFitness = 0.0;
	minFitness = 0.0;
	avgFitness = 0.0;
	genes.resize(conf.population);

	for (int i = 0; i < getPopulation(); i++)
		genes[i].setGaIndividual(conf.gtypeLength, conf.gtypeMax, conf.ptypeMin, conf.ptypeMax, conf.numOfVariable);
}

void GaGenaration::evaluation() {
	// calc ptype from gtype
	for (int i = 0; i < getPopulation(); i++)
		genes[i].ptype.decodeGtype(genes[i].gtype);
}

void GaGenaration::selection() {
	// calc fitness of each genes
	avgFitness = 0.0;
	for (int i = 0; i < getPopulation(); i++) {
		genes[i].calcFitness(genes[i].ptype);
		avgFitness += genes[i].getFitness();
	}
	avgFitness /= genes.size();

	// sort genes by fitness
	std::sort(genes.begin(), genes.end(), GaIndividual::compareGaIndividualPredicate);

	// identify genes
	for (int i = 0; i < getPopulation(); i++)
		genes[i].setRank(i);

	// calc min and max of fitness
	maxFitness = genes[0].getFitness();
	minFitness = genes[getPopulation() - 1].getFitness();
}

void GaGenaration::coutResult(int num) {

	std::cout << "----------------------------------------------------" << std::endl;
	std::cout << "###  parent   xsite  gtype                    ptype   fitness" << std::endl;
	int count = 0;
	for (auto gene : genes) {
		std::cout << std::setw(3) << count << " (" << std::setw(3) << gene.getParent1() << "," << std::setw(3) << gene.getParent2() << ")  " << std::setw(3) << gene.getCrossPoint() << "   "
			<< gene.gtype << " " << std::setw(6) << gene.ptype << "    " << std::setw(6) << gene.getFitness() << std::endl;
		count++;
	}
	std::cout << std::endl << "Mutation: " << mutateCount << std::endl
		<< num << ", " << maxFitness << ", " << avgFitness << ", " << minFitness
		<< ", " << genes[0].ptype << ", "
		<< genes[0].gtype << std::endl;
}

void GaGenaration::crossover() {
	// Reincarnation of Elite
	for (int i = 0; i < getPopulation()*ELITE_RATE; i++) {
		genes[i].setParent1(0);
		genes[i].setParent2(0);
		genes[i].setCrossPoint(0);
	}

	// decide parent
	std::vector<GaIndividual> parent1;
	std::vector<GaIndividual> parent2;
	GaIndividual p2Candidate;
	for (int i = 0; i < getPopulation()*(1.0 - ELITE_RATE); i++) {
		parent1.push_back(decideParentRoulette(genes));
		p2Candidate = decideParentRoulette(genes);
		while (p2Candidate.getRank() == parent1[i].getRank()) {
			p2Candidate = decideParentRoulette(genes);
		}
		parent2.push_back(p2Candidate);
	}

	// cross over parent
	std::vector<GaIndividual> children;
	for (int i = 0; i < getPopulation()*(1.0 - ELITE_RATE); i++)
		children.push_back(getChildCrossOver(parent1[i], parent2[i]));

	// modify generation
	int numOfChildren = static_cast<int>(children.size());
	for (int i = 0; i < numOfChildren; i++)
		genes[getPopulation() - children.size() + i] = children[i];
}

void GaGenaration::mutation() {
	mutateCount = 0;
	for (int i = 0; i < getPopulation()*(1.0 - ELITE_RATE); i++)
		mutateCount += genes[static_cast<int>(getPopulation()*ELITE_RATE + i)].gtype.mutateGtype();
}

GaIndividual GaGenaration::decideParentRoulette(std::vector<GaIndividual>& candidate) {
	double r = (double)rand() / RAND_MAX; // 0~1
	std::vector<double> fitnessCategory;
	double totalFitness;
	int numOfCandidate = static_cast<int>(candidate.size());

	// calc fitness category and total fitness
	fitnessCategory.push_back(0.0);
	for (int i = 0; i < numOfCandidate; i++)
		fitnessCategory.push_back(fitnessCategory[i] + genes[i].getFitness());
	totalFitness = fitnessCategory[numOfCandidate];

	// decide parent
	for (int i = 0; i < numOfCandidate; i++)
		if ((r >= fitnessCategory[i] / totalFitness) && (r < fitnessCategory[i + 1] / totalFitness))
			return candidate[i];

	return candidate[numOfCandidate-1];
}

GaIndividual GaGenaration::getChildCrossOver(GaIndividual& parent1, GaIndividual& parent2) {
	GaIndividual child;
	child = parent1;

	// set parent data
	child.setParent1(parent1.getRank());
	child.setParent2(parent2.getRank());

	// decide intersection ramdomely
	int intersection;
	intersection = (int)(parent1.gtype.getGCodeLength() * rand() / RAND_MAX);
	child.setCrossPoint(intersection);

	// calc new gtype
	child.gtype.crossoverGtype(parent2.gtype, intersection);

	return child;
}

#endif // GA_UTIL_H