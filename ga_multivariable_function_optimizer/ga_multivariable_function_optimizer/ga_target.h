#ifndef GA_TARGET_H
#define GA_TARGET_H

#include <cmath>

/*
 * GA target params configuration
 */
struct GaConfiguration {
	int numOfVariable;
	int gtypeLength;
	int gtypeMax;
	int population;
	std::vector<double> ptypeMax;
	std::vector<double> ptypeMin;
};

/*
 * GA Farget Function
 */
double fx(std::vector<double> x) {
	double res = 0.0;

	for (int i = 0; i < static_cast<int>(x.size()); i++)
		res += (x[i] - 1.0) * (x[i] - 1.0);

	return res;
}

// GA Fitness Function
double gy(double y) {
	return 1.0 / (1.0 + fabs(y));
}

#endif // GA_TARGET_H