#pragma once
#include <vector>
#include "defs.h"

using std::vector;


namespace Means 
{
	double mean(vector<double> & data);

	double mean(vector<double> & data, std::vector<double> & weights);

	void mean_sd(vector<double> & data, double & mean, double & sd);


	bool circular_mean_360(vector<double> & data, double & mean);

	bool circular_mean_360(vector<double> & data, double & mean, double & x_sd, double & y_sd);

	bool circular_mean_360(vector<double> & data, vector<double> & weights, double & mean);


	bool circular_mean_180(vector<double> & data, double & mean);

	bool circular_mean_180(vector<double> & data, vector<double> & weights, double & mean);
}