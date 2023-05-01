#include "means.h"

using std::vector;


double Means::mean(vector<double> & data)
{
	double m = 0;
	for (uint i = 0; i < data.size(); i++) m += data[i];
	m /= data.size();
	return m;
}



double Means::mean(vector<double> & data, vector<double> & weights)
{
	double m = 0, sum_weights = 0;
	for (uint i = 0; i < data.size(); i++) {
		m += weights[i] * data[i];
		sum_weights += weights[i];
	}
	m /= sum_weights;
	return m;
}



void Means::mean_sd(vector<double> & data, double & mean, double & sd)
{
	mean = sd = 0;
	for (uint i = 0; i < data.size(); i++) mean += data[i];
	mean /= data.size();

	for (uint i = 0; i < data.size(); i++) {
		double d = (data[i] - mean);
		sd += d * d;
	}
	sd = sqrt(sd / data.size());
}



bool Means::circular_mean_360(vector<double> & data, double & mean)
{
	mean = 0;
	uint n = data.size();
	vector<double> x(n, 0);
	vector<double> y(n, 0);
	for (uint i = 0; i < data.size(); i++) {
		x[i] = cos(data[i]);
		y[i] = sin(data[i]);
	}
	double x_mean = Means::mean(x);
	double y_mean = Means::mean(y);

	if (x_mean == 0 && y_mean == 0) {
		return false;
	} else {
		mean = 180 * atan2(y_mean, x_mean) / PI;
		return true;
	}
}



bool Means::circular_mean_360(vector<double> & data, double & mean, double & x_sd, double & y_sd)
{
	mean = 0;
	uint n = data.size();
	vector<double> x(n, 0);
	vector<double> y(n, 0);
	for (uint i = 0; i < data.size(); i++) {
		x[i] = cos(data[i]);
		y[i] = sin(data[i]);
	}
	double x_mean = Means::mean(x);
	double y_mean = Means::mean(y);

	if (x_mean == 0 && y_mean == 0) return false;

	mean = atan2(y_mean, x_mean);
	
	vector<double> x_rot(n, 0), y_rot(n, 0);
	for (uint i = 0 ; i < data.size() ; i++) {
		x_rot[i] = cos(-mean) * x[i] - sin(-mean) * y[i];
		y_rot[i] = sin(-mean) * x[i] + cos(-mean) * y[i];
	}
	double x_rot_mean, x_rot_sd, y_rot_mean, y_rot_sd;
	Means::mean_sd(x_rot, x_rot_mean, x_rot_sd);
	Means::mean_sd(y_rot, y_rot_mean, y_rot_sd);

	mean = 180 * mean / PI;
	x_sd = x_rot_sd;
	y_sd = y_rot_sd;

	return true;
}



bool Means::circular_mean_360(vector<double> & data, vector<double> & weights, double & mean)
{
	mean = 0;
	uint n = data.size();
	vector<double> x(n, 0);
	vector<double> y(n, 0);
	for (uint i = 0; i < data.size(); i++) {
		x[i] = cos(data[i]);
		y[i] = sin(data[i]);
	}
	double x_mean = Means::mean(x, weights);
	double y_mean = Means::mean(y, weights);

	if (x_mean == 0 && y_mean == 0) {
		return false;
	} else {
		mean = 180 * atan2(y_mean, x_mean) / PI;
		return true;
	}
}



bool Means::circular_mean_180(vector<double> & data, double & mean)
{
	mean = 0;
	uint n = data.size();
	vector<double> x(n, 0);
	vector<double> y(n, 0);
	for (uint i = 0; i < data.size(); i++) {
		x[i] = cos(2 * data[i]);
		y[i] = sin(2 * data[i]);
	}
	double x_mean = Means::mean(x);
	double y_mean = Means::mean(y);

	if (x_mean == 0 && y_mean == 0) {
		return false;
	} else {
		mean = 180 * atan2(y_mean, x_mean) / PI;
		if (mean < 0) {
			mean = 180 + mean / 2;
		} else {
			mean /= 2;
		}
		return true;
	}
}



bool Means::circular_mean_180(vector<double> & data, vector<double> & weights, double & mean)
{
	mean = 0;
	uint n = data.size();
	vector<double> x(n, 0);
	vector<double> y(n, 0);
	for (uint i = 0; i < data.size(); i++) {
		x[i] = cos(2 * data[i]);
		y[i] = sin(2 * data[i]);
	}
	double x_mean = Means::mean(x, weights);
	double y_mean = Means::mean(y, weights);

	if (x_mean == 0 && y_mean == 0) {
		return false;
	} else {
		mean = 180 * atan2(y_mean, x_mean) / PI;
		if (mean < 0) {
			mean = 180 + mean / 2;
		} else {
			mean /= 2;
		}
		return true;
	}
}