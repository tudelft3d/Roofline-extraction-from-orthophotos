#ifndef LSD_WORKER_H
#define LSD_WORKER_H

#include "kinetic_model.h"

namespace LSD 
{
    void find_segments(Kinetic_Model* model);

	void find_regions(int n, int* & reg_img, int reg_x, int reg_y, std::vector<std::list<std::pair<int, int> > > & regions);

	void compute_reference_gradient_from_lsd(int n, double* & lsd_modulus, double* & lsd_angles, int reg_x, std::vector<std::list<std::pair<int, int> > > & regions, std::vector<double> & means, std::vector<double> & angles);
#if NOT_MEASURING_PERFORMANCES
	void compute_reference_gradient(Matrix<double> & I_grad_m, Matrix<double> & I_grad_t, double lsd_x1, double lsd_y1, double lsd_x2, double lsd_y2, double & m_mean, double & t_mean);
#endif

	void compute_true_coordinates(double lsd_x1, double lsd_y1, double lsd_x2, double lsd_y2, double _width, int rows, double & x1, double & x2, double & y1, double & y2);

#if NOT_MEASURING_PERFORMANCES
    void make_layer(Kinetic_Model* model);
#endif
}

#endif // LSD_WORKER_H
