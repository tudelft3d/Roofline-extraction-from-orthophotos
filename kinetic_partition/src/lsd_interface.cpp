#include "lsd_interface.h"
#include "lsd.h"
#include "segment_ray.h"
#include "means.h"
#include "trace.h"



void LSD::find_regions(int n, int* & reg_img, int reg_x, int reg_y, vector<list<pair<int, int> > > & regions)
{
	// Initializes a vector of n lists
	regions = vector<list<pair<int, int> > >(n, list<pair<int, int> >());

	// Loops on the map of regions
	// When it meets a pixel that belongs to a certain area r (r > 1),
	// we add its coordinates to the list found in entry r - 1
	for (int i = 0 ; i < reg_y ; i++) {
		for (int j = 0 ; j < reg_x ; j++) {
			int r = reg_img[j + i * reg_x];
			if (r > 0) {
				regions[r - 1].push_back(std::make_pair(i, j));
			}
		}
	}
}


void LSD::compute_reference_gradient_from_lsd(int n, double* & lsd_modulus, double* & lsd_angles, int reg_x, vector<list<pair<int, int> > > & regions, vector<double> & means, vector<double> & angles)
{
	means = vector<double>(n, 0);
	angles = vector<double>(n, 0);
	for (int r = 0 ; r < n ; r++) {

		vector<double> v_modulus;
		vector<double> v_theta;
		v_modulus.reserve(regions[r].size());
		v_theta.reserve(regions[r].size());

		for (list<pair<int, int> >::iterator it_p = regions[r].begin() ; it_p != regions[r].end() ; it_p++) {
			int i = it_p->first;
			int j = it_p->second;
			v_modulus.push_back(lsd_modulus[j + i * reg_x]);
			v_theta.push_back(lsd_angles[j + i * reg_x]);
		}

		double mean_modulus, mean_theta;
		mean_modulus = Means::mean(v_modulus);
		Means::circular_mean_360(v_theta, v_modulus, mean_theta);

		means[r] = mean_modulus;
		angles[r] = mean_theta;
	}
}


void LSD::find_segments(Kinetic_Model *model)
{
	clock_t t_begin = clock();

    int lsd_segments_size = 0;
    double* lsd = NULL;
    std::vector<Segment *> & segments = model->segments;

    if (segments.size() != 0) {
        model->tree->delete_parallel_nodes();
        Segment::clear_segments(segments);
    }

	//int *reg_img = NULL;
	//int reg_x, reg_y;
	//double* lsd_modulus = NULL;
	//double* lsd_angles = NULL;
    lsd = LineSegmentDetection(&lsd_segments_size, model->I_data,
                                        model->I_data_cols, model->I_data_rows,
                                        model->params->lsd_scale, model->params->lsd_sigma_scale, model->params->lsd_quant, model->params->lsd_angle, model->params->lsd_log_eps, model->params->lsd_density,
                                        1024, NULL, NULL, NULL);

	//free(reg_img);
	//free(lsd_modulus);
	//free(lsd_angles);

	/*model->I_grad_t_uchar = Matrix<uchar>(model->I.rows, model->I.cols, 3);
	for (uint i = 0 ; i < model->I.rows ; i++) {
		for (uint j = 0 ; j < model->I.cols ; j++) {
			if (reg_img[j + i * reg_x] > 0) {
				for (uint c = 0; c < 3; c++) model->I_grad_t_uchar(i, j, c) = 255;
			} else {
				for (uint c = 0; c < 3; c++) model->I_grad_t_uchar(i, j, c) = 0;
			}
		}
	}*/


    segments.reserve(lsd_segments_size);
    for (int i = 0 ; i < lsd_segments_size ; i++) {

		// For lsd: The coordinates origin is at the center of pixel (0,0). with y-axis pointing downwards. 
		// We work in a new coordinates system where y-axis points upwards, and the coordinates origin at the center of pixel (rows-1 , 0)
		double x1 = lsd[7 * i + 0];
		double y1 = int(model->I_data_rows)-1 - lsd[7 * i + 1];
		double x2 = lsd[7 * i + 2];
		double y2 = int(model->I_data_rows)-1 - lsd[7 * i + 3];
		double width = lsd[7 * i + 4];
		double log_nfa = lsd[7 * i + 6];

		segments.push_back(new Segment(i, x1, y1, x2, y2, width, log_nfa, false, true));
    }

	free(lsd);

	clock_t t_end = clock();
	trace(model->params->verbose_level, 5, "** Found " + std::to_string(model->segments.size()) + " segments : " + std::to_string(float(t_end - t_begin) / CLOCKS_PER_SEC) + " s.");
	model->switch_to_base_matrix();
}

void LSD::compute_true_coordinates(double lsd_x1, double lsd_y1, double lsd_x2, double lsd_y2, double lsd_width, int rows, double & x1, double & x2, double & y1, double & y2)
{
	Point2d P_1 = Point2d(lsd_x1, rows-1 - lsd_y1);
	Point2d P_2 = Point2d(lsd_x2, rows-1 - lsd_y2);
	Point2d B = (P_1 + P_2) / 2;

	// P_1 and P_2 are opposite corners of a region, identified by the LSD, that contains many pixels with similar angles
	// They divide this rectangular region, into two opposite right triangles. These triangles are made of three vertices, 
	// whose angles are PI / 2, alpha and beta. Let alpha be the smaller of the two angles.
	double dx = (P_2.x - P_1.x), dy = (P_2.y - P_1.y);
	double alpha = asin(lsd_width / sqrt(dx * dx + dy * dy));
	double beta = PI / 2 - alpha;

	// Theta is the angle made by the region, relatively to the (Ox) axis.
	double theta = atan2(P_2.y - P_1.y, P_2.x - P_1.x);
	double theta_deg = theta * 180 / PI;

	// Transforms P_1 and P_2's coordinates to the local frame, centered in B and rotated by -theta
	Point2d P_1_tr = Point2d(cos(-theta) * (P_1.x - B.x) - sin(-theta) * (P_1.y - B.y), 0.0);
	Point2d P_2_tr = Point2d(cos(-theta) * (P_2.x - B.x) - sin(-theta) * (P_2.y - B.y), 0.0);

	// Obtains the coordinates of the central line
	Point2d Q_1_tr = Point2d(P_1_tr.x - cos(beta) * lsd_width, P_1_tr.y - sin(beta));
	Point2d Q_2_tr = Point2d(P_2_tr.x + cos(beta) * lsd_width, P_2_tr.y + sin(beta));

	// Inverse transformation
	Point2d Q_1 = Point2d(B.x + cos(theta) * Q_1_tr.x - sin(theta) * Q_1_tr.y, B.y + sin(theta) * Q_1_tr.x + cos(theta) * Q_1_tr.y);
	Point2d Q_2 = Point2d(B.x + cos(theta) * Q_2_tr.x - sin(theta) * Q_2_tr.y, B.y + sin(theta) * Q_2_tr.x + cos(theta) * Q_2_tr.y);

	// Finally, obtains the desired variables
	x1 = Q_1.x, y1 = Q_1.y;
	x2 = Q_2.x, y2 = Q_2.y;
}