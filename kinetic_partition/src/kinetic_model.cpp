#include "kinetic_model.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include "svg.h"
#include "gdal_priv.h"
#include "cpl_conv.h"

Kinetic_Model::Kinetic_Model()
{
	GDALAllRegister();

    params = new Parameters();

    basename = "";
	time_string = "";

    I = Matrix<uchar>();
    I_grad_m = Matrix<double>();
    I_grad_t = Matrix<double>();
	used = Matrix<bool>();

	label_id_to_class = map<Label_ID, int>();
	I_lab = Matrix<uchar>();
	I_prob = Matrix<double>();
	I_lab_grad_m = Matrix<double>();
	I_lab_grad_t = Matrix<double>();

    I_data_size = 0;
    I_data = NULL;

    segments = vector<Segment *>();
    tree = NULL;

    rays = vector<Ray *>();
    schedule = NULL;
    graph = NULL;
	initial_graph = NULL;

    elapsed_time_building_graph = 0;
	elapsed_time_grouping = 0;
}


Kinetic_Model::Kinetic_Model(const Kinetic_Model & m)
{
    generator = m.generator;
    params = m.params;

    basename = m.basename;
	time_string = m.time_string;

    I = m.I;
	I_grad_m = m.I_grad_m;
	I_grad_t = m.I_grad_t;
	used = m.used;
	label_id_to_class = m.label_id_to_class;
	I_lab = m.I_lab;
	I_prob = m.I_prob;
	I_lab_grad_m = m.I_lab_grad_m;;
	I_lab_grad_t = m.I_lab_grad_t;

    I_data_size = m.I_data_size;
    I_data = m.I_data;

    segments = std::vector<Segment *>(m.segments.begin(), m.segments.end());
    tree = m.tree;

    rays = vector<Ray *>(m.rays.begin(), m.rays.end());
    schedule = m.schedule;
    graph = m.graph;
	initial_graph = m.initial_graph;
}


Kinetic_Model::~Kinetic_Model()
{
	delete params;
    clear();
}


void Kinetic_Model::add_line(list<LineItem *> &container, double x1, double y1, double x2, double y2, uchar r, uchar g, uchar b)
{
	LineItem* line = new LineItem(x1, y1, x2, y2, r, g, b);
    container.push_back(line);
}


void Kinetic_Model::set_path_input_image(const std::string & _path)
{
    params->path_input_image = _path;
}

void Kinetic_Model::set_max_label()
{
	params->max_label = I_prob.channels;
}

void Kinetic_Model::set_path_input_prediction(const std::string & _path)
{
	params->path_input_prediction = _path;
}

void Kinetic_Model::set_min_length(double proportion) {
	params->split_params.length_thresh = proportion*sqrt(I.rows*I.rows + I.cols*I.cols);
}

void Kinetic_Model::set_time_string()
{
	// Get time identifier
    time_t rawtime;
    struct tm * timeinfo;
    char _time_string [18];
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (_time_string, 18, "%Y-%m-%d-%H%M%S", timeinfo);
	time_string = std::string(_time_string);
}

void Kinetic_Model::set_angle_hist_index() {
	bin_to_angle_index = vector<int>(90, -1);
	double alpha_bin_width = 180. / double(bin_to_angle_index.size());
	int i = 0;
	for (map<double, Node_Parallel_Segments*>::iterator it_m = tree->parallel_segments.begin() ; it_m != tree->parallel_segments.end() ; it_m++) {
		double theta = it_m->first;
		reg_angles[i++] = theta;
	}

	if (!reg_angles.empty()) {
		// Look for the closest regularized angle
		for (int i = 0; i < bin_to_angle_index.size(); ++i) {
			double alpha = i*alpha_bin_width;
			double min_dist = 180;
			int closest_angle_index = 0;
			for (map<int, double>::iterator it_m = reg_angles.begin(); it_m != reg_angles.end(); it_m++) {
				for (int k = -1; k <= 1; k++) {
					double dist = fabs(it_m->second - alpha + k * 180);
					if (dist < min_dist) {
						min_dist = dist;
						closest_angle_index = it_m->first;
					}
				}
			}
			if (min_dist < 2.5) bin_to_angle_index[i] = closest_angle_index;
		}
	}
}


void Kinetic_Model::set_basename()
{
    boost::filesystem::path path_processed_image (params->path_input_image);
    basename = path_processed_image.stem().string();
}


void Kinetic_Model::set_gradient_maps()
{
    I_grad_m = Matrix<double>(I.rows, I.cols, 1);
    I_grad_t = Matrix<double>(I.rows, I.cols, 1);

    // 1. Smoothes I and computes luminance at the same time
    Matrix<double> I_smoothed = Matrix<double>(I.rows, I.cols, 3);
    int radius = 2; // sigma 0.95 is best for removing 2-staircases
    Matrix<double> G;
    Matrix<double>::gaussian_matrix(G, radius, 0.95);

    for (uint i = 0 ; i < I.rows ; i++) {
        for (uint j = 0 ; j < I.cols ; j++) {
            double s_r = 0, s_g = 0, s_b = 0;
            for (int k = -radius ; k <= radius ; k++) {
                for (int l = -radius ; l <= radius ; l++) {
                    double g_kl = G(k + radius, l + radius);
                    uint i_kl = uint(jclamp(0, int(i) + k, int(I.rows) - 1));
                    uint j_kl = uint(jclamp(0, int(j) + l, int(I.cols) - 1));
                    s_r += g_kl * I(i_kl, j_kl, 0);
                    s_g += g_kl * I(i_kl, j_kl, 1);
                    s_b += g_kl * I(i_kl, j_kl, 2);
                }
            }
            I_smoothed(i, j, 0) = jclamp(0, s_r, 255);
            I_smoothed(i, j, 1) = jclamp(0, s_g, 255);
            I_smoothed(i, j, 2) = jclamp(0, s_b, 255);
        }
    }

    // 2. Computes luminance (between 0 and 255), and Cb, Cr channels (between -127 and 128)
    Matrix<double> I_lum = Matrix<double>(I.rows, I.cols, 1);
	Matrix<double> I_cb = Matrix<double>(I.rows, I.cols, 1);
	Matrix<double> I_cr = Matrix<double>(I.rows, I.cols, 1);
    for (uint i = 0 ; i < I.rows ; i++) {
        for (uint j = 0 ; j < I.cols ; j++) {
            I_lum(i, j) = (0.2126 * I_smoothed(i, j, 0) + 0.7152 * I_smoothed(i, j, 1) + 0.0722 * I_smoothed(i, j, 2));
			I_cb(i, j) = 0.5*(I_smoothed(i, j, 2) - I_lum(i, j)) / (1 - 0.0722);
			I_cr(i, j) = 0.5*(I_smoothed(i, j, 0) - I_lum(i, j)) / (1 - 0.2126);
        }
    }
    I_smoothed.release();

    // 3 Computes gradient
	/*
	Norm 2 computation using 2x2 pixel window:
	A B
	C D
	where A is pixel (i,j).
	Then
	gx = B+D - (A+C)   horizontal difference
	gy = A+B - (C+D)   vertical difference
	*/


	/* gradient magnitude threshold */
	double rho = params->lsd_quant / sin(params->split_params.angle);

	double max_grad = 0;
	int mask_x[2][2] = { {-1, 1}, {-1, 1} };
	int mask_y[2][2] = { {1, 1}, {-1, -1} };
    for (uint i = 0 ; i < I.rows ; i++) {
        for (uint j = 0 ; j < I.cols ; j++) {
            double g_x = 0;
            double g_y = 0;
			double bg_x = 0;
			double bg_y = 0;
			double rg_x = 0;
			double rg_y = 0;

			for (int k = 0 ; k <= 1 ; k++) {
				for (int l = 0 ; l <= 1 ; l++) {
					uint i_kl = uint(jclamp(0, int(i) + k, int(I.rows) - 1));
                    uint j_kl = uint(jclamp(0, int(j) + l, int(I.cols) - 1));
                    g_x += mask_x[k][l] * I_lum(i_kl, j_kl);
                    g_y += mask_y[k][l] * I_lum(i_kl, j_kl);
					bg_x += mask_x[k][l] * I_cb(i_kl, j_kl);
					bg_y += mask_y[k][l] * I_cb(i_kl, j_kl);
					rg_x += mask_x[k][l] * I_cr(i_kl, j_kl);
					rg_y += mask_y[k][l] * I_cr(i_kl, j_kl);
				}
			}
			double norm = sqrt((g_x * g_x + g_y * g_y) / 4. + (bg_x * bg_x + bg_y * bg_y) / 4. + (rg_x * rg_x + rg_y * rg_y) / 4.);
			I_grad_m(i, j) = norm;

            if (norm < rho) {
                I_grad_t(i, j) = NOTDEF;
            } else {
				I_grad_t(i, j) = atan2(g_y + bg_y + rg_y, g_x + bg_x + rg_x);
            }
			/* look for the maximum of the gradient */
			if (norm > max_grad) max_grad = norm;
        }
    }
	/*normalize by the maximum of the gradient*/
	for (uint i = 0; i < I.rows; i++) {
		for (uint j = 0; j < I.cols; j++) {
			I_grad_m(i, j) /= max_grad;
		}
	}
	
	I_lum.release();
	I_cb.release();
	I_cr.release();
}

/*consider scaling before computing gradients to avoid staircase artifacts*/
void Kinetic_Model::set_gradient_label_maps() {
	I_lab_grad_m = Matrix<double>(I.rows, I.cols, 1);
	I_lab_grad_t = Matrix<double>(I.rows, I.cols, 1);

	// 1. Smoothes I_lab
	Matrix<double> I_lab_smoothed = Matrix<double>(I.rows, I.cols, 1);
	int radius = 2;
	Matrix<double> G;
	Matrix<double>::gaussian_matrix(G, radius, 0.95);  // 0.95 is best for removing 2-staircases
	double max_value = 0;
	for (uint i = 0; i < I.rows; i++) {
		for (uint j = 0; j < I.cols; j++) {
			double s = 0;
			for (int k = -radius; k <= radius; k++) {
				for (int l = -radius; l <= radius; l++) {
					double g_kl = G(k + radius, l + radius);
					uint i_kl = uint(jclamp(0, int(i) + k, int(I.rows) - 1));
					uint j_kl = uint(jclamp(0, int(j) + l, int(I.cols) - 1));
					s += I_lab(i_kl, j_kl) != INVALID_LABEL ? g_kl * I_lab(i_kl, j_kl) : 0;
				}
			}
			I_lab_smoothed(i, j) = s;
			if (s > max_value) max_value = s;
		}
	}
	/*set range to [0, 255]*/
	for (uint i = 0; i < I.rows; i++) {
		for (uint j = 0; j < I.cols; j++) {
			I_lab_smoothed(i, j) *= 255.0 / max_value;
		}
	}

	// 3 Computes gradient
	/*
	Norm 2 computation using 2x2 pixel window:
	A B
	C D
	where A is pixel (i,j).
	Then
	gx = B+D - (A+C)   horizontal difference
	gy = A+B - (C+D)   vertical difference
	*/

	/* gradient magnitude threshold */
	double rho = params->lsd_quant / sin(params->split_params.angle);

	double max_grad = 0;
	int mask_x[2][2] = { { -1, 1 },{ -1, 1 } };
	int mask_y[2][2] = { { 1, 1 },{ -1, -1 } };
	for (uint i = 0; i < I.rows; i++) {
		for (uint j = 0; j < I.cols; j++) {
			double g_x = 0;
			double g_y = 0;

			for (int k = 0; k <= 1; k++) {
				for (int l = 0; l <= 1; l++) {
					uint i_kl = uint(jclamp(0, int(i) + k, int(I.rows) - 1));
					uint j_kl = uint(jclamp(0, int(j) + l, int(I.cols) - 1));
					g_x += mask_x[k][l] * I_lab_smoothed(i_kl, j_kl);
					g_y += mask_y[k][l] * I_lab_smoothed(i_kl, j_kl);
				}
			}
			double norm = sqrt((g_x * g_x + g_y * g_y) / 4.);
			I_lab_grad_m(i, j) = norm;
			
			if (norm < rho) {
				I_lab_grad_t(i, j) = NOTDEF;
			}
			else {
				I_lab_grad_t(i, j) = atan2(g_y, g_x);
			}

			/* look for the maximum of the gradient */
			if (norm > max_grad) max_grad = norm;
		}
	}
	/*normalize by the maximum of the gradient*/
	for (uint i = 0; i < I.rows; i++) {
		for (uint j = 0; j < I.cols; j++) {
			I_lab_grad_m(i, j) /= max_grad;
		}
	}
	I_lab_smoothed.release();
}

void Kinetic_Model::hsv_to_rgb(double & h, double & s, double & v, uchar & r, uchar & g, uchar & b)
{
    int t_i = int(h / 60) % 6;
    double f = h / 60.0 - t_i;
    double l = v * (1 - s);
    double m = v * (1 - f * s);
    double n = v * (1 - (1 - f) * s);
    switch (t_i) {
    case 0: r = uchar(255 * v); g = uchar(255 * n); b = uchar(255 * l); break;
    case 1: r = uchar(255 * m); g = uchar(255 * v); b = uchar(255 * l); break;
    case 2: r = uchar(255 * l); g = uchar(255 * v); b = uchar(255 * n); break;
    case 3: r = uchar(255 * l); g = uchar(255 * m); b = uchar(255 * v); break;
    case 4: r = uchar(255 * n); g = uchar(255 * l); b = uchar(255 * v); break;
    case 5: r = uchar(255 * v); g = uchar(255 * l); b = uchar(255 * m); break;
    }
}


void Kinetic_Model::clear_line_items(list<LineItem *> & objects)
{
    for (list<LineItem *>::iterator it_l = objects.begin() ; it_l != objects.end() ; it_l++) {
        delete (*it_l);
    }
    objects.clear();
}


void Kinetic_Model::set_lsd_scale(double z)
{
    params->lsd_scale = z;
}

void Kinetic_Model::set_lsd_sigma_scale(double z)
{
    params->lsd_sigma_scale = z;
}

void Kinetic_Model::set_lsd_quant(double z)
{
    params->lsd_quant = z;
}

void Kinetic_Model::set_lsd_angle(double z)
{
    params->lsd_angle = z;
}

void Kinetic_Model::set_lsd_log_eps(double z)
{
    params->lsd_log_eps = z;
}

void Kinetic_Model::set_lsd_density(double z)
{
    params->lsd_density = z;
}

void Kinetic_Model::set_rega_epsilon(double z)
{
	params->rega_epsilon = z;
}

void Kinetic_Model::set_rega_ms_sigma(double z)
{
    params->rega_ms_sigma = z;
}

void Kinetic_Model::set_rega_ms_epsilon(double z)
{
    params->rega_ms_epsilon = z;
}

void Kinetic_Model::set_rega_ms_distance(int z)
{
    params->rega_ms_distance = z;
}

void Kinetic_Model::set_rega_ms_min_terms(int z)
{
    params->rega_ms_min_terms = z;
}

void Kinetic_Model::set_rega_ms_smooth_dist(double z)
{
    params->rega_ms_smooth_dist = z;
}

void Kinetic_Model::set_rega_angle_function(int z)
{
    params->rega_angle_function = z;
}

void Kinetic_Model::set_rega_angle_const(double z)
{
    params->rega_angle_const = z;
}

void Kinetic_Model::set_rega_angle_offset(double z)
{
    params->rega_angle_offset = z;
}

void Kinetic_Model::set_regp_ms_sigma(double z)
{
    params->regp_ms_sigma = z;
}

void Kinetic_Model::set_regp_ms_epsilon(double z)
{
    params->regp_ms_epsilon = z;
}

void Kinetic_Model::set_regp_ms_distx(int z)
{
    params->regp_ms_distx = z;
}

void Kinetic_Model::set_regp_ms_disty(double z)
{
    params->regp_ms_disty = z;
}

void Kinetic_Model::set_regp_ms_min_terms(int z)
{
    params->regp_ms_min_terms = z;
}

void Kinetic_Model::set_regp_ms_smooth_dist(double z)
{
    params->regp_ms_smooth_dist = z;
}

void Kinetic_Model::set_regp_trans_const(double z)
{
	params->regp_trans_const = z;
}

void Kinetic_Model::set_prop_policy(int z)
{
    params->prop_policy = z;
}

void Kinetic_Model::set_prop_ttl(int z)
{
    params->prop_ttl = z;
}

void Kinetic_Model::set_prop_distance(int z)
{
    params->prop_distance = z;
}

void Kinetic_Model::set_prop_min_edge(double z)
{
    params->prop_min_edge = z;
}

void Kinetic_Model::set_prop_range(int z)
{
    params->prop_range = double(z);
}


void Kinetic_Model::reinit()
{
    clear();
    init();
}


void Kinetic_Model::clear()
{
    if (graph != NULL) {
        delete graph;
        graph = NULL;
    }
	if (initial_graph != NULL) {
		delete initial_graph;
		initial_graph = NULL;
	}
    Ray::clear_rays(rays);
    //IndexedEvent::clear_schedule(schedule);
	schedule = NULL;

	label_id_to_class.clear();
    I.release();
	I_lab.release();
	I_prob.release();
	I_grad_m.release();
	I_grad_t.release();
	used.release();
	I_lab_grad_m.release();
	I_lab_grad_t.release();
    delete_byte_array(I_data_size, I_data);
    if (tree != NULL) {
        delete tree;
        tree = NULL;
    }
    Segment::clear_segments(segments);
}


void Kinetic_Model::switch_to_base_matrix()
{
	delete_byte_array(I_data_size, I_data);

	GDALDataset *poDataset = (GDALDataset *)GDALOpen(params->path_input_image.c_str(), GA_ReadOnly);
	if (poDataset != NULL) {
		I = Matrix<uchar>(poDataset);
	} else {
		std::cerr << "Error : GDAL couldn't load the specified file" << std::endl;
		exit(-1);
	}
	GDALClose(poDataset);
}


void Kinetic_Model::init()
{
	// Fills I_data
	GDALDataset* im_dataset = (GDALDataset *)GDALOpen(params->path_input_image.c_str(), GA_ReadOnly);
	if (im_dataset == NULL) {
		throw std::logic_error("Error : a problem occured while loading the image");
	}
	I_data_cols = im_dataset->GetRasterXSize();
	I_data_rows = im_dataset->GetRasterYSize();
	int channels = im_dataset->GetRasterCount();

	reallocate_byte_array<double>(I_data_rows * I_data_cols, I_data_size, I_data);
	if (channels == 1) {
		// Case of a grayscale image
		// We read all the rows of the image, and then we sequentially copy these rows to I_data
		GDALRasterBand* im_grayscale_band = im_dataset->GetRasterBand(1);
		double* buffer = (double*)CPLMalloc(I_data_cols * sizeof(double));
		for (int i = 0 ; i < I_data_rows ; i++) {
			im_grayscale_band->RasterIO(GF_Read, 0, i, I_data_cols, 1, buffer, I_data_cols, 1, GDT_Float64, 0, 0);
			memcpy(&(I_data[i * I_data_cols]), buffer, I_data_cols * sizeof(double));
		}
		CPLFree(buffer);

	} else if (channels == 3 || channels == 4) {
		// Case of a RGB image or RGBA image
		// It is almost the same process as before, except that we must compute the luminance
		GDALRasterBand* im_r_band = im_dataset->GetRasterBand(1);
		GDALRasterBand* im_g_band = im_dataset->GetRasterBand(2);
		GDALRasterBand* im_b_band = im_dataset->GetRasterBand(3);
		double* buffer_r = (double*)CPLMalloc(I_data_cols * sizeof(double));
		double* buffer_g = (double*)CPLMalloc(I_data_cols * sizeof(double));
		double* buffer_b = (double*)CPLMalloc(I_data_cols * sizeof(double));
		double* buffer = (double*)CPLMalloc(I_data_cols * sizeof(double));
		for (int i = 0 ; i < I_data_rows ; i++) {
			im_r_band->RasterIO(GF_Read, 0, i, I_data_cols, 1, buffer_r, I_data_cols, 1, GDT_Float64, 0, 0);
			im_g_band->RasterIO(GF_Read, 0, i, I_data_cols, 1, buffer_g, I_data_cols, 1, GDT_Float64, 0, 0);
			im_b_band->RasterIO(GF_Read, 0, i, I_data_cols, 1, buffer_b, I_data_cols, 1, GDT_Float64, 0, 0);
			for (int j = 0 ; j < I_data_cols ; j++) {
				buffer[j] = 0.2126 * buffer_r[j] + 0.7152 * buffer_g[j] + 0.0722 * buffer_b[j];
			}
			memcpy(&(I_data[i * I_data_cols]), buffer, I_data_cols * sizeof(double));
		}
		CPLFree(buffer_r);
		CPLFree(buffer_g);
		CPLFree(buffer_b);
		CPLFree(buffer);

	} else {
		throw std::logic_error("Error : the image should have one or three or four channels");
	}
	GDALClose(im_dataset);

	GDALDataset *poDataset = (GDALDataset *)GDALOpen(params->path_input_image.c_str(), GA_ReadOnly);
	if (poDataset != NULL) {
		int channels = poDataset->GetRasterCount();
		if (channels == 1) { // make I to have 3 channels
			Matrix<uchar>I_gray(poDataset);
			I = Matrix<uchar>(I_gray.rows, I_gray.cols, 3);
			for (uint i = 0; i < I.rows; i++) {
				for (uint j = 0; j < I.cols; j++) {
					for (uint c = 0; c < I.channels; c++) {
						I(i, j, c) = I_gray(i,j);
					}
				}
			}
			I_gray.release();
		}
		else {
			I = Matrix<uchar>(poDataset);
		}
	}
	else {
		std::cerr << "Error : GDAL couldn't load the specified file" << std::endl;
		exit(-1);
	}
	GDALClose(poDataset);

	set_min_length(0.02);

	used = Matrix<bool>(I.rows, I.cols, 1, false);

	set_gradient_maps();

	set_time_string();
    set_basename();

    segments = vector<Segment *>();
    tree = new Segment_Regularization_Tree();
	reg_angles = map<int, double>();
	bin_to_angle_index = vector<int>(90, -1);

    applied_regularization_angles = false;
    applied_regularization_ordinates = false;

    rays = vector<Ray *>();
    schedule = NULL;
    graph = NULL;
	initial_graph = NULL;
    elapsed_time_building_graph = 0;
	elapsed_time_grouping = 0;
}


void Kinetic_Model::set_label_maps(const set<pair<int, int>> & strokes_0, const set<pair<int, int>> & strokes_1) {
	I_lab.release();
	I_prob.release();
	label_id_to_class.clear();
	if (params->semantic_type == LABEL && !params->path_input_prediction.empty()) {
		GDALDataset *predDataset = (GDALDataset *)GDALOpen(params->path_input_prediction.c_str(), GA_ReadOnly);
		if (predDataset != NULL) {
			I_lab = Matrix<uchar>(predDataset);
			map<uchar, Label_ID> classes_added;
			uint label_id = 0;
			for (uint i = 0; i < I.rows; ++i) {
				for (uint j = 0; j < I.cols; j++) {
					auto added = classes_added.emplace(pair<uchar, Label_ID>(I_lab(i, j), label_id));
					if (added.second) {
						label_id_to_class.insert(pair<Label_ID, int>(label_id, I_lab(i, j)));
						I_lab(i, j) = label_id++; //re-assign label
					}
					else {
						I_lab(i, j) = added.first->second; //re-assign label
					}
				}
			}
			I_prob = Matrix<double>(I.rows, I.cols, label_id, 0);
			for (uint i = 0; i < I.rows; ++i) {
				for (uint j = 0; j < I.cols; j++) {
					int label = I_lab(i, j);
					if (label > 0) I_prob(i, j, label - 1) = 1;
				}
			}
		}
		else {
			std::cerr << "Error : GDAL couldn't load the specified file" << std::endl;
			exit(-1);
		}
		GDALClose(predDataset);
	}
	else if (params->semantic_type == PROBABILITY && !params->path_input_prediction.empty()) {
		GDALDataset *predDataset = (GDALDataset *)GDALOpen(params->path_input_prediction.c_str(), GA_ReadOnly);
		if (predDataset != NULL) {
			
			Matrix<uchar> I_prob_0(predDataset);
			double max_T = 255;
			I_prob = Matrix<double>(I_prob_0.rows, I_prob_0.cols, I_prob_0.channels, 0);
			double max_fg_prob = 0, min_fg_prob = max_T;
			Matrix<double> G;
			int radius = 3;
			Matrix<double>::gaussian_matrix(G, radius, 0.95);
			for (uint i = 0; i < I.rows; i++) {
				for (uint j = 0; j < I.cols; j++) {
					double foreground_prob = 0;
					for (int m = -radius; m <= radius; m++) {
						for (int l = -radius; l <= radius; l++) {
							double g_ml = G(m + radius, l + radius);
							uint i_ml = uint(jclamp(0, int(i) + m, int(I.rows) - 1));
							uint j_ml = uint(jclamp(0, int(j) + l, int(I.cols) - 1));
							for (uint k = 0; k < I_prob_0.channels; ++k) foreground_prob += g_ml * I_prob_0(i_ml, j_ml, k);
							if (foreground_prob > max_fg_prob) max_fg_prob = foreground_prob;
							if (foreground_prob < min_fg_prob) min_fg_prob = foreground_prob;
						}
					}
				}
			}
			if (max_fg_prob < 1e-5) max_fg_prob = 1; // set a minimal value in case of absence of objects

			I_lab = Matrix<uchar>(I_prob.rows, I_prob.cols, 1);
			label_id_to_class.insert(pair<Label_ID, int>(0, 0));
			for (uint k = 0; k < I_prob.channels; ++k) {
				label_id_to_class.insert(pair<Label_ID, int>(k+1, k+1));
			}
			for (uint i = 0; i < I_prob.rows; ++i) {
				for (uint j = 0; j < I_prob.cols; ++j) {
					uint max_k = 0;
					double max_p = 0, sum_p = 0;
					for (uint k = 0; k < I_prob.channels; ++k) {
						// set range to [0, 1]
						I_prob(i, j, k) = jclamp(0.0, double(I_prob_0(i, j, k) - min_fg_prob) / (max_fg_prob - min_fg_prob),1.0);
						double p = I_prob(i, j, k);
						sum_p += p;
						if (p > max_p) { max_k = k; max_p = p; }
					}
					I_lab(i, j) = max_p > (1 - sum_p) ? max_k + 1 : 0; // label id = 1 + channel index
				}
			}
			I_prob_0.release();
		}
		else {
			std::cerr << "Error : GDAL couldn't load the specified file" << std::endl;
			exit(-1);
		}
		GDALClose(predDataset);
	}

	set_max_label();
}
