#include "parameters.h"

Parameters::Parameters()
{
	reset();
}


Parameters::~Parameters()
{

}

void Parameters::reset()
{
	verbose_level = 3;
	path_input_image = "";
	path_input_directory = "";

	path_input_prediction = "";
	path_input_prediction_directory = "";

	path_input_partition = "";

	path_output_directory = ".";
	suffix = "";

	// LSD default parameters
	reset_lsd();

	// Regularization default parameters (1)
	reset_rega();

	// Regularization parameters (2)
	reset_regp();

	// Propagation parameters
	reset_prop();

	// Merge and split default parameters
	semantic_type = NONE;
	max_iters = 5000;
	max_label = 0;
	fidelity_beta = 0.0002;
	prior_lambda = 0.000001;

	// Segment detection parameters
	reset_split();
}

void Parameters::reset_lsd()
{
	lsd_scale = 0.8;
	lsd_sigma_scale = 0.6;
	lsd_quant = 2.0;
	lsd_angle = 22.5;
	lsd_log_eps = -3;
	lsd_density = 0.7;
}


void Parameters::reset_rega()
{
	rega_regp_enabled = true;
	
	rega_epsilon = 1.5;

	rega_ms_sigma = 4;
	rega_ms_epsilon = 0.25;
	rega_ms_distance = 200;
	rega_ms_min_terms = 6;
	rega_ms_smooth_dist = 3;

	rega_angle_function = 0;
	rega_angle_const = 3;
	rega_angle_offset = 1.8;
}


void Parameters::reset_regp()
{
	regp_ms_sigma = 4;
	regp_ms_epsilon = 0.25;
	regp_ms_min_terms = 2;
	regp_ms_distx = 4500;
	regp_ms_disty = 10;
	regp_ms_smooth_dist = 3;

	regp_trans_const = 1.5;
}

void Parameters::reset_prop()
{
    prop_K = 2;
	prop_policy = 1;
	prop_ttl = 1;
	prop_distance = 50;
	prop_min_edge = 1e-6;
	prop_corner_eps = 1e-6;
	prop_range = 50;
}

void Parameters::reset_split()
{
	split_params.m_value_thresh = 0;
	split_params.lab_m_value_thresh = 0;
	split_params.distance_to_edge = 2.9;
	split_params.bin_width = 0.002;
	split_params.length_thresh = 30.0; // is reset by kinetic_model init()
	split_params.gradient_percentage = 0.3;
	split_params.angle = 0.3926; // +-22.5 degrees
	split_params.min_dist_angle_score = 0.5; // half maximum of the two normal distributions
}
