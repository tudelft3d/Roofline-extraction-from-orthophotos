#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>

typedef enum {LABEL, PROBABILITY, NONE} Semantic_Type;

struct Split_Params {
	double m_value_thresh;
	double lab_m_value_thresh;
	double distance_to_edge;
	double bin_width;
	double length_thresh;
	double gradient_percentage;
	double angle;
	double min_dist_angle_score;
};

class Parameters
{
public:
    Parameters();

    ~Parameters();

    void reset();

	void reset_lsd();

	void reset_rega();

	void reset_regp();

	void reset_prop();

	void reset_split();

public:
    /** Level 0 : user inputs */

	int verbose_level;

    std::string path_input_image;
	std::string path_input_directory;
	std::string path_input_prediction;
	std::string path_input_prediction_directory;
	std::string path_input_partition;

	std::string path_output_directory;
	std::string suffix;

    /** Level 1 : LSD parameters */

    double lsd_scale;
    double lsd_sigma_scale;
	double lsd_quant;
    double lsd_angle;
    double lsd_log_eps;
    double lsd_density;

    /** Level 2 : Regularization parameters (1) */

    bool rega_regp_enabled;

	double rega_epsilon;
    double rega_ms_sigma;
    double rega_ms_epsilon;
    int rega_ms_distance;
    int rega_ms_min_terms;
    double rega_ms_smooth_dist;

    int rega_angle_function;
    double rega_angle_const;
    double rega_angle_offset;

    /** Level 3 : Regularization parameters (2) */

    double regp_ms_sigma;
    double regp_ms_epsilon;
    int regp_ms_distx;
    double regp_ms_disty;
    int regp_ms_min_terms;
    double regp_ms_smooth_dist;

	double regp_trans_const;

    /** Level 4 : Propagation */
    int prop_K;
    int prop_policy;
    int prop_ttl;
    int prop_distance;
    double prop_min_edge;
	double prop_corner_eps;
    double prop_range;

	bool prop_check_m_enabled;
	int prop_m_compare_to;
	double prop_m_factor;
	double prop_m_fixed_magnitude;

	bool prop_check_t_enabled;
	int prop_t_compare_to;
	double prop_t_tolerance;

	/** Level 5 : Segmentation */
	Semantic_Type semantic_type;
	int max_iters;
	Split_Params split_params;
	int max_label;
	double fidelity_beta;
	double prior_lambda;
};

#endif // PARAMETERS_H
