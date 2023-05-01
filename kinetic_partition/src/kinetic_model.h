#ifndef MODEL_H
#define MODEL_H

#include <random>
#include <vector>
#include <opencv2/core/core.hpp>
#include <set>

using std::default_random_engine;
using std::vector;
using std::set;
using cv::Mat;

#include "parameters.h"
#include "matrix.h"
#include "segment_ray.h"
#include "segment_tree.h"
#include "indexed_event.h"
#include "partition.h"
#include "line_item.h"

class Kinetic_Model
{
public:
    Kinetic_Model();

    Kinetic_Model(const Kinetic_Model & m);

    virtual ~Kinetic_Model();

    void set_path_input_image(const std::string & _path);

	void set_path_input_prediction(const std::string & _path);

	void switch_to_base_matrix();

	void set_max_label();

	void set_min_length(double proportion);

	void set_angle_hist_index();

	void set_time_string();

	void set_basename();

    void reinit();

	void set_label_maps(const set<pair<int, int>> & strokes_1, const set<pair<int, int>> & strokes_2);

	void add_line(list<LineItem *> & container, double x1, double y1, double x2, double y2, uchar r, uchar g, uchar b);

	void clear_line_items(list<LineItem *> & objects);

private:
    void clear();

    void init();

	void create_clusters();

	void set_gradient_maps();

	void set_gradient_label_maps();

	void hsv_to_rgb(double & h, double & s, double & v, uchar & r, uchar & g, uchar & b);

public:
    void set_lsd_scale(double z);
    void set_lsd_sigma_scale(double z);
    void set_lsd_quant(double z);
    void set_lsd_angle(double z);
    void set_lsd_log_eps(double z);
    void set_lsd_density(double z);

	void set_rega_epsilon(double z);
    void set_rega_ms_sigma(double z);
    void set_rega_ms_epsilon(double z);
    void set_rega_ms_distance(int z);
    void set_rega_ms_min_terms(int z);
    void set_rega_ms_smooth_dist(double z);

	void set_rega_angle_function(int z);
    void set_rega_angle_const(double z);
    void set_rega_angle_offset(double z);

    void set_regp_ms_sigma(double z);
    void set_regp_ms_epsilon(double z);
    void set_regp_ms_distx(int z);
    void set_regp_ms_disty(double z);
    void set_regp_ms_min_terms(int z);
    void set_regp_ms_smooth_dist(double z);
	void set_regp_trans_function(int z);
    void set_regp_trans_const(double z);

    void set_prop_policy(int z);
    void set_prop_ttl(int z);
    void set_prop_distance(int z);
    void set_prop_min_edge(double z);
    void set_prop_range(int z);

	void set_prop_extra_enabled(bool z);

	void set_prop_region_length(double z);
	void set_prop_region_width(double z);
	void set_prop_sub_region_length(double z);
	void set_prop_sub_region_width(double z);
	void set_prop_compared_quantity(int z);
	void set_prop_ratio(double z);
	void set_prop_dot_th(double z);

	void set_prop_check_m_enabled(bool z);
	void set_prop_m_compare_to(int z);
	void set_prop_m_factor(double z);
	void set_prop_m_fixed_magnitude(double z);

	void set_prop_check_t_enabled(bool z);
	void set_prop_t_compare_to(int z);
	void set_prop_t_tolerance(double z);

private:
    template <typename T>
    void reallocate_byte_array(uint size, uint & byte_array_size, T* & byte_array)
    {
        if (size != byte_array_size) {
            if (byte_array != NULL) delete[] byte_array;
            byte_array_size = size;
            byte_array = new T[byte_array_size];
        }
    }

    template <typename T>
    void delete_byte_array(uint & byte_array_size, T* & byte_array)
    {
        if (byte_array_size != 0) {
            delete[] byte_array;
            byte_array = NULL;
            byte_array_size = 0;
        }
    }

public:
    default_random_engine generator;
    Parameters* params;

    /** Level 0 : user inputs */

    std::string basename;
	std::string time_string;

    /** Level 1 : LSD */

    Matrix<uchar> I;
    Matrix<double> I_grad_m;
    Matrix<double> I_grad_t;
	Matrix<bool> used;

    uint I_data_size;
    double* I_data; //luminance, in [0, 255]
	int I_data_rows;
	int I_data_cols;

    vector<Segment *> segments;
    Segment_Regularization_Tree* tree;
	map<int, double> reg_angles;
	vector<int> bin_to_angle_index;

    /** Level 2 : Regularization (1) */

    bool applied_regularization_angles;

    /** Level 3 : Regularization (2) */

    bool applied_regularization_ordinates;

    /** Level 4 : Propagation */

    vector<Ray *> rays;
    IndexedEvent* schedule;
    Partition* graph;
	Partition* initial_graph;
    list<LineItem *> L_prop;
    list<LineItem *> L_next_operation;
    list<LineItem *> L_split_segments;

    /** Level 5 : Evolving cells */

	map<Label_ID, int> label_id_to_class; // 0 for background, classes are re-assigned labels starting from 1, and +1 for each class that exists in the image
	Matrix<uchar> I_lab; // label matrix. 
	Matrix<double> I_prob; // probability matrix. One channel for each non-background class, i.e. channel i corresponds to label i+1
	Matrix<double> I_lab_grad_m;
	Matrix<double> I_lab_grad_t;

	double elapsed_time_building_graph;
	double elapsed_time_grouping;
};

#endif // MODEL_H
