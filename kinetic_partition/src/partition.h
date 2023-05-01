#pragma once
#include "parameters.h"
#include "partition_elements.h"
#include "indexed_event.h"
#include "quadtree.h"
#include "matrix.h"
#include <set>
#include <random>
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

using std::list;
using std::pair;
using std::vector;
using std::set;

using cv::Vec3b;
using cv::Vec3d;

typedef enum { MERGE, SPLIT, REMOVE, STOP } Operator_Type;


struct Boundary_Edge_Info {
	Edge* e;
	Label_ID label;
	Boundary_Edge_Info(Edge* e_, Label_ID label_) : e(e_), label(label_) {}
};

struct Vertex_Property {
	Face* f;
	double t;
};

struct Edge_Property {
	double t;
	Label_ID merged_label;
};

typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, Vertex_Property, Edge_Property> Dual_Graph;


struct Merge_Graph_Vertex_Property {
	list<Face*> superfacet;
	list<Boundary_Edge_Info> boundary_edges; /*edges that overlap with support facet*/
	size_t n_edges;
	int label;
	vector<double> semantic_probabilities;
	int n_pixels;
	double perimeter;
};

struct Merge_Graph_Edge_Property {
	list<Edge *> shared_edges;
	double t = 0.;
	Label_ID merged_label;
};

typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, Merge_Graph_Vertex_Property, Merge_Graph_Edge_Property> Dual_Merge_Graph;

class Partition
{
public:
	Partition();

    Partition(uint _rows, uint _cols, Parameters* _params);

	Partition(const Partition & P);

	Partition & operator= (const Partition & P);

	~Partition();
	
private:
	void clear();
	
	void copy(const Partition & P);

public:
	int & get_id_vertices();

	void init_edges_polygon_boundary(vector<list<HalfEdge *>> & edges, vector<Ray *> & rays, map<HalfEdge *, list<Segment *>> & boundary_collinear_segments, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
		vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges, vector<Vertex *> & initial_vertices_of_rays, vector<Vertex *> & initial_outer_vertices);

	void build_edge(vector<Ray *> & rays, IndexedEvent* current_event, bool & r_i_continue, Vertex* intersection_point, bool vertex_is_corner, bool vertex_is_new, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
		vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges, vector<Vertex *> & initial_vertices_of_rays);

	void build_faces();

	void merge_containers(list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices, vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges);

	void build_merge_graph(Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t);

	void build_dual_graph(Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t);

	void naive_labeling();

	pair<double, double> highest_merge_split_gains();

	Operator_Type realize_next_operation(Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, 
		Matrix<double> & I_m, Matrix<double> & I_t, Matrix<double> & I_lab_m, Matrix<double> & I_lab_t, map<int, double> & angles, vector<int> & bin_to_angle_index);

	void mark_edges_next_operation();
	
	void find_pixels_inside_all_facets();

	void compute_feature_all_facets(Matrix<double> & I_prob, Matrix<uchar> & I);

	void weight_all_edges(Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t);

	pair<double, double> joint_grad_cost_and_length(const list<Edge *> & edges);

	void find_active_gradients_all_facets(Matrix<double> & I_m, Matrix<double> & I_t, bool is_color_image);

	void detect_line_segments_all_facets(Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t, bool is_color_image);

	void estimate_cdf_grad_m(Matrix<double> & I_m, bool is_color_image);

	void extend_detected_segments_all_facets(Matrix<uchar> & I, Matrix<double> & I_prob,
		Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t, map<int, double> & angles, vector<int> & bin_to_angle_index);

	void remove_bivalent_vertiex(Vertex* v, Matrix<double> & I_m, Matrix<double> & I_t, Matrix<double> & I_lab_m, Matrix<double> & I_lab_t, 
		Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, map<int, double> & angles, vector<int> & bin_to_angle_index);

	Edge* merge_collinear_edges(Vertex* v, list<Vertex*>& vertices_erased, bool destroy=true);

	void seal_and_remove_bivalent_vertices(int n_boundary_vertices);

	void reset_indices();

	bool check_removal_criteria(Vertex* v, Matrix<double> & I_m, const Matrix<double> & I_t, Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used);

	vector<Face *> split_facet(Face* f, Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t, Matrix<double> & I_lab_m, Matrix<double> & I_lab_t);

	void compute_energy();

	void save_svg(const std::string & directory, const std::string &basename, bool show_details, std::string suffix = "");

	void draw_faces_svg(const std::string & directory, const std::string & fileNameWithoutExtension, list<Face *> f_list);

	void subpartition_info_to_svg(const std::string & directory, const vector<Segment *> & segments, Face* f = nullptr);

	void save_ply(const std::string & directory, const std::string & basename);

	int get_n_iterations() const {return n_iterations;}

private:
	Face* third_facet(Vertex* v_i, HalfEdge* h_ij, Face* f_i, Face* f_j);
	
	void propose_next_operation();

	void split_operator(Matrix<uchar> & I, Matrix<double> & I_prob,
		Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t, Matrix<double> & I_lab_m, Matrix<double> & I_lab_t, 
		vector<Dual_Graph::vertex_descriptor> & graph_vertices, map<int, double> & angles, vector<int> & bin_to_angle_index);

	void merge_gain(const list<Dual_Merge_Graph::vertex_descriptor> & face_vertices_list, const map<list<Edge *>, bool>& shared_edges, 
	Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t, double & gain, Label_ID & merged_label, bool ignore_boundary_conditions, bool verbose=false);

	void merge_gain(Face* f_1, Face* f_2, Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t, 
		double & gain, Label_ID & merged_label, bool verbose = false);

	void merge_gain(const list<vector<double>>& faces_prob_hists, const list<vector<Boundary_Edge_Info>> & adj_face_edges_1ist, 
		const list<int>& labels, const list<int>& num_pixels, const map<list<Edge *>, bool>& shared_edges, 
		double & gain, Label_ID & merged_label, bool ignore_boundary_conditions, bool verbose);

	int detect_line_segments(Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t, Face* f, bool is_color_image);

	void extend_detected_segments(Face* f, Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t, 
		map<int, double> & angles, vector<int> & bin_to_angle_index, double & gain);
	
	void merge_operator(Matrix<uchar> & I, Matrix<double> & I_prob,
		Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t, 
		Matrix<double> & I_lab_m, Matrix<double> & I_lab_t, vector<Dual_Graph::vertex_descriptor> & vi_vj, map<int, double> & angles, vector<int> & bin_to_angle_index);

	Face* merge_two_facets(map<int, list<Face *>::iterator> & map_facets, Face* f_i, Face* f_j, 
		Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t, Label_ID& label);

	bool loop_inside_sub_facet(HalfEdge* h_e, vector<list<HalfEdge*>> & sub_facet, vector<list<HalfEdge*>> & parent_facet);

	void build_split_faces(const set<Edge*>& split_edges, Face* parent_face, vector<Face *> & split_facets);

	void split_edge_and_update_facets(HalfEdge* h, Vertex* v, bool destroy = TRUE);

public:
    //void fill_matrix_of_iterators();

    int count_nonsimple_faces();

    void save_edges(Matrix<uchar> & I, const std::string & directory, const std::string & basename, double f = 1.0);

	void save_labels(const std::string &directory, const std::string &basename, const map<Label_ID, int> & label_id_to_class);

	void save_faces(const std::string & directory, const std::string & basename);

    void save_graph_definition(const std::string & directory, const std::string & basename);

    void save_boundaries(const std::string & directory, const std::string & basename);

private:
	Vertex* erase(Vertex* v);

	Edge* erase(Edge* e, bool destroy);

	void erase(list<Edge *> & l_e, bool destroy);

	void push_back(Vertex* v);

	void push_back(Edge* e);

public:

	void draw_highlighted_facets(Matrix<uchar> & J);

	void draw_active_gradients(Matrix<uchar> & J, double alpha = 0.6);

    void draw_edges(Matrix<uchar> & I, Matrix<uchar> & J, double f = 1.0);

	void draw_labels(Matrix<uchar> & J, const map<Label_ID, int> & label_id_to_class, bool color);

	void draw_faces(Matrix<uchar> & J);

	void draw_intermediate_graph(Matrix<uchar> & I, Matrix<uchar> & J, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
		vector<list<Edge *> > & outer_edges, vector<list<Edge *> > & inner_edges);

public:
	int rows;
	int cols;
	int area;
	double diag_length;

	// pixel to vertex coordinate offset (gradient is computed on the pixel's bottom-right corner)
	Point2d offset = Point2d(0.5, -0.5);

    Parameters* params;

	Quadtree* quadtree;

	int id_vertices;
	int id_edges;
	int id_faces;
	int id_segments;

	Vertex* vertices_head;
	Vertex* vertices_tail;
	uint v_size;

	Edge* edges_head;
	Edge* edges_tail;
	uint e_size;

	list<Face *> faces;

	vector<double> cdf_grad_m;
	vector<double> cdf_lab_grad_m;

	Operator_Type next_op;
	vector<Dual_Graph::vertex_descriptor> vs; // graph vertices to be modified in the next operation
	
	// interface
	vector<vector<pair<Point2d, Point2d>>> edges_next_operation;
	vector<Face *> faces_next_operation;
	double energy;
	std::default_random_engine generator;

private:
	set<pair<Dual_Graph::edge_descriptor, double>, queuecomp> queue_merge;
	set<pair<Dual_Graph::vertex_descriptor, double>, queuecomp> queue_split;
	Dual_Graph dual_graph;
	map<Face*, Dual_Graph::vertex_descriptor> face_to_graph;

	set< pair<Vertex*, double> , bivalentcomp > queue_bivalent; // from highest to lowest
	set< pair<Vertex*, double>, bboxcomp_by_x> bboxes_by_x; // from highest to lowest
	set< pair<Vertex*, double>, bboxcomp_by_y> bboxes_by_y; // from highest to lowest

	set<pair<Dual_Merge_Graph::edge_descriptor, double>, queuecomp> queue_merge_graph;
	Dual_Merge_Graph merge_graph;
	int n_iterations;
};
