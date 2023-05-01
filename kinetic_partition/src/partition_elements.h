#pragma once
#include <list>
#include <vector>
#include <set>

#include "comparators.h"
#include "indexed_event.h"
#include "quadtree_point.h"
#include "segment_ray.h"
#include "parameters.h"

using cv::Point2d;
using cv::Point2i;
using std::list;
using std::vector;
using std::pair;
using std::array;
using std::set;

inline double compute_pixel_score(double & grad_score, double & angle_diff) {
	/* sigma = 10 degrees */
	double d_angle_thresh = fabs(angle_diff) - 0.087; // allow +-5 degrees
	if (d_angle_thresh < 0) d_angle_thresh = 0;
	return grad_score*exp(-d_angle_thresh*d_angle_thresh / 0.06092348);
}

inline double compute_pixel_weight(double & transverse_dist) {
	/* sigma_d = 1.0 */
	return exp(-transverse_dist*transverse_dist / 2);
}

inline double compute_unit_grad_cost(const double & unit_edge_cost, const double & mean_angle_diff, const double & mean_grad_center_offset) {
	return unit_edge_cost > 1e-4 ? unit_edge_cost : 1e-4;
}

inline float fidelity_edge(bool consistent_label, const double & unit_grad_cost, const double & length) {
	//return consistent_label ? 0 : float(unit_grad_cost) * float(length);
	return float(unit_grad_cost) * float(length);
}

class Vertex;

class Edge;

class HalfEdge;

class Face;

typedef enum { INNER_EDGE, OUTER_EDGE, ARTIFICIAL_EDGE } Edge_Type;
typedef pair<int, list<pair<int, int>>> HalfEdges_Location;
typedef int Label_ID;

struct Adj_Edge_Info {
	Edge* e;
	Label_ID label;
	int n_overlap;
	int n_nonoverlap;
};

class Edge
{

protected:
	Edge(int & _id_edge, Vertex* _v1, Vertex* _v2, double _width, double & alpha_1, double & alpha_2, bool check = true);

public:
	virtual ~Edge();

public:
	void init_vertices_and_halfedges(Vertex* _v1, double & _alpha_1, Vertex* _v2, double & _alpha_2, bool check);

	double get_alpha_1();

	double get_alpha_2();

    virtual void line_equation(double & a, double & b, double & c) = 0;

	double unit_grad_cost();

	void disable();

	void weight(Matrix<bool> & used, const Matrix<double> & I_m, const Matrix<double> & I_t, Point2d & offset, Split_Params & split_params, const vector<double> & cdf_grad_m);

	static void weight(const Point2d & A, const Point2d & B, double & width, vector<Point2i> & region, Vec2d & normal, double & grad_center_offset,
		Matrix<bool> & used, const Matrix<double> & I_m, const Matrix<double> & I_t, Point2d & offset, Split_Params & split_params, 
		const vector<double> & cdf_grad_m, double & unit_edge_cost, double & mean_angle_diff, double & grad_weight, bool check_used);

public:
	int id_edge;
	double length;
	double width; // width of the support rectangle
	double unit_edge_cost;
	double mean_angle_diff;
	double grad_center_offset;
	double grad_weight;
	bool enabled;

	Edge_Type type;
	
	Vertex* v1;
	Vertex* v2;

	Vec2d normal;

	vector<Point2i> region;

private:
	double alpha_1;
	double alpha_2;

public:
	HalfEdge* v1_v2;
	HalfEdge* v2_v1;

	Edge* e_prev;
	Edge* e_next;

	//int tag;
};


class Artificial_Edge : public Edge
{
public:
	Artificial_Edge(int & _id_edge, double _a, double _b, double _c, Vertex* _v1, Vertex* _v2, double _width, double _alpha_1, double _alpha_2);

	//Artificial_Edge(const Artificial_Edge & E);

	~Artificial_Edge();

	void line_equation(double &a, double &b, double &c);
private:
	double a;
	double b;
	double c;
};


class Inner_Edge : public Edge
{
public:
	Inner_Edge(int & _id_edge, Ray* _r, int _tag, Vertex* _v1, Vertex* _v2, double _width, double alpha_1, double alpha_2, double _t1 = 0, double _t2 = 0, bool check=true);

	Inner_Edge(int & _id_edge, const set<Ray*> & _support, int _tag, Vertex* _v1, Vertex* _v2, double _width, double alpha_1, double alpha_2, bool check = true);

	//Inner_Edge(const Inner_Edge & E);

	~Inner_Edge();

public:
	void add_supporting_ray(Ray* s);
	
	void add_supporting_rays(set<Ray *> & s);

	void get_supporting_rays(set<Ray *> & s);

	void get_supporting_segments(set<Segment *> & s);

	Ray* get_front_supporting_ray();

	Segment* get_front_supporting_segment();

	void line_equation(double &a, double &b, double &c);

	void time_range(Ray* r, double &t_1, double &t_2);

public:
	set<Ray *> rays;
	int tag; // tag = overlap_with_outer_edge ? -1 : (overlap_with_inner_edge? -2 : 0);
};



class Outer_Edge : public Edge
{
public:
	Outer_Edge(int & _id_edge, Image_Boundary _boundary, Vertex* _v1, Vertex* _v2, double _width, double alpha_1, double alpha_2, bool _v1_2, double _t1 = 0, double _t2 = 0);

	Outer_Edge(int & _id_edge, Image_Boundary _boundary, Face* _outer_face, pair<int, int> _init_boundary_location, uint _support_boundary_index, 
		Vertex* _v1, Vertex* _v2, double _width, double alpha_1, double alpha_2, bool _v1_2, Vec2d _dir_v1v2, double _t1, double _t2);

	Outer_Edge(int & _id_edge, Image_Boundary _boundary, Face* _outer_face, pair<int, int> _init_boundary_location, uint _support_boundary_index,
		Vertex* _v1, Vertex* _v2, double _width, double alpha_1, double alpha_2, bool _v1_2, double _t1, double _t2);

	~Outer_Edge();

	void line_equation(double &a, double &b, double &c);

	void time_range(double &t_1, double &t_2, const Point2d & O, int intersected_index);
public:
	Image_Boundary boundary;

	pair<int, int> init_boundary_location;

	uint support_boundary_index;

	Vec2d dir_v1v2; // we associate to the outer edge an imaginary ray that propagates from v1 to v2 that reaches v2 at time = 0.

	bool v1_2;

	Face* outer_face;
};



class HalfEdge
{
public:
	HalfEdge(Edge* _e, bool _v1_v2);

	~HalfEdge();

	void set(Face* _f = nullptr);

	HalfEdge* opposite();

	static void intersects_if_extended(HalfEdge* h, list<HalfEdge *> & intersectable_halfedges, Point2d & intersection_point, HalfEdge* & intersected);

	static int is_clockwise(list<HalfEdge *> & _edges);

public:
	Edge* e;
	bool v1_v2;
	Face* f;
};


class Vertex : public Quadtree_Point
{
public:
	Vertex(int & _id_vertex, double x, double y);

	Vertex(int & _id_vertex, IndexedEvent* _event, double x, double y);

	Vertex(int & _id_vertex, IndexedEvent* _event, Point2d & _pt);

	//Vertex(const Vertex & V);

	~Vertex();

	Point2d keypoint() { return pt; }

	double incidence_angle(const Edge* e);
	
	void add(double alpha, HalfEdge *h, bool check);

	void add(Face *f);

	void remove(HalfEdge *h);

	void remove(Face *f);

	int connectivity();

	bool outer_vertex_is_very_close(Point2d & _pt);

	bool inner_vertex_created_by_colinear_ray(Ray* r_intersectant, Ray* r_intersected, vector<Ray *> & rays, bool verbose);

	bool boundary_vertex_on_colinear_ray(Ray* r_intersectant, vector<Ray *> & rays);

	bool same_intersection_point(Ray* r_intersectant, Ray* r_intersected);

	bool has_no_duplicated_event();

	static bool right_turn(HalfEdge* h_prev, Vertex* v_curr, HalfEdge* h_next);

	static bool right_turn(Vertex* v_prev, Vertex* v_curr, Vertex* v_next);

	static bool approximate_coordinates(Ray* r_i, Ray* r_j, double t_i, int rows, int cols, Point2d & pt, double & corner_eps);

	static bool approximate_coordinates_polygon_corners(IndexedEvent* upcoming_event, Ray* r_i, Ray* r_j, double t_i, vector<list<Outer_Edge *> > & outer_edges, Point2d & pt, Vertex* & intersection_vertex, double min_edge_len, double corner_eps);

private:
	static bool round_coordinates(Point2d & pt, int rows, int cols, double eps);

public:
	int id_vertex;
	list<IndexedEvent*> events;
	Point2d pt;
	vector<pair<double, HalfEdge*> > directions; // incident angles from v, alpha from -PI to PI, small to large
	set<Face *> faces;

	double energy_gain = -10000;
	pair<double, double> bbox = pair<double, double>(0., 0.);

	Vertex* v_prev;
	Vertex* v_next;
};


typedef enum {
	UNDETERMINED = -1,
	CORNER = 0, // vertices which define a corner of the facet and not a flat angle
	NON_CORNER_TRIVALENT = 1,
	NON_CORNER_BIVALENT = 2
} Vertex_Type;


struct Split_Edge_Property {
	double length;
	double width;
	double unit_edge_cost;
	double mean_angle_diff;
	double grad_center_offset;
	double grad_weight;
	vector<Point2i> region;
	Vec2d normal;
};


class Face
{
public:
	Face(int & _id_face, vector<list<HalfEdge *>> & _edges);

	Face(int & _id_face, list<HalfEdge *> & _edges);

	~Face();

	bool is_simple(bool verbose);

	bool is_valid();

	void find_pixels_inside_facet(const Point2d & offset = Point2d(0.5, -0.5));

	void find_active_gradients(Matrix<double> & I_m, Matrix<double> & I_t, Point2d & offset, Split_Params & split_params, bool is_color_image);
	
	void discard_grad_pixels_on_edges(vector<set<Point2i, pointcomp>> & bins, Matrix<double> & I_m, double value_thresh, double bin_width);

	void delete_invalid_segments();

	void print_details();

	void compute_feature(Matrix<double> & I_prob, Matrix<uchar> & I, const Point2d & offset, int & max_label);

	void list_vertices();

private:

	void find_leftmost_and_rightmost_edges(int contour_id, double & x_min, double & x_max, double & y_min, double & y_max, list<HalfEdge *>::iterator & it_l, list<HalfEdge *>::iterator & it_r);

	void update_pixel_grid(int contour_id, cv::Mat & pixel_grid, int pixel_grid_x_0, int pixel_grid_y_1, double x_min, double x_max, double y_min, double y_max, const Point2d & offset, list<HalfEdge *>::iterator it_l, list<HalfEdge *>::iterator it_r, bool verbose = false);

public:
	void get_supporting_segments(set<Segment *> & supporting_segments);

	void get_neighboring_faces(set<Face *> & neighboring_faces);

	void classify();

	void add_non_corner(HalfEdge* h, Vertex* v, HalfEdge* h_1, HalfEdge* h_2, double min_edge_len);

	void remove_non_corners(Vertex* v1, Edge* e1, Vertex* v2, HalfEdge* h); // include the edge for v1

	void set_as_bivalent(Vertex* v);

	void set_as_trivalent(Vertex* v);

	void estimate_label(double alpha, int area, double diag_length); // alpha = weight on regularization, from 0 to 1

	void estimate_label(double & alpha, int & area, double & diag_length, bool ignore_boundary_conditions);

	static void histogram(Matrix<double> & I_prob, Matrix<uchar> & I, vector<Point2i> & pixels,
		vector<double> & hist, int & max_label);

	static void estimate_label(vector<double> & hist, int & n_pixels, int & label,
		double & alpha, int & area, double & diag_length, vector<vector<Adj_Edge_Info>*> & adj_edges_info, bool ignore_boundary_conditions);

	static void compute_merged_histogram(const list<vector<double>> & hists, vector<double> & merged_hist);

	static void merged_halfedges(Face* f_1, Face* f_2, vector<HalfEdges_Location>  & l_12, vector<HalfEdges_Location>  & l_21, vector<list<HalfEdge *>> & merged);

	static int init_start_halfedge(vector<HalfEdge *> & edges, list<pair<int, int>> & h, vector<list<HalfEdge *>> & merged_contours);

	static void merged_contour(int start_pos, Face* f_1, Face* f_2, vector<HalfEdge *> & edges_1, vector<HalfEdge *> & edges_2, HalfEdges_Location & hl_12, HalfEdges_Location & hl_21, vector<list<HalfEdge *>> & merged_contours);

public:
	int id_face;
	int label;
	vector<list<HalfEdge *>> edges; //edges[0] is the clockwise outer contour. For i>=1, edges[i] is a counter-clockwise inner contour
	vector<list<pair<Vertex *, Vertex_Type> >> vertices;
	bool check_split;

	double split_gain;
	vector<double> semantic_probabilities;
	vector<Point2i> pixels;
	vector<set<Point2i, pointcomp>> active_gradients;  // histogram of active gradient pixels
	vector<set<Point2i, pointcomp>> active_lab_gradients; // histogram of active I_lab gradient pixels
	set<Segment *, costcomp> detected_segments;

	vector<vector<Vertex *>> split_vertex_chains;
	map<pair<Vertex *, Vertex *>, Split_Edge_Property> split_edges_property;
	map<Vertex*, pair<int, int>> split_vertex_outer_edge_location; // for undocked split vertex on an outer edge
};

