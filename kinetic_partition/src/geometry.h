#pragma once
#include <utility>
#include <vector>
#include "segment_ray.h"
#include "segment_tree.h"
#include "comparators.h"

using std::vector;
using std::pair;
using std::make_pair;

class Edge;

class HalfEdge;

class Outer_Edge;

struct rect
{
	double x1, y1, x2, y2;  /* first and second point of the line segment */
	double width, w_max, w_min;        /* rectangle width */
	double x, y;          /* center of the rectangle */
	double theta;        /* angle */
	double dx, dy;        /* (dx,dy) is vector oriented as the line segment */
	double prec;         /* tolerance angle */
	double p;            /* probability of a point with angle within 'prec' */
};

namespace Geometry
{
	int dcmp(double x);

    double mes(Vec2d & _a, Vec2d & _b);

    double distance_initial_coordinates(Segment *s, Segment *t);

	double distance(Point2d & p, Segment *s);

    //bool are_neighbors_euclidian(Segment *s, Segment *t, double d_lim, double & d);

	//bool bounding_boxes_intersect(Segment *s, Segment *t, double d_lim);

	//bool bounding_boxes_intersect(Segment *s, Segment *t);

	void build_tree_polygon_boundaries(vector<Segment *> & segments, vector<list<HalfEdge *>> & polygon_edges, Segment_Regularization_Tree * tree, 
		map<int, double> & angles, vector<int> & bin_to_angle_index, map<HalfEdge *, list<Segment *>> & boundary_collinear_segments);

	void build_tree_when_regularization_disabled(vector<Segment *> & segments, Segment_Regularization_Tree * tree);

    //void find_groups_of_parallel_segments(Segment_Regularization_Tree* tree, vector<vector<int> > & groups);

	//void find_groups_of_parallel_segments(vector<Segment *> & segments, vector<vector<int> > & groups);

	//void find_groups_of_colinear_segments(Segment_Reguarization_Tree* tree, vector<Segment *> & segments, list<pair<int, int> > & colinear_indices);

	//void find_groups_of_colinear_segments(vector<Segment *> & segments, list<pair<int, int> > & colinear_indices);

	void merge_overlapped_segments(Segment_Regularization_Tree *tree);

	bool OnSegment(Point2d P1,Point2d P2,Point2d Q);

	bool is_inside_contour(Point2d & pt, list<HalfEdge *> & contour_edges);

	bool is_inside_contour(Point2d & pt, list<Point2d*> & contour_pts);

	bool is_between(double alpha, double alpha_1, double alpha_2);

	void disable_segments_outside_boundaries(vector<Segment *> & segments, int rows, int cols);

	void disable_segments_outside_polygon_boundaries(vector<Segment *> & segments, vector<list<HalfEdge *>> & polygon_edges);

	void disable_segments_overlapping_polygon_boundaries(vector<Segment *> & segments, vector<list<HalfEdge *>> & polygon_edges, double eps);

	//void merge_overlapped_segments(vector<Segment *> & segments, list<pair<int, int> > & colinear_indices);

	void intersection_polygon_boundary(Ray* r_i, vector<list<Outer_Edge *> > & boundary_edges, vector<double> & t_i, vector<Image_Boundary> & r_j, vector<int> & outer_edges_index, vector<double> & t_j, double min_edge_len, double corner_eps);

	void intersection_boundary(Ray* r_i, int rows, int cols, double & t_i, Image_Boundary & r_j, double & t_j);

	void direct_intersection(Ray* r_i, Image_Boundary b_j, double rows, double cols, double & t_i, double & t_j);

	void direct_intersection(Ray* r_i, Ray* r_j, double & t_i, double & t_j);

	bool intersection(Ray *s, Ray *t, double max_time_s0, double max_time_s1, unsigned int & s_index, unsigned int & t_index, double & s_time, double & t_time);

	bool intersection_colinear(Ray* s, Ray* t, double & s_time, double & t_time);

	bool are_overlapping(Segment *s, Segment *t, Point2d & P, Point2d & Q);

	void region_grow(Point2i seed, Matrix<double> & I_m, Matrix<double> & I_t, vector<Point2i> & reg, double & grad_angle,
		vector<set<Point2i, pointcomp>> & bins, double value_thresh, double bin_width, double prec);

	void region2rect(vector<Point2i> & reg, Matrix<double> & I_m, double reg_angle, double prec, double p, rect * rec);

	double get_rec_theta(vector<Point2i> & reg, double x, double y, Matrix<double> & I_m, double reg_angle, double prec);

	bool refine(vector<Point2i> & reg, int * reg_size, Matrix<double> & I_m, vector<set<Point2i, pointcomp>> & bins, double value_thresh, double bin_width,
		double & reg_angle, double & grad_angle, double prec, double p, rect * rec, Matrix<double> & I_t, double density_th);

	bool reduce_rect_width(vector<Point2i> & reg, int * reg_size, Matrix<double> & I_m, vector<set<Point2i, pointcomp>> & bins, double value_thresh, double bin_width, 
		double & reg_angle, double & grad_angle, double prec, double p, rect * rec, Matrix<double> & I_t, double density_th);

	bool reduce_region_radius(vector<Point2i> & reg, int * reg_size, Matrix<double> & I_m, vector<set<Point2i, pointcomp>> & bins, double value_thresh, double bin_width,
		double & reg_angle, double & grad_angle, double prec, double p, rect * rec, Matrix<double> & I_t, double density_th);

    double angle_diff_signed(double & a, double & b);

	double angle_diff_abs(double & a, double & b);

	double dist(double x1, double y1, double x2, double y2);

	void rect_copy(struct rect * in, struct rect * out);
};

bool inline in(double a, double x, double b) { return (a <= x && x <= b); }

inline bool sort_by_proximity_to_segment(pair<HalfEdge *, pair<double, double>> & h_i, pair<HalfEdge *, pair<double, double>> & h_j) {
	return (h_i.second.first < h_j.second.first);
}