#pragma once
#include <opencv2/core.hpp>
#include "parameters.h"
#include "kinetic_model.h"
#include "segment_ray.h"
#include "indexed_event.h"
#include "partition_elements.h"
#include "quadtree.h"

using cv::Mat;
using cv::Size2i;

using std::list;
using std::pair;
using std::vector;

class Partition;

namespace Propagation
{
    void propagate_image_domain(Kinetic_Model* model);

	void propagate_inside_polygon(vector<list<HalfEdge *>>& boundary_edges, Partition* graph, Segment_Regularization_Tree* tree, map<int, double> & angles, vector<int> & bin_to_angle_index, vector<Ray *> & rays, vector<Segment *> & segments, int &n_init_vertices, int k_hit = 2);

#if NOT_MEASURING_PERFORMANCES
	void print_histogram(Kinetic_Model *model);
#endif

    void build_rtree(vector<Segment *> & segments, double & D, Boost_RTree & rtree_segments);

    void search_neighborhood(Segment* s, double & D, Boost_RTree & rtree_segments, vector<Segment *> & segments, list<Segment *> & neighborhood);

	void schedule_events_at_polygon_boundaries(vector<Ray *> & rays, vector<list<Outer_Edge *> > & outer_edges, vector<double> & maximal_sizes, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, double min_edge_len, double corner_eps);

	void schedule_events_between_colinear_segments(Segment_Regularization_Tree* tree, vector<IndexedEvent *> & events_colinear);

    void schedule_events_between_non_colinear_segments(vector<Segment *> & segments, vector<Ray *> & rays, vector<double> & maximal_sizes, IndexedEvent** intersectants, IndexedEvent** intersected, vector<IndexedEvent *> & T, double &t_0, double &t_1);

    void schedule_events_between_non_colinear_segments_with_rtree(vector<Segment *> & segments, vector<Ray *> & rays, vector<double> & maximal_sizes, IndexedEvent** intersectants, IndexedEvent** intersected, vector<IndexedEvent *> & T, double &t_0, double &t_1);

	void schedule_events(vector<Segment *> & segments, vector<Ray *> & rays, IndexedEvent** schedule, IndexedEvent** intersectants, IndexedEvent** intersected, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, vector<double> & maximal_sizes, vector<IndexedEvent *> & events_colinear, double t_0, double t_1);

	void insert_indexed_event(IndexedEvent* event, IndexedEvent** intersectants, IndexedEvent** intersected);

	void insert_indexed_event(IndexedEvent* event, IndexedEvent** intersectants, IndexedEvent** intersected, IndexedEvent** intersectants_tails, IndexedEvent** intersected_tails);

	void remove_indexed_event(IndexedEvent** schedule, IndexedEvent* indexed_event, IndexedEvent** intersectants, IndexedEvent** intersected, bool destroy);

	void remove_indexed_event(IndexedEvent** schedule, IndexedEvent* event_at_boundary, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, double t_0, double t_1, bool destroy);

	void remove_references_to_intersectant_ray_inside_polygon(IndexedEvent* upcoming_event, int r_index, IndexedEvent** schedule, IndexedEvent** intersectants, IndexedEvent** intersected, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, double t_0, double t_1);

	void prune_events_and_build_edges(vector<Ray *> & rays, Partition* & graph, vector<list<HalfEdge *>> & boundary_edges, IndexedEvent** first_schedule, IndexedEvent** simplified_schedule, IndexedEvent** last_event,
		IndexedEvent** intersectants, IndexedEvent** intersected, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, Quadtree* & quadtree, double t_0, double t_1, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
		vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges, vector<Vertex *> & initial_vertices_of_rays, int & active_rays);

	IndexedEvent* pop_event(vector<Ray *> & rays, Partition* & graph, vector<list<Outer_Edge *> > & outer_edges, vector<Vertex *> & initial_vertices_of_rays, IndexedEvent** schedule,
		Quadtree* & quadtree, Vertex* & intersection_point, bool & is_corner, bool & is_new_vertex);

	void update_schedule(vector<Ray *> & rays, bool & r_i_continue, IndexedEvent** schedule, IndexedEvent** intersectants, IndexedEvent** intersected, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, double t_0, double t_1, int & active_rays);

	void evaluate_primary_condition(Ray* r_i, Ray* r_j, double t_i, Parameters* params);

	void define_region(Ray* r_i, Point2d & P, int width, double length, Matrix<double> & I, vector<vector<pair<int, int> > > & regions);

#if NOT_MEASURING_PERFORMANCES
	void evaluate_secondary_condition(Ray* r_i, Matrix<double> & I_m, Matrix<double> & I_t, Point2d & P, Parameters* params);

	bool evaluate_secondary_condition_at_region_scale(Ray* r_i, vector<pair<int, int> > & region, Matrix<double> & M, Matrix<double> & T, Parameters* params);
	
	bool evaluate_secondary_condition_double_detector(Ray* r_i, vector<pair<int, int> > & region, Matrix<double> & M, Matrix<double> & T, Parameters* params);

	bool evaluate_secondary_condition_dot_product(Ray* r_i, vector<pair<int, int> > & region, Matrix<double> & M, Matrix<double> & T, Parameters* params);
#endif
}
