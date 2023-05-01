#include "propagation.h"
#include "geometry.h"
#include <ctime>
#include <queue>
#include "defs.h"
#include "trace.h"
#include "means.h"
#include <iostream>
#include <fstream>
#include "svg.h"
#include "partition.h"

void Propagation::propagate_inside_polygon(vector<list<HalfEdge *>>& boundary_edges, Partition* graph, Segment_Regularization_Tree* tree, map<int, double> & angles, vector<int> & bin_to_angle_index, vector<Ray *> & rays, vector<Segment *> & segments, int & n_init_vertices, int k_hit)
{
	// Propagate rays from segments inside a polygon with holes. Rays stop when meeting the outer contour of the polygon, 
	// but continue propagation when meeting the inner contours. This eliminates ambiguity in resulting facets. 

	// Deletes data from previous execution
	Ray::clear_rays(rays);

	// set parameter
	double min_edge_length = graph->params->prop_min_edge;
	double corner_eps = graph->params->prop_corner_eps;
	int K = k_hit; // intersecting times

	map<HalfEdge *, list<Segment *>> boundary_collinear_segments;
	Geometry::build_tree_polygon_boundaries(segments, boundary_edges, tree, angles, bin_to_angle_index, boundary_collinear_segments);

	Geometry::merge_overlapped_segments(tree);
	Geometry::disable_segments_overlapping_polygon_boundaries(segments, boundary_edges, 1.0);
	Geometry::disable_segments_outside_polygon_boundaries(segments, boundary_edges);

	Ray::build_rays(segments, rays, K, graph->params->split_params.length_thresh);

	size_t nb_rays = rays.size();

	IndexedEvent* first_schedule = NULL;
	IndexedEvent* simplified_schedule = NULL;
	IndexedEvent* last_event_of_simplified_schedule = NULL;

	IndexedEvent** intersectants = new IndexedEvent*[nb_rays];
	for (size_t i = 0; i < nb_rays; i++) {
		intersectants[i] = NULL;
	}
	IndexedEvent** intersected = new IndexedEvent*[nb_rays];
	for (size_t i = 0; i < nb_rays; i++) {
		intersected[i] = NULL;
	}
	vector<vector<IndexedEvent*>> events_at_polygon_boundaries(nb_rays, vector<IndexedEvent*>());

	list<Vertex *> outer_vertices, inner_vertices;
	vector<list<Outer_Edge *> > outer_edges;
	vector<list<Inner_Edge *> > inner_edges;
	vector<Vertex *> initial_vertices_of_rays;
	vector<Vertex *> initial_outer_vertices;

	// Defines the first vertices and edges of the graph :
	// Vertices corresponding to the corners of the polygon, N vertices for the centers of the N segments
	// We also add edges that correspond to the boundaries of the polygon
	graph->init_edges_polygon_boundary(boundary_edges, rays, boundary_collinear_segments, outer_vertices, inner_vertices, outer_edges, inner_edges, initial_vertices_of_rays, initial_outer_vertices);
	n_init_vertices = outer_vertices.size();

	// The computation of events is performed by a loop : while there remain active rays, we compute all the possible
	// intersections happening within a certain range of time [t_0, t_1] and update the graph of intersections.
	vector<double> maximal_sizes;
	vector<IndexedEvent *> events_colinear;
	schedule_events_at_polygon_boundaries(rays, outer_edges, maximal_sizes, events_at_polygon_boundaries, min_edge_length, corner_eps);
	schedule_events_between_colinear_segments(tree, events_colinear);
	int active_rays = nb_rays;
	int iteration = 0;
	double range = graph->params->prop_range;
	double t_0 = -FLT_MAX, t_1 = range;

	while (active_rays > 0) {

		// Computes all events happening between t_0 and t_1
		schedule_events(segments, rays, &first_schedule, intersectants, intersected, events_at_polygon_boundaries, maximal_sizes, events_colinear, t_0, t_1);

		// Simplifies this list of events
		prune_events_and_build_edges(rays, graph, boundary_edges, &first_schedule, &simplified_schedule, &last_event_of_simplified_schedule,
			intersectants, intersected, events_at_polygon_boundaries, graph->quadtree, t_0, t_1,
			outer_vertices, inner_vertices, outer_edges, inner_edges, initial_vertices_of_rays, active_rays);

		// Updates the temporal range
		++iteration;
		t_0 = iteration * range;
		if (active_rays > 250) {
			t_1 = (iteration + 1) * range;
		}
		else {
			// If there are only a few rays left (less than 5%) we can afford to accept any kind of events
			t_1 = FLT_MAX;
		}
	}

	delete[] intersectants;
	delete[] intersected;
	for (int i = 0; i < nb_rays; ++i) {
		for (int j = 0; j < events_at_polygon_boundaries[i].size(); ++j) delete events_at_polygon_boundaries[i][j];
		events_at_polygon_boundaries[i].clear();
	}
	events_at_polygon_boundaries.clear();
	initial_vertices_of_rays.clear();
	

	graph->merge_containers(outer_vertices, inner_vertices, outer_edges, inner_edges);
	
	graph->build_faces();

	// validity check
	bool valid_f = true, valid_e = true, valid_v = true;
	for (list<Face *>::iterator it_f = graph->faces.begin(); it_f != graph->faces.end(); ++it_f) {
		if (HalfEdge::is_clockwise((*it_f)->edges[0]) == 0) {
			valid_f = false;
			std::cout << "invalid facet : " << (*it_f)->id_face << std::endl;
		}
	}

	Edge* e = graph->edges_head;
	while (e != NULL) {
		if (e->v1_v2->f != nullptr || e->v2_v1->f != nullptr) {
			if (e->v1_v2->f == e->v2_v1->f || e->v1 == e->v2 || e->length == 0) {
				valid_e = false;
				std::cout << "invalid edge: " << e->v1->id_vertex << "-->" << e->v2->id_vertex << std::endl;
				if (e->v1_v2->f != nullptr) {
					std::cout << "Facet : " << e->v1_v2->f->id_face << std::endl;
				}
			}
		}
		e = e->e_next;
	}

	Vertex* v = graph->vertices_head;
	while (v != NULL) {
		if (v->connectivity() <= 1) {
			valid_v = false;
		}
		set<Vertex *> connected_vertices;
		for (int i = 0; i < v->directions.size(); ++i) {
			Edge* e_i = v->directions[i].second->e;
			Vertex* connect_v = e_i->v1 == v ? e_i->v2 : e_i->v1;
			if (connected_vertices.count(connect_v) != 0) {
				valid_v = false;
				std::cout << "repeated vertex connection : " << v->id_vertex << "-->" << connect_v->id_vertex << std::endl;
			}
			connected_vertices.insert(connect_v);
		}
		v = v->v_next;
	}

	assert(valid_f);
	assert(valid_e);
	assert(valid_v);

	graph->seal_and_remove_bivalent_vertices(n_init_vertices);

	// Enables again segments whose state has been modified while calling the merge_overlapped segments
	Segment::enable(segments);
}


void Propagation::propagate_image_domain(Kinetic_Model *model)
{

	// Data to be accessed and modified
    vector<Ray *> & rays = model->rays;
    IndexedEvent* & schedule = model->schedule;
    Partition* & graph = model->graph;
	if (graph != NULL) delete graph;
    graph = new Partition(model->I.rows, model->I.cols, model->params);
	
	int rows = graph->rows;
	int cols = graph->cols;

	// To init the graph, we first build a set of temporary vertices that represent the corners of the image and edges that represent the image boundary.
	int temp_id_corners = -4;
	Vertex* top_left = new Vertex(temp_id_corners, -0.5, rows - 0.5);
	Vertex* top_right = new Vertex(temp_id_corners, cols - 0.5, rows - 0.5);
	Vertex* bottom_right = new Vertex(temp_id_corners, cols - 0.5, -0.5);
	Vertex* bottom_left = new Vertex(temp_id_corners, -0.5, -0.5);

	// temporary edges
	int temp_id_boundaries = -4;
	double width = 1.5;
	Outer_Edge* top = new Outer_Edge(temp_id_boundaries, TOP_IMAGE, top_left, top_right, width, 0, PI, true, 0, cols);
	Outer_Edge* bottom = new Outer_Edge(temp_id_boundaries, BOTTOM_IMAGE, bottom_left, bottom_right, width, 0, PI, false, 0, cols);
	Outer_Edge* left = new Outer_Edge(temp_id_boundaries, LEFT_IMAGE, bottom_left, top_left, width, PI / 2, -PI / 2, true, 0, rows);
	Outer_Edge* right = new Outer_Edge(temp_id_boundaries, RIGHT_IMAGE, bottom_right, top_right, width, PI / 2, -PI / 2, false, 0, rows);

	list<HalfEdge*> outer_contour;
	outer_contour.push_back(top->v1_v2);
	outer_contour.push_back(right->v2_v1);
	outer_contour.push_back(bottom->v2_v1);
	outer_contour.push_back(left->v1_v2);
	vector<list<HalfEdge*>> boundary_edges = {outer_contour};

	int n_init_vertices = 4;
	propagate_inside_polygon(boundary_edges, model->graph, model->tree, model->reg_angles, model->bin_to_angle_index, model->rays, model->segments, n_init_vertices, model->params->prop_K);
	graph->reset_indices();
	// free memory
	delete top;
	delete bottom;
	delete left;
	delete right;
	delete top_left;
	delete top_right;
	delete bottom_right;
	delete bottom_left;
}


void Propagation::build_rtree(vector<Segment *> & segments, double & D, Boost_RTree & rtree_segments)
{
    double x_min, x_max, y_min, y_max;
    for (uint i = 0 ; i < segments.size() ; i++) {

        // For each segment, we create a box that represents the bounding box of this segment
        Segment* s_i = segments[i];
        if (!s_i->is_disabled) {
            if (s_i->finalEnd1.x < s_i->finalEnd2.x) {
                x_min = s_i->finalEnd1.x; x_max = s_i->finalEnd2.x;
            } else {
                x_min = s_i->finalEnd2.x; x_max = s_i->finalEnd1.x;
            }
            if (s_i->finalEnd1.y < s_i->finalEnd2.y) {
                y_min = s_i->finalEnd1.y; y_max = s_i->finalEnd2.y;
            } else {
                y_min = s_i->finalEnd2.y; y_max = s_i->finalEnd1.y;
            }

            // We insert this box in the r-tree
            Boost_Box s_box(Boost_Point(x_min - D, y_min - D), Boost_Point(x_max + D, y_max + D));
            rtree_segments.insert(std::make_pair(s_box, i));
        }
    }
}



void Propagation::search_neighborhood(Segment* s, double & D, Boost_RTree & rtree_segments, vector<Segment *> & segments, list<Segment *> & neighborhood)
{
    neighborhood.clear();

    // We search for segments whose bounding boxes' edges are located at a distance lower than D from the bounding box of s
    double x_min, x_max, y_min, y_max;
    vector<Boost_Value> possible_neighbors;

    if (s->finalEnd1.x < s->finalEnd2.x) {
        x_min = s->finalEnd1.x; x_max = s->finalEnd2.x;
    } else {
        x_min = s->finalEnd2.x; x_max = s->finalEnd1.x;
    }
    if (s->finalEnd1.y < s->finalEnd2.y) {
        y_min = s->finalEnd1.y; y_max = s->finalEnd2.y;
    } else {
        y_min = s->finalEnd2.y; y_max = s->finalEnd1.y;
    }
    Boost_Box query(Boost_Point(x_min - D, y_min - D), Boost_Point(x_max + D, y_max + D));
    rtree_segments.query(bgi::intersects(query), std::back_inserter(possible_neighbors));

    for (uint i = 0 ; i < possible_neighbors.size() ; i++) {
        Segment* s_i = segments[possible_neighbors[i].second];
        neighborhood.push_back(s_i);
    }
}


void Propagation::schedule_events_at_polygon_boundaries(vector<Ray *> & rays, vector<list<Outer_Edge *> > & outer_edges, vector<double> & maximal_sizes, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, double min_edge_len, double corner_eps)
{
	uint nb_rays = uint(rays.size());
	maximal_sizes = vector<double>(nb_rays, 0.0);
	vector<Image_Boundary> boundaries;
	vector<int> outer_edges_index;
	vector<double> t_i = { FLT_MAX }, t_j = { FLT_MAX };

	for (uint i = 0; i < nb_rays; i++) {
		Geometry::intersection_polygon_boundary(rays[i], outer_edges, t_i, boundaries, outer_edges_index, t_j, min_edge_len, corner_eps);
		maximal_sizes[i] = t_i[0]; // as the time that intersects the polygon outer contour
		uint nb_intersected = uint(t_i.size());
		events_at_polygon_boundaries[i].resize(nb_intersected);
		for (uint j = 0; j < nb_intersected; j++) {
			assert(boundaries[j] < 0);
			events_at_polygon_boundaries[i][j] = new IndexedEvent(i, boundaries[j], outer_edges_index[j], t_i[j], t_j[j]);
		}
	}
}


void Propagation::schedule_events_between_colinear_segments(Segment_Regularization_Tree* tree, vector<IndexedEvent *> & events_colinear)
{

	// Accesses to colinear segments
	for (map<double, Node_Parallel_Segments*>::iterator it_m1 = tree->parallel_segments.begin() ; it_m1 != tree->parallel_segments.end() ; it_m1++) {
		Node_Parallel_Segments* node_parallel = it_m1->second;
		for (map<double, Node_Colinear_Segments*>::iterator it_m2 = node_parallel->colinear_segments.begin() ; it_m2 != node_parallel->colinear_segments.end() ; it_m2++) {
			Node_Colinear_Segments* node_colinear = it_m2->second;

			if (node_colinear->colinear_segments.size() == 1) {
				continue;
			} else {
				list<Segment *>::iterator it_s1, it_s2;
				for (it_s1 = node_colinear->colinear_segments.begin() ; it_s1 != node_colinear->colinear_segments.end() ; it_s1++) {
					Segment* s1 = (*it_s1);
					if (s1->is_disabled) {
						continue;
					}
					it_s2 = it_s1;
					while (++it_s2 != node_colinear->colinear_segments.end()) {
						Segment* s2 = (*it_s2);
						if (s2->is_disabled) {
							continue;
						}

						bool loop = true;
						for (int i = 0; i < 2 && loop; i++) {
							Ray* r_i = (i == 0 ? s1->rays.first : s1->rays.second);
							for (int j = 0; j < 2 && loop; j++) {
								Ray* r_j = (j == 0 ? s2->rays.first : s2->rays.second);
								double t_i, t_j;
								bool exists = Geometry::intersection_colinear(r_i, r_j, t_i, t_j);
								if (exists) {
									// We choose to duplicate the event
									// Later, when events will be pruned, every time we will consider this kind of event we will
									// also has to take into account the "twin" event
									IndexedEvent* event_ij = new IndexedEvent(int(r_i->index), int(r_j->index), t_i, t_j, true);
									IndexedEvent* event_ji = new IndexedEvent(int(r_j->index), int(r_i->index), t_j, t_i, true);
									// Insert both events in the table
									events_colinear.push_back(event_ij);
									events_colinear.push_back(event_ji);
									loop = false;
								}
							}
						}
					}
				}
			}
		}
	}

	std::sort(events_colinear.begin(), events_colinear.end(), before_indexed_event);
}


void Propagation::schedule_events_between_non_colinear_segments(vector<Segment *> & segments, vector<Ray *> & rays, vector<double> & maximal_sizes, IndexedEvent** intersectants, IndexedEvent** intersected, vector<IndexedEvent *> & T, double & t_0, double & t_1)
{
    // We use a structure of sparse matrix to store indexed events : but we also know in which order the events are created
    // Indeed, each new element is added at the end of the current row, or at the end of the current column
    // That's why we use the following tables : to get a quick reference on the last element of each row and the last element of each column
    uint nb_rays = uint(rays.size());
    IndexedEvent** intersectant_tails = new IndexedEvent*[nb_rays];
    IndexedEvent** intersected_tails = new IndexedEvent*[nb_rays];
    for (unsigned int i = 0; i < nb_rays; i++) {
        intersectant_tails[i] = NULL;
        intersected_tails[i] = NULL;
    }

    // The algorithm to compute the intersection of two rays r_i and r_j takes the shape of a double loop on i and j,
    // where i and j represent the indices of the segments. We suppose that i is the intersectant ray.
    uint nb_segments = uint(segments.size());
    for (uint i = 0 ; i < nb_segments ; i++) {
        Segment* s_i = segments[i];

        // If the i-th segment has been discarded, or if both rays related to it have already been stopped, jump to the next segment
        if (s_i->is_disabled || s_i->rays.first == NULL || (s_i->rays.first->has_been_stopped() && s_i->rays.second->has_been_stopped())) continue;
        Ray* r_i = s_i->rays.first;

        // For every segment likely to be intersected by the segment i-th
        for (uint j = 0 ; j < nb_segments ; j++) {
            if (i == j) continue;
			Segment* s_j = segments[j];

            // If segment j is disabled, we jump to the next segment
            if (s_j->is_disabled || s_j->rays.first == NULL) continue;
            Ray* r_j = s_j->rays.first;

            uint ind_i = -1, ind_j = -1;
            double t_i = -FLT_MAX, t_j = -FLT_MAX;
            bool exists = Geometry::intersection(r_i, r_j, maximal_sizes[r_i->index], maximal_sizes[r_i->index + 1], ind_i, ind_j, t_i, t_j);
            if (exists) {
                uint index_intersectant, index_intersected;
                double t_intersectant, t_intersected;
                if (t_i >= t_j) {
                    index_intersectant = ind_i;
                    index_intersected = ind_j;
                    t_intersectant = t_i;
                    t_intersected = t_j;
                } else {
                    index_intersectant = ind_j;
                    index_intersected = ind_i;
                    t_intersectant = t_j;
                    t_intersected = t_i;
                }
                // If the intersectant ray actually belongs to the i-th segment and if it is still propagating
                if (rays[index_intersectant]->parent->index == s_i->index && !rays[index_intersectant]->has_been_stopped() && t_0 <= t_intersectant && t_intersectant < t_1) {
                    // There are two possibilites for an a priori valid intersection :
                    // - the intersected ray is still active
                    // - the intersected ray has been stopped but it is intersected in a point that actually exists (check time)
                    if (!rays[index_intersected]->has_been_stopped() || (rays[index_intersected]->has_been_stopped() && t_intersected <= rays[index_intersected]->t)) {
                        IndexedEvent* event = new IndexedEvent(index_intersectant, index_intersected, t_intersectant, t_intersected, false);
                        // Insert it in the table
                        insert_indexed_event(event, intersectants, intersected, intersectant_tails, intersected_tails);
                        // Inserts it in a list, which will be later sorted
                        T.push_back(event);
                    }
                }
            }
        }
    }
    delete[] intersectant_tails;
    delete[] intersected_tails;
}


void Propagation::schedule_events_between_non_colinear_segments_with_rtree(vector<Segment *> & segments, vector<Ray *> & rays, vector<double> & maximal_sizes, IndexedEvent** intersectants, IndexedEvent** intersected, vector<IndexedEvent *> & T, double &t_0, double & t_1)
{
    Boost_RTree rtree_segments;
    build_rtree(segments, t_1, rtree_segments);

    // Initializes tails
    uint nb_rays = uint(rays.size());
    IndexedEvent** intersectant_tails = new IndexedEvent*[nb_rays];
    IndexedEvent** intersected_tails = new IndexedEvent*[nb_rays];
    for (unsigned int i = 0; i < nb_rays; i++) {
        intersectant_tails[i] = NULL;
        intersected_tails[i] = NULL;
    }

    uint nb_segments = uint(segments.size());
    for (uint i = 0 ; i < nb_segments ; i++) {
        Segment* s_i = segments[i];

        if (s_i->is_disabled || s_i->rays.first == NULL || (s_i->rays.first->has_been_stopped() && s_i->rays.second->has_been_stopped())) continue;
        Ray* r_i = s_i->rays.first;

        // Searches for all the segments which are likely to be intersected by s_i
        list<Segment *> neighbors;
        search_neighborhood(s_i, t_1, rtree_segments, segments, neighbors);
        vector<IndexedEvent*> local_events;

        // For every segment likely to be intersected by the i-th segment
        for (list<Segment *>::iterator it_s = neighbors.begin() ; it_s != neighbors.end() ; it_s++) {

            Segment* s_j = (*it_s);
            if (s_i == s_j || s_j->rays.first == NULL) continue;
            Ray* r_j = s_j->rays.first;

            uint ind_i = -1, ind_j = -1;
            double t_i = -FLT_MAX, t_j = -FLT_MAX;
            bool exists = Geometry::intersection(r_i, r_j, maximal_sizes[r_i->index], maximal_sizes[r_i->index + 1], ind_i, ind_j, t_i, t_j);
            if (exists) {
                uint index_intersectant, index_intersected;
                double t_intersectant, t_intersected;
                if (t_i >= t_j) {
                    index_intersectant = ind_i;
                    index_intersected = ind_j;
                    t_intersectant = t_i;
                    t_intersected = t_j;
                } else {
                    index_intersectant = ind_j;
                    index_intersected = ind_i;
                    t_intersectant = t_j;
                    t_intersected = t_i;
                }

                if (rays[index_intersectant]->parent->index == s_i->index && !rays[index_intersectant]->has_been_stopped() && t_0 <= t_intersectant && t_intersectant < t_1) {
                    if (!rays[index_intersected]->has_been_stopped() || (rays[index_intersected]->has_been_stopped() && t_intersected <= rays[index_intersected]->t)) {
                        IndexedEvent* event = new IndexedEvent(index_intersectant, index_intersected, t_intersectant, t_intersected, false);
                        // We do not insert this event in the matrix of events for the moment
                        // Instead we put it in a vector, sorted on index_intersected
                        local_events.push_back(event);
                    }
                }
            }
        }

        // Sorts the table of events, and inserts elements in the chronological list of events, as well as in the matrix of events
        std::sort(local_events.begin(), local_events.end(), sorted_by_index_of_intersected_segment);
        for (uint l = 0 ; l < local_events.size() ; l++) {
            insert_indexed_event(local_events[l], intersectants, intersected, intersectant_tails, intersected_tails);
            T.push_back(local_events[l]);
        }
    }

    rtree_segments.clear();
    delete[] intersectant_tails;
    delete[] intersected_tails;
}


void Propagation::schedule_events(vector<Segment *> & segments, vector<Ray *> & rays, IndexedEvent** schedule, IndexedEvent** intersectants, IndexedEvent** intersected, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, vector<double> & maximal_sizes, vector<IndexedEvent *> & events_colinear, double t_0, double t_1)
{
	vector<IndexedEvent *> T = vector<IndexedEvent *>();
	uint nb_rays = uint(rays.size());

	// First step : selects the possible intersections of the rays with the boundaries
	for (uint i = 0; i < nb_rays; i++) {
		vector<IndexedEvent*> events = events_at_polygon_boundaries[i];
		uint nb_intersected = uint(events.size());
		for (uint j = 0; j < nb_intersected; ++j) {
			if (!rays[i]->has_been_stopped() && events[j]->t_intersectant >= t_0 && events[j]->t_intersectant < t_1) {
				T.push_back(events[j]);
			}
		}
	}

	// Second, computes the intersections between all the non-colinear rays
	// Our algorithm takes the shape of a double loop on segments
	// To reduce the complexity we may use a structure of R-Tree
	if (t_1 < FLT_MAX) {
		schedule_events_between_non_colinear_segments_with_rtree(segments, rays, maximal_sizes, intersectants, intersected, T, t_0, t_1);
	}
	else {
		schedule_events_between_non_colinear_segments(segments, rays, maximal_sizes, intersectants, intersected, T, t_0, t_1);
	}

	// Third step : take into account all the intersections between colinear rays
	if (events_colinear.size() > 0) {
		vector<IndexedEvent *>::iterator it = events_colinear.begin();
		IndexedEvent* event = NULL;
		while (it != events_colinear.end() && (*it)->t_intersectant < t_1) {
			event = (*it);
			if (!rays[event->intersectant]->has_been_stopped() && !rays[event->intersected]->has_been_stopped()) {
				insert_indexed_event(event, intersectants, intersected);
				T.push_back(event);
			}
			else {
				delete event;
			}
			it++;
		}
		events_colinear.erase(events_colinear.begin(), it);
	}

	// Final step : obtains a list of chronologically sorted events and deletes the auxiliary vector
	//std::cout << T.size() << " events scheduled" << std::endl;
	std::sort(T.begin(), T.end(), before_indexed_event);

	if (T.size() > 0) {
		*schedule = T[0];
	}
	else {
		*schedule = NULL;
	}
	for (int i = 1; i < int(T.size()); i++) T[i]->previous = T[i - 1];
	for (int i = 0; i < int(T.size()) - 1; i++) T[i]->next = T[i + 1];
	T.clear();
}


void Propagation::insert_indexed_event(IndexedEvent* event, IndexedEvent** intersectants, IndexedEvent** intersected)
{
	int i = event->intersectant;
	int j = event->intersected;
	// Inserts this event at the appropriate location horizontally
	IndexedEvent* previous_event = NULL;
	IndexedEvent* current_event = intersectants[i];
	while (current_event != NULL && current_event->intersected < j) {
		previous_event = current_event;
		current_event = current_event->next_j;
	}
	if (current_event != NULL) {
		event->next_j = current_event;
		current_event->prev_j = event;
	}
	if (previous_event != NULL) {
		event->prev_j = previous_event;
		previous_event->next_j = event;
	} else {
		intersectants[i] = event;
	}
	// Inserts this event at the appropriate location vertically
	previous_event = NULL;
	current_event = intersected[j];
	while (current_event != NULL && current_event->intersectant < i) {
		previous_event = current_event;
		current_event = current_event->next_i;
	}
	if (current_event != NULL) {
		event->next_i = current_event;
		current_event->prev_i = event;
	}
	if (previous_event != NULL) {
		event->prev_i = previous_event;
		previous_event->next_i = event;
	} else {
		intersected[j] = event;
	}
}


void Propagation::insert_indexed_event(IndexedEvent* event, IndexedEvent** intersectants, IndexedEvent** intersected, IndexedEvent** intersectants_tails, IndexedEvent** intersected_tails) {
	int i = event->intersectant;
	int j = event->intersected;
	IndexedEvent* intersectant_tail = intersectants_tails[i];
	IndexedEvent* intersected_tail = intersected_tails[j];
	if (intersectant_tail == NULL) {
		intersectants[i] = event;
	} else {
		intersectant_tail->next_j = event;
		event->prev_j = intersectant_tail;
	}
	intersectants_tails[i] = event;
	if (intersected_tail == NULL) {
		intersected[j] = event;
	} else {
		intersected_tail->next_i = event;
		event->prev_i = intersected_tail;
	}
	intersected_tails[j] = event;
}


void Propagation::remove_indexed_event(IndexedEvent** schedule, IndexedEvent* indexed_event, IndexedEvent** intersectants, IndexedEvent** intersected, bool destroy) {
	IndexedEvent* previous = indexed_event->previous;
	IndexedEvent* next = indexed_event->next;
	int r = indexed_event->intersectant;
	int s = indexed_event->intersected;
	// Updates schedule
	if (previous != NULL) {
		previous->next = next;
	} else {
		*schedule = next;
	}
	if (next != NULL) {
		next->previous = previous;
	}
	// Updates sparse matrix
	IndexedEvent* previous_i = indexed_event->prev_i;
	IndexedEvent* previous_j = indexed_event->prev_j;
	IndexedEvent* next_i = indexed_event->next_i;
	IndexedEvent* next_j = indexed_event->next_j;
	if (previous_j == NULL) {
		intersectants[r] = next_j;
	} else {
		previous_j->next_j = next_j;
	}
	if (next_j != NULL) {
		next_j->prev_j = previous_j;
	}
	if (previous_i == NULL) {
		intersected[s] = next_i;
	} else {
		previous_i->next_i = next_i;
	}
	if (next_i != NULL) {
		next_i->prev_i = previous_i;
	}
	// Calls destructor if needed
	if (destroy) delete indexed_event;
}


void Propagation::remove_indexed_event(IndexedEvent** schedule, IndexedEvent* event_at_boundary, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, double t_0, double t_1, bool destroy) {
	IndexedEvent* previous = event_at_boundary->previous;
	IndexedEvent* next = event_at_boundary->next;
	int r = event_at_boundary->intersectant;
	// Updates schedule if the event is included in it
	if (t_0 <= event_at_boundary->t_intersectant && event_at_boundary->t_intersectant < t_1) {
		if (previous != NULL) {
			previous->next = next;
		}
		else {
			*schedule = next;
		}
		if (next != NULL) {
			next->previous = previous;
		}
	}

	// Updates table of events at boundary
	vector<IndexedEvent*>::iterator it = std::find(events_at_polygon_boundaries[r].begin(), events_at_polygon_boundaries[r].end(), event_at_boundary);
	assert(it != events_at_polygon_boundaries[r].end());
	events_at_polygon_boundaries[r].erase(it);
	
	// Calls destructor if needed
	if (destroy) delete event_at_boundary;
}


void Propagation::remove_references_to_intersectant_ray_inside_polygon(IndexedEvent* upcoming_event, int index, IndexedEvent** schedule, IndexedEvent** intersectants, IndexedEvent** intersected, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, double t_0, double t_1)
{
	// Removes all the events where ray 'ray_index' plays the role of intersectant
	IndexedEvent* current_event = intersectants[index];
	IndexedEvent* next_event = NULL;
	while (current_event != NULL) {
		next_event = current_event->next_j;
		bool destroy = (current_event != upcoming_event);
		remove_indexed_event(schedule, current_event, intersectants, intersected, destroy);
		current_event = next_event;
	}

	// Removes some events where ray r plays the role of intersected (those after the upcoming event, or whose intersectant has already been intersected in the upcoming event)
	current_event = intersected[index];
	next_event = NULL;
	while (current_event != NULL) {
		next_event = current_event->next_i;
		if (current_event->t_intersected > upcoming_event->t_intersectant || current_event->intersectant == upcoming_event->intersected) {
			remove_indexed_event(schedule, current_event, intersectants, intersected, current_event != upcoming_event);
		}
		current_event = next_event;
	}

	// Pops events for the ray involving boundaries
	for (vector<IndexedEvent*>::reverse_iterator it_r = events_at_polygon_boundaries[index].rbegin(); it_r != events_at_polygon_boundaries[index].rend(); ++it_r) {
		IndexedEvent* event = *it_r;
		bool destroy = (event != upcoming_event);
		remove_indexed_event(schedule, event, events_at_polygon_boundaries, t_0, t_1, destroy);
	}
}


void Propagation::prune_events_and_build_edges(vector<Ray *> & rays, Partition* & graph, vector<list<HalfEdge *>> & boundary_edges, IndexedEvent** first_schedule, IndexedEvent** simplified_schedule, IndexedEvent** last_event, 
	IndexedEvent** intersectants, IndexedEvent** intersected, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, Quadtree* & quadtree, double t_0, double t_1, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
	vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges, vector<Vertex *> & initial_vertices_of_rays, int & active_rays)
{
	Vertex* intersection_point = NULL;
	bool is_corner = false;
	bool is_new_vertex = true;
	bool r_i_continue = true;

	// Pops the first event of the raw schedule
	// If last_event == NULL, the processed event is the first to be processed
	// If last_event != NULL, some events have already be processed at a previous iteration, which means that the behaviour
	// of the algorithm depends on the value of curr_event : if it is not NULL, we append this event to the simplified schedule
	// with the help of last_event, otherwise we do nothing.

	IndexedEvent* curr_event = pop_event(rays, graph, outer_edges, initial_vertices_of_rays, first_schedule, quadtree, intersection_point, is_corner, is_new_vertex);
	if (*last_event == NULL) {
		*simplified_schedule = curr_event;
		if (*simplified_schedule != NULL) {
			graph->build_edge(rays, curr_event, r_i_continue, intersection_point, is_corner, is_new_vertex, outer_vertices, inner_vertices, outer_edges, inner_edges, initial_vertices_of_rays);
			update_schedule(rays, r_i_continue, first_schedule, intersectants, intersected, events_at_polygon_boundaries, t_0, t_1, active_rays);
			(*simplified_schedule)->previous = NULL;
			*last_event = curr_event;
		}
	}
	else if (curr_event != NULL) {
		(*last_event)->next = curr_event;
		if (curr_event != NULL) {
			graph->build_edge(rays, curr_event, r_i_continue, intersection_point, is_corner, is_new_vertex, outer_vertices, inner_vertices, outer_edges, inner_edges, initial_vertices_of_rays);
			update_schedule(rays, r_i_continue, first_schedule, intersectants, intersected, events_at_polygon_boundaries, t_0, t_1, active_rays);
			curr_event->previous = *last_event;
			*last_event = curr_event;
		}
	}
	
	// While there are events to process, we append them to the simplified schedule using last_event.
	// Naturally, we update this variable at each iteration.
	while (curr_event != NULL) {
		curr_event = pop_event(rays, graph, outer_edges, initial_vertices_of_rays, first_schedule, quadtree, intersection_point, is_corner, is_new_vertex);
		(*last_event)->next = curr_event;
		if (curr_event != NULL) {
			graph->build_edge(rays, curr_event, r_i_continue, intersection_point, is_corner, is_new_vertex, outer_vertices, inner_vertices, outer_edges, inner_edges, initial_vertices_of_rays);
			update_schedule(rays, r_i_continue, first_schedule, intersectants, intersected, events_at_polygon_boundaries, t_0, t_1, active_rays);
			curr_event->previous = *last_event;
			*last_event = curr_event;
		}
	}
}


IndexedEvent* Propagation::pop_event(vector<Ray *> & rays, Partition* & graph, vector<list<Outer_Edge *> > & outer_edges, vector<Vertex *> & initial_vertices_of_rays, IndexedEvent** schedule, 
	Quadtree* & quadtree, Vertex* & intersection_point, bool & is_corner, bool & is_new_vertex)
{
	if (*schedule == NULL) {
		return NULL;
	}
	else {
		Parameters* params = graph->params;

		// Gets the first element of the schedule
		IndexedEvent* upcoming_event = *schedule;
		int i = upcoming_event->intersectant;
		int j = upcoming_event->intersected;
		double t_i = upcoming_event->t_intersectant, t_j = upcoming_event->t_intersected;

		// This upcoming event defines a new vertex in the graph, unless multiple rays intersect in the same point
		// So, let's find or create such a vertex
		Ray *r_i = rays[i], *r_j = (j >= 0 ? rays[j] : nullptr);

		Point2d pt;
		is_corner = Vertex::approximate_coordinates_polygon_corners(upcoming_event, r_i, r_j, t_i, outer_edges, pt, intersection_point, params->prop_min_edge, params->prop_corner_eps);

		bool meets_colinear_ray = false;
		if (is_corner) {
			intersection_point->events.push_back(upcoming_event);
			is_new_vertex = false;
			if (intersection_point->boundary_vertex_on_colinear_ray(r_i, rays)) {
				meets_colinear_ray = true;
			}
		}
		else {
			bool should_create_vertex = true;

			list<Vertex *> neighbors;
			quadtree->search(pt, 1, neighbors);
			if (!neighbors.empty()) {
				if (r_j == nullptr) {
					// If the ray intersects a boundary. 
					// We check: 
					// 1. the distance to the closest neighbor vertex that falls on the intersected boundary.
					// 2. if the intersection point coincides an existing point on r_i.
					// 3. if the intersection point coincides with an existing point on a collinear ray that has already intersected this boundary.
					bool is_vertex_on_ri = false;
					list<Outer_Edge *> & edges_for_this_corner = outer_edges[upcoming_event->intersected_index];
					double dist_pt_v = 1;
					Vertex* v = nullptr;
					for (list<Vertex *>::iterator it_v = neighbors.begin(); it_v != neighbors.end(); ++it_v) {
						double next_dist = norm(pt - (*it_v)->pt);
						if (next_dist < dist_pt_v) {
							bool on_intersected = false;
							for (int i = 0; i < (*it_v)->directions.size(); ++i) {
								Edge* e = (*it_v)->directions[i].second->e;
								if (e->type == OUTER_EDGE) {
									Outer_Edge* oe = static_cast<Outer_Edge *>(e);
									if (oe->init_boundary_location == edges_for_this_corner.front()->init_boundary_location) {
										dist_pt_v = next_dist;
										v = (*it_v);
										on_intersected = true;
										break;
									}
								}
							}
						}
						
						if (next_dist < 1e-6) {
							for (int i = 0; i < (*it_v)->directions.size(); ++i) {
								Edge* e = (*it_v)->directions[i].second->e;
								if (e->type == INNER_EDGE) {
									Inner_Edge * ie = static_cast<Inner_Edge *>(e);
									if (ie->rays.count(r_i) > 0) {
										dist_pt_v = next_dist;
										v = (*it_v);
										is_vertex_on_ri = true; break;
									}
								}
							}
						}
					}

					if ( (dist_pt_v < params->prop_min_edge && v != nullptr) || is_vertex_on_ri) {
						intersection_point = v;
						is_new_vertex = false;
						should_create_vertex = false;
						if (intersection_point->boundary_vertex_on_colinear_ray(r_i, rays)) {
							meets_colinear_ray = true;
						}
					}
				}
				else {
					// The ray intersects another ray.
					// If the intersectant ray is colinear to another ray,
					// we determine if the point where it meets the intersectant ray is issued from the intersection of a colinear ray.

					for (list<Vertex *>::iterator it_v = neighbors.begin(); it_v != neighbors.end(); ++it_v) {
						Vertex* v = (*it_v);
						if (v->inner_vertex_created_by_colinear_ray(r_i, r_j, rays, false)) {
							intersection_point = v;
							is_new_vertex = false;
							should_create_vertex = false;
							meets_colinear_ray = true;
							break;
						}
						else { // we check if the vertex has already been created by the intersection between the intersected/intersectant ray with another ray, 
							// or it coincides with the initial vertex of the intersected/intersectant ray
							double next_dist = cv::norm(v->pt - pt);
							if (next_dist < params->prop_min_edge) {
								bool is_vertex_on_rj = false;
								bool is_vertex_on_ri = false;
								for (int i = 0; i < v->directions.size(); ++i) {
									Edge* e = v->directions[i].second->e;
									if (e->type == INNER_EDGE) {
										Inner_Edge * ie = static_cast<Inner_Edge *>(e);
										if (ie->rays.count(r_j) > 0) {
											is_vertex_on_rj = true; break;
										}
										else if (ie->rays.count(r_i) > 0) {
											is_vertex_on_ri = true; break;
										}
									}
								}
								if (v == initial_vertices_of_rays[j]) is_vertex_on_rj = true;
								if (v == initial_vertices_of_rays[i]) is_vertex_on_ri = true;

								if (is_vertex_on_rj || is_vertex_on_ri) {
									// If two rays intersect sufficiently close to an existing vertex on r_j, we do not create a new vertex.
									intersection_point = v;
									is_new_vertex = false;
									should_create_vertex = false;
								}
							}

						}
					}

				}
			}

			//if (r_i->parent->node_colinear != NULL) {
			//list<Vertex *> neighbor_vertices;
			//quadtree->search(pt, 1, neighbor_vertices);
			//if (!neighbor_vertices.empty()) {
			//intersection_point = NULL;
			//for (list<Vertex *>::iterator it_v = neighbor_vertices.begin(); it_v != neighbor_vertices.end(); it_v++) {
			//if ((*it_v)->created_by_colinear_ray(r_i, r_i->parent->node_colinear, rays)) {
			//intersection_point = (*it_v);
			//is_new_vertex = false;
			//intersection_point->events.push_back(upcoming_event);
			//meets_colinear_ray = true;
			//break;
			//}
			//}
			//}
			//}

			// If none of the conditions above could be fulfilled we create a new vertex
			if (should_create_vertex) {
				intersection_point = new Vertex(graph->get_id_vertices(), upcoming_event, pt);
				quadtree->add(intersection_point);
				is_new_vertex = true;
			}
			else { // Add event to the existing intersection point
				intersection_point->events.push_back(upcoming_event);
				assert(intersection_point->has_no_duplicated_event());
			}
		}
		// The ray r stops propagating if one of the four conditions is filled :
		// - it has interescted the polygon outer boundary (does not stop at inner boundaries)
		// - it has intersected a colinear ray
		// - it has intersected a non-colinear ray in a point where it meets a colinear ray
		// We finally define the stop condition
		bool force_stop = (j == POLYGON_OUTER) || upcoming_event->is_colliding_colinear || meets_colinear_ray;
		if (force_stop) {
			r_i->stop();
		}
		else {
			if (r_i->primary_condition) {
				// keep the primary condition when intersecting a polygon inner boundary
				if (j != POLYGON_INNER) {
					evaluate_primary_condition(r_i, r_j, t_i, params);
				}
			}
			r_i->should_ray_be_stopped();
		}
		return upcoming_event;
	}
}

void Propagation::update_schedule(vector<Ray *> & rays, bool & r_i_continue, IndexedEvent** schedule, IndexedEvent** intersectants, IndexedEvent** intersected, vector<vector<IndexedEvent*>> & events_at_polygon_boundaries, double t_0, double t_1, int & active_rays)
{
	if (*schedule == NULL) return;

	IndexedEvent* current_event = *schedule;
	int i = current_event->intersectant;
	int j = current_event->intersected;
	double t_i = current_event->t_intersectant, t_j = current_event->t_intersected;

	// This upcoming event defines a new vertex in the graph, unless multiple rays intersect in the same point
	// So, let's find or create such a vertex
	Ray *r_i = rays[i], *r_j = (j >= 0 ? rays[j] : nullptr);

	if (r_i->has_been_stopped() && !r_i_continue) {
		// If the ray stops its propagation, we erase :
		// - all future events actively involving rays[i]
		// - all events passively involving rays[i] which are posterior to the current time
		// - the event on the boundary, if it still exists
		remove_references_to_intersectant_ray_inside_polygon(current_event, i, schedule, intersectants, intersected, events_at_polygon_boundaries, t_0, t_1);

		// If, in addition to this, if the upcoming event is a collision between two colinear rays,
		// we set the time to leave of the other ray to 0 and remove all references to that ray
		if (current_event->is_colliding_colinear) {
			r_j->stop();
			remove_references_to_intersectant_ray_inside_polygon(current_event, j, schedule, intersectants, intersected, events_at_polygon_boundaries, t_0, t_1);
			r_j->set_time(t_j);
			active_rays--;
		}
		r_i->set_time(t_i);
		active_rays--;
	}
	else {
		// Otherwise, if the ray doesn't stop its propagation, we simply erase the reference to this event in the matrix of events
		if (j != POLYGON_INNER && j != POLYGON_OUTER) {
			remove_indexed_event(schedule, current_event, intersectants, intersected, false);
		}
		else {
			remove_indexed_event(schedule, current_event, events_at_polygon_boundaries, t_0, t_1, false);
		}
	}

	// Updates the pointers of the structure, and finally returns the event
	if (*schedule != NULL) (*schedule)->previous = NULL;
	current_event->prev_i = NULL;
	current_event->prev_j = NULL;
	current_event->next_i = NULL;
	current_event->next_j = NULL;
}



void Propagation::evaluate_primary_condition(Ray* r_i, Ray* r_j, double t_i, Parameters* params)
{
	if (t_i < 0) {
		r_i->primary_condition = true;
	} else {
		switch (params->prop_policy) {
		case 1:
			// The ray should be stopped if it intersects another ray for the k-th time
			r_i->ttl--;
			r_i->primary_condition = (r_i->ttl > 0);
			break;
		case 2:
			// The ray should be stopped if it intersects another ray for the k-th time
			// or if it intersects another non orthogonal ray
			r_i->ttl--;
			r_i->primary_condition = (r_i->ttl > 0 && r_i->intersects_orthogonal_ray(r_j));
			break;
		case 3:
			// The ray should be stopped if it has lived more than t units of time
			r_i->primary_condition = (t_i <= params->prop_distance);
			break;
		case 4:
			// The ray should be stopped if it has lived more than t units of time
			// or if it intersects another non orthogonal ray
			r_i->primary_condition = (t_i <= params->prop_distance && r_i->intersects_orthogonal_ray(r_j));
			break;
		}
	}
}


void Propagation::define_region(Ray* r_i, Point2d & P, int width, double length, Matrix<double> & I, vector<vector<pair<int, int> > > & region)
{
	typedef pair<int, int> Pixel;
	uint rows = I.rows, cols = I.cols;

	// A region of interest is defined by lines of pixels going from P to Q, Q being located at a distance r_length from P.
	// However if r_width > 1, P and Q move along the normal axis to the ray r_i.

	Vec2d r_dir = r_i->OA;
	Vec2d r_nor = Vec2d(-r_dir[1], r_dir[0]);

	Point2d Q = Point2d(P.x + length * r_dir[0], P.y + length * r_dir[1]);

	uint p_x = uint(jclamp(0, P.x, cols - 1)), p_y = uint(jclamp(0, rows - P.y, rows - 1));
	uint q_x = uint(jclamp(0, Q.x, cols - 1)), q_y = uint(jclamp(0, rows - Q.y, rows - 1));

	uint pu_x = uint(jclamp(0, P.x + (width / 2) * r_nor[0], cols - 1));
	uint pu_y = uint(jclamp(0, rows - P.y - (width / 2) * r_nor[1], rows - 1));
	uint pl_x = uint(jclamp(0, P.x - (width / 2) * r_nor[0], cols - 1));
	uint pl_y = uint(jclamp(0, rows - P.y + (width / 2) * r_nor[1], rows - 1));

	list<Pixel> normal_axis;
	I.line(pu_y, pu_x, pl_y, pl_x, normal_axis);

	region.clear();
	region.reserve(normal_axis.size());
	for (list<Pixel>::iterator it_p = normal_axis.begin() ; it_p != normal_axis.end() ; it_p++) 
	{
		// From the translated P, finds the translated Q
		uint pc_y = it_p->first, pc_x = it_p->second;
		int dx = int(pc_x) - int(p_x), dy = int(pc_y) - int(p_y);
		uint qc_x = uint(int(q_x) + dx), qc_y = uint(int(q_y) + dy);

		// Finds the line between the two points, and eliminates pixels that do not fall in the image
		list<Pixel> line;
		I.line(pc_y, pc_x, qc_y, qc_x, line);

		list<Pixel>::iterator it_p1 = line.begin();
		Pixel front = line.front();
		if (front.first < 0 || front.first >= int(rows) || front.second < 0 || front.second >= int(cols)) {

			// If the first pixel of the line is not in the image it is needless to include this line in the region
			// since it means that, obviously, we are too far away from the initial pixel, or there is not enough pixels
			// that we can measure
			continue;

		} else {
			
			// Removes the last pixels that don't fall within the image
			Pixel last = line.back();
			if (last.first < 0 || last.first >= int(rows) || last.second < 0 || last.second >= int(cols)) {
				line.reverse();
				list<Pixel>::iterator it_p2 = line.begin();
				while (it_p2->first < 0 || it_p2->first >= int(rows) || it_p2->second < 0 || it_p2->second >= int(cols)) ++it_p2;
				line.erase(line.begin(), it_p2);
				line.reverse();
			}

			// Adds this line
			region.push_back(vector<Pixel>(line.begin(), line.end()));
		}
	}
}

