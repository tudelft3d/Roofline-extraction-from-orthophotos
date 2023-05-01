#include "partition.h"
#include <iostream>
#include <fstream>
#include <queue>
#include <ctime>
#include <set>
#include <map>
#include <numeric>
#include "defs.h"
#include "trace.h"
#include "geometry.h"
#include <cassert>
#include "svg.h"
#include <sstream>
#include <Eigen/Dense>
#include "propagation.h"

using std::set;
using std::map;
using std::min;
using std::max;

Partition::Partition() : 
	rows (0), 
	cols (0), 
	params (nullptr), 
	quadtree (nullptr), 
	id_vertices(0), 
	id_edges (0), 
	id_faces(0), 
	id_segments(0),
	vertices_head (nullptr), 
	vertices_tail(nullptr), 
	v_size(0), 
	edges_head(nullptr), 
	edges_tail(nullptr), 
	e_size(0),
	next_op(STOP),
	n_iterations(0)
{
	area = rows*cols;
	diag_length = sqrt(double(rows*rows + cols*cols));
	faces = list<Face *>();
	dual_graph = Dual_Graph();
	merge_graph = Dual_Merge_Graph();
	cdf_grad_m = vector<double>();
	cdf_lab_grad_m = vector<double>();
	energy = 0;
	generator = std::default_random_engine(time(0));
}


Partition::Partition(uint _rows, uint _cols, Parameters *_params) :
	rows (_rows),
	cols (_cols),
    params(_params)
{
	quadtree = new Quadtree(-0.5, cols-0.5, -0.5, rows-0.5);
	
	area = rows*cols;
	diag_length = sqrt(double(rows*rows + cols*cols));

	id_vertices = 0;
	id_edges = 0;
	id_faces = 0;
	id_segments = 0;

	vertices_head = NULL;
	vertices_tail = NULL;
	v_size = 0;

	edges_head = NULL;
	edges_tail = NULL;
	e_size = 0;

    faces = list<Face *>();
	dual_graph = Dual_Graph();
	merge_graph = Dual_Merge_Graph();

	cdf_grad_m = vector<double>();
	cdf_lab_grad_m = vector<double>();
	energy = 0;

	next_op = STOP;
	generator = std::default_random_engine(time(0));

	n_iterations = 0;
}


Partition::~Partition()
{
	clear();
}


void Partition::clear()
{
	// Erase all facets
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		// set all segments in it_f to be deleted
		for (auto it_s = (*it_f)->detected_segments.begin(); it_s != (*it_f)->detected_segments.end(); ++it_s)  (*it_s)->to_delete = true;
		delete (*it_f);
	}
	faces.clear();
	
	// Erase dual graphs
	dual_graph.clear();
	merge_graph.clear();
	face_to_graph.clear();
	// Erase the queue
	queue_split.clear();
	queue_merge.clear();
	queue_merge_graph.clear();
	queue_bivalent.clear();
	bboxes_by_x.clear();
	bboxes_by_y.clear();

	// Erases all edges
	Edge *e_curr = edges_head, *e_next = nullptr;
	while (e_curr != nullptr) {
		e_next = e_curr->e_next;
		delete e_curr;
		e_curr = e_next;
	}
	edges_head = edges_tail = nullptr;
	e_size = 0;
	
	// Erases all vertices and deletes quadtree
	Vertex *v_curr = vertices_head, *v_next = nullptr;
	while (v_curr != nullptr) {
		v_next = v_curr->v_next;
		delete v_curr;
		v_curr = v_next;
	}
	vertices_head = vertices_tail = nullptr;
	v_size = 0;

	if (quadtree != nullptr) delete quadtree;
	quadtree = nullptr;

	// params is set to nullptr (not destroyed here)
	params = nullptr;

	// Finally, resets the number of rows and columns
	rows = cols = 0;
	area = diag_length = 0;
	id_segments = 0;
	id_vertices = 0;
	id_edges = 0;
	id_faces = 0;

	// reset cdf
	cdf_grad_m.clear();
	cdf_lab_grad_m.clear();

	energy = 0;
	edges_next_operation.clear();
	faces_next_operation.clear();

	n_iterations = 0;
}


void Partition::copy(const Partition & P)
{
	rows = P.rows;
	cols = P.cols;
	area = P.area;
	diag_length = P.diag_length;
	offset = P.offset;
	params = P.params;

	// Copies vertices
	// First initializes a quadtree
	quadtree = new Quadtree(0, cols, 0, rows);

	// Initializes a vector of vertices that we will later use to initialize our edges
	vector<Vertex *> v_vertices;
	v_vertices.reserve(P.v_size);
	map<int, int> p_v_id_to_v_id;
	id_vertices = 0;
	vertices_head = vertices_tail = nullptr;
	v_size = 0;

	Vertex* p_v = P.vertices_head;
	while (p_v != nullptr) {
		// Deep copy of V
		Vertex* v = new Vertex(id_vertices, p_v->pt.x, p_v->pt.y);
		for (list<IndexedEvent *>::iterator it_e = p_v->events.begin() ; it_e != p_v->events.end() ; it_e++) {
			v->events.push_back(new IndexedEvent(**it_e));
		}
		quadtree->add(v);
		push_back(v);
		v_vertices.push_back(v);
		p_v_id_to_v_id.insert(pair<int, int>(p_v->id_vertex, v->id_vertex));
		// Iterates
		p_v = p_v->v_next;
	}

	// Initializes a vector of edges that we will later use to build copies of P's facets
	vector<Edge *> v_edges;
	v_edges.reserve(P.e_size);

	id_edges = 0;
	edges_head = edges_tail = nullptr;
	e_size = 0;

	Edge* p_e = P.edges_head;
	while (p_e != nullptr) {
		double alpha_v1 = p_e->v1->incidence_angle(p_e);
		double alpha_v2 = p_e->v2->incidence_angle(p_e);
		Vertex* v1 = v_vertices[p_v_id_to_v_id.at(p_e->v1->id_vertex)];
		Vertex* v2 = v_vertices[p_v_id_to_v_id.at(p_e->v2->id_vertex)];
		Edge* e = nullptr;
		if (p_e->type == INNER_EDGE) {
			Inner_Edge* casted_p_e = static_cast<Inner_Edge *>(p_e);
			set<Ray *> support;
			casted_p_e->get_supporting_rays(support);
			e = new Inner_Edge(id_edges, support, 0, v1, v2, p_e->width, alpha_v1, alpha_v2);
		} else if (p_e->type == OUTER_EDGE) {
			Outer_Edge* casted_p_e = static_cast<Outer_Edge *>(p_e);
			Image_Boundary boundary = casted_p_e->boundary;
			e = new Outer_Edge(id_edges, boundary, v1, v2, p_e->width, alpha_v1, alpha_v2, casted_p_e->v1_2);
		}
		e->normal = p_e->normal;
		e->grad_center_offset = p_e->grad_center_offset;
		e->unit_edge_cost = p_e->unit_edge_cost;
		e->mean_angle_diff = p_e->mean_angle_diff;
		e->grad_weight = p_e->grad_weight;
		e->region = p_e->region;
		push_back(e);
		v_edges.push_back(e);
		p_e = p_e->e_next;
	}

	// And now here come the facets
	id_faces = 0;
	for (list<Face *>::const_iterator it_f = P.faces.begin() ; it_f != P.faces.end() ; it_f++) {
		Face* p_f = (*it_f);
		vector<list<HalfEdge *>> & p_f_halfedges = p_f->edges;

		vector<list<HalfEdge *>> f_halfedges(p_f_halfedges.size());
		for (int id = 0; id < p_f_halfedges.size(); ++id) {
			for (list<HalfEdge *>::iterator it_h = p_f_halfedges[id].begin(); it_h != p_f_halfedges[id].end(); ++it_h) {
				HalfEdge* p_h = (*it_h);
				HalfEdge* h = (p_h->v1_v2 ? v_edges[p_h->e->id_edge]->v1_v2 : v_edges[p_h->e->id_edge]->v2_v1);
				f_halfedges[id].push_back(h);
			}
		}
		
		faces.push_back(new Face(id_faces, f_halfedges));
		faces.back()->pixels = vector<Point2i>((*it_f)->pixels);
		faces.back()->semantic_probabilities = vector<double>((*it_f)->semantic_probabilities);
		faces.back()->active_gradients = vector<set<Point2i, pointcomp>>((*it_f)->active_gradients);
		faces.back()->active_lab_gradients = vector<set<Point2i, pointcomp>>((*it_f)->active_lab_gradients);
		faces.back()->split_gain = (*it_f)->split_gain;
		faces.back()->label = (*it_f)->label;
	}

	v_edges.clear();
	v_vertices.clear();
	p_v_id_to_v_id.clear();

	cdf_grad_m = vector<double>(P.cdf_grad_m);
	cdf_lab_grad_m = vector<double>(P.cdf_lab_grad_m);
	energy = P.energy;
	id_segments = P.id_segments;

	// set iteration counter to 0
	n_iterations = 0;
}


Partition::Partition(const Partition & P)
{
	copy(P);
}


Partition & Partition::operator= (const Partition & P)
{
	if (&P != this) 
	{
		clear();
		copy(P);
	}
	return *this;
}


int & Partition::get_id_vertices()
{
	return id_vertices;
}



Vertex* Partition::erase(Vertex* v)
{	
	Vertex *v_prev = v->v_prev, *v_next = v->v_next;
	if (v_prev != nullptr) {
		v_prev->v_next = v_next;
	} else {
		vertices_head = v_next;
	}
	if (v_next != nullptr) {
		v_next->v_prev = v_prev;
	} else {
		vertices_tail = v_prev;
	}
	--v_size;
	return v_next;
}


Edge* Partition::erase(Edge* e, bool destroy)
{
	Edge *e_prev = e->e_prev, *e_next = e->e_next;
	if (e_prev != nullptr) {
		e_prev->e_next = e_next;
	} else {
		edges_head = e_next;
	}
	if (e_next != nullptr) {
		e_next->e_prev = e_prev;
	} else {
		edges_tail = e_prev;
	}
	if (destroy) delete e;
	e = nullptr;
	--e_size;
	return e_next;
}


void Partition::erase(list<Edge *> & l_e, bool destroy)
{
	for (list<Edge *>::iterator it_e = l_e.begin(); it_e != l_e.end(); it_e++) erase(*it_e, destroy);
	l_e.clear();
}



void Partition::push_back(Vertex* v)
{
	if (vertices_head == nullptr) {
		vertices_head = v;
		vertices_tail = v;
		v->v_prev = v->v_next = nullptr;
	} else {
		Vertex* v_last = vertices_tail;
		v_last->v_next = v;
		v->v_prev = v_last;
		vertices_tail = v;
	}
	++v_size;
}


void Partition::push_back(Edge* e)
{
	if (edges_head == nullptr) {
		edges_head = e;
		edges_tail = e;
		e->e_prev = e->e_next = nullptr;
	} else {
		Edge* e_last = edges_tail;
		e_last->e_next = e;
		e->e_prev = e_last;
		edges_tail = e;
	}
	++e_size;
}

void Partition::init_edges_polygon_boundary(vector<list<HalfEdge *>> & boundary_edges, vector<Ray *> & rays, map<HalfEdge *, list<Segment *>> & boundary_collinear_segments, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
	vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges, vector<Vertex *> & initial_vertices_of_rays, vector<Vertex *> & initial_outer_vertices)
{
	// To init the graph, we first add a set of vertices that represent the corners of the polygon, and build edges to connect such vertices. 
	// For each ray, we record the indices of collinear boundaries. 
	int nb_outer_edges = 0;
	for (int id = 0; id < boundary_edges.size(); ++id) nb_outer_edges += boundary_edges[id].size();
	outer_edges.reserve(nb_outer_edges);
	
	// Special case: sometimes a vertex appears multiple times in one or several boundaries
	map<Vertex*, Vertex*> corner_to_vertex;
	for (int id = 0; id < boundary_edges.size(); ++id) {
		list<HalfEdge *>::iterator it_h = boundary_edges[id].begin();
		// add the first vertex
		Vertex* corner = (*it_h)->v1_v2 ? (*it_h)->e->v1 : (*it_h)->e->v2;
		Vertex* v_start = nullptr;
		if (corner_to_vertex.count(corner) == 0) {
			v_start = new Vertex(id_vertices, corner->pt.x, corner->pt.y);
			corner_to_vertex[corner] = v_start;
			outer_vertices.push_back(v_start);
			initial_outer_vertices.push_back(v_start);
		}
		else {
			v_start = corner_to_vertex.at(corner);
		}
		int h_id = 0;
		Vertex* v_curr = v_start, *v_next = nullptr;
		for (list<HalfEdge *>::iterator it_h = boundary_edges[id].begin(); it_h != boundary_edges[id].end(); ++it_h) {
			if (*it_h != boundary_edges[id].back()) {
				corner = (*it_h)->v1_v2 ? (*it_h)->e->v2 : (*it_h)->e->v1;
				if (corner_to_vertex.count(corner) == 0) {
					v_next = new Vertex(id_vertices, corner->pt.x, corner->pt.y);
					corner_to_vertex[corner] = v_next;
					outer_vertices.push_back(v_next);
					initial_outer_vertices.push_back(v_next);
				}
				else {
					v_next = corner_to_vertex.at(corner);
				}
			}
			else {
				v_next = v_start;
			}
			Image_Boundary boundary = (id == 0) ? POLYGON_OUTER : POLYGON_INNER;
			double alpha_1 = (*it_h)->v1_v2 ? (*it_h)->e->get_alpha_1() : (*it_h)->e->get_alpha_2();
			double alpha_2 = (*it_h)->v1_v2 ? (*it_h)->e->get_alpha_2() : (*it_h)->e->get_alpha_1();
			pair<int, int> init_boundary_location(id, h_id);
			uint support_boundary_index = outer_edges.size();
			Face* f_adj = (*it_h)->opposite()->f;
			Outer_Edge* oe = new Outer_Edge(id_edges, boundary, f_adj, init_boundary_location, support_boundary_index,
				v_curr, v_next, (*it_h)->e->width, alpha_1, alpha_2, true, 0, (*it_h)->e->length);
			outer_edges.push_back(list<Outer_Edge *>(1, oe));
			v_curr = v_next;
			// record collinearity between edge and ray
			if (boundary_collinear_segments.count(*it_h) > 0) {
				list<Segment *> & collinear_segments = boundary_collinear_segments.at(*it_h);
				for (list<Segment *>::iterator it_s = collinear_segments.begin(); it_s != collinear_segments.end(); ++it_s) {
					if (!(*it_s)->is_disabled) {
						Ray* r_OA = (*it_s)->rays.first, *r_OB = (*it_s)->rays.second;
						if (r_OA != NULL) r_OA->collinear_boundaries.insert(oe->id_edge);
						if (r_OB != NULL) r_OB->collinear_boundaries.insert(oe->id_edge);
					}
				}
			}
			++h_id;
		}
	}
	for (list<Vertex *>::iterator it_v = outer_vertices.begin(); it_v != outer_vertices.end(); it_v++) quadtree->add(*it_v);

	// Then we build initial vertices that correspond to the centers of the segments.
	inner_edges.reserve(rays.size());
	initial_vertices_of_rays.reserve(rays.size());
	for (int i = 0; i < rays.size() / 2; i++) {
		Vertex* v = nullptr;
		Point2d & O = rays[2 * i]->O;
		list<Vertex *> neighbors;
		quadtree->search(O, 0.5, neighbors);
		if (!neighbors.empty()) {
			double min_distance = 100000;
			for (list<Vertex *>::iterator it_v = neighbors.begin(); it_v != neighbors.end(); ++it_v) {
				if ((*it_v)->id_vertex >= initial_outer_vertices.size()) {
					double distance = norm(Vec2d(O.x - (*it_v)->pt.x, O.y - (*it_v)->pt.y));
					if (distance <= min_distance) {
						min_distance = distance;
						v = (*it_v);
					}
				}
				
			}
			if (min_distance > params->prop_min_edge) v = nullptr;
			else {
				Point2d shift = v->pt - O;
				rays[2 * i]->A += shift;
				rays[2 * i + 1]->A += shift;
				O = v->pt;
				rays[2 * i + 1]->O = v->pt;
				rays[2 * i]->parent->finalBarycenter += shift;
				rays[2 * i]->parent->finalEnd1 += shift;
				rays[2 * i]->parent->finalEnd2 += shift;
				rays[2 * i]->parent->c = -rays[2 * i]->parent->a * rays[2 * i]->parent->finalBarycenter.x - rays[2 * i]->parent->b * rays[2 * i]->parent->finalBarycenter.y;
			}
		}
		if (v == nullptr) {
			v = new Vertex(id_vertices, O.x, O.y);
			//v->events.push_back(new IndexedEvent(2 * i, 2 * i + 1, -rays[2 * i]->initial_length, -rays[2 * i + 1]->initial_length, true));
			quadtree->add(v);
			inner_vertices.push_back(v);
		}
		// Then we insert two empty lists in the vector of edges.
		// Indeed, each entry of this vector is associated to a ray.
		inner_edges.push_back(list<Inner_Edge *>());
		inner_edges.push_back(list<Inner_Edge *>());
		// Remember the vertices we created in a vector
		initial_vertices_of_rays.push_back(v);
		initial_vertices_of_rays.push_back(v);
	}
}


void Partition::build_edge(vector<Ray *> & rays, IndexedEvent* current_event, bool & r_i_continue, Vertex* vertex, bool vertex_is_corner, bool vertex_is_new, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
	vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges, vector<Vertex *> & initial_vertices_of_rays)
{
	Ray* r_i = rays[current_event->intersectant];
	Ray* r_j = (current_event->intersected >= 0 ? rays[current_event->intersected] : NULL);
	r_i_continue = !r_i->has_been_stopped();

	// We skip the intersection of two rays on their collinear boundary (the stopping criteria still apply)
	if (current_event->is_colliding_ray) {
		if (r_i->t_collinear_boundaries.size() > 0) {
			for (list<pair<double, double>>::iterator it = r_i->t_collinear_boundaries.begin(); it != r_i->t_collinear_boundaries.end(); ++it) {
				double t1 = it->first, t2 = it->second;
				if (current_event->t_intersectant >= t1 && current_event->t_intersectant <= t2) {
					if (vertex_is_new) {
						quadtree->remove(vertex, true);
					}
					return;
				}
			}
		}
		if (r_j->t_collinear_boundaries.size() > 0) {
			for (list<pair<double, double>>::iterator it = r_j->t_collinear_boundaries.begin(); it != r_j->t_collinear_boundaries.end(); ++it) {
				double t1 = it->first, t2 = it->second;
				if (current_event->t_intersected >= t1 && current_event->t_intersected <= t2) {
					if (vertex_is_new) {
						quadtree->remove(vertex, true);
					}
					return;
				}
			}
		}
	}

	// First, the intersected ray or boundary : we split the corresponding edge if a new vertex is to be added
	double abscissa;
	
	bool overlap_with_outer_edge = false, overlap_with_inner_edge = false;

	if (!current_event->is_colliding_ray && !vertex_is_corner) {
		
		// On one hand, if the interesected ray is an image boundary
		if (vertex_is_new) outer_vertices.push_back(vertex);

		// First, we determine which edge must be split
		Image_Boundary intersected_boundary = Image_Boundary(current_event->intersected);
		
		list<Outer_Edge *> & edges_for_this_corner = outer_edges[current_event->intersected_index];
		// set the reference point O for time = 0, which is the end vertex of the original outer edge
		Point2d & O = edges_for_this_corner.back()->v2->pt;

		if (intersected_boundary == TOP_IMAGE || intersected_boundary == BOTTOM_IMAGE) {
			abscissa = vertex->pt.x;
		} else if (intersected_boundary == LEFT_IMAGE || intersected_boundary == RIGHT_IMAGE) {
			abscissa = vertex->pt.y;
		}
		else if (intersected_boundary == POLYGON_OUTER || intersected_boundary == POLYGON_INNER) {
			abscissa = -cv::norm(vertex->pt - O);
		}

		// We loop on the list of edges associated to this boundary
		list<Outer_Edge *>::iterator it1 = edges_for_this_corner.begin();
		double t1, t2;
		while (it1 != edges_for_this_corner.end()) {
			if (intersected_boundary == POLYGON_OUTER || intersected_boundary == POLYGON_INNER) {
				(*it1)->time_range(t1, t2, O, current_event->intersected_index);
			}
			else {
				(*it1)->time_range(t1, t2, Point2d(0,0), current_event->intersected_index);
			}
			if (it1 == edges_for_this_corner.begin()) {
				t1 -= 0.1;
			}
			else if (it1 == std::prev(edges_for_this_corner.end())) {
				t2 += 0.1;
			}
			if (t1 <= abscissa && abscissa < t2) break;
			it1++;
		}
		// The loop has stopped on the edge to split
		Outer_Edge* split_edge = *it1;
		Vertex* v1 = split_edge->v1, *v2 = split_edge->v2;

		if (vertex != v1 && vertex != v2) {
			double alpha_1 = v1->incidence_angle(split_edge), alpha_2 = v2->incidence_angle(split_edge);
			bool create_v1_v = true, create_v_v2 = true;
			//for (auto it_dir = vertex->directions.begin(); it_dir != vertex->directions.end(); it_dir++) {
			//	if (it_dir->second->e->v1 != v1 && it_dir->second->e->v2 != v1) create_v1_v = false;
			//	if (it_dir->second->e->v1 != v2 && it_dir->second->e->v2 != v2) create_v_v2 = false;
			//}
			it1 = edges_for_this_corner.erase(it1);
			if (create_v1_v) {
				Outer_Edge* v1_v = new Outer_Edge(id_edges, intersected_boundary, split_edge->outer_face, split_edge->init_boundary_location, current_event->intersected_index,
					v1, vertex, split_edge->width, alpha_1, alpha_2, split_edge->v1_2, split_edge->dir_v1v2, t1, abscissa);
				edges_for_this_corner.insert(it1, v1_v);
			}
			if (create_v_v2) {
				Outer_Edge* v_v2 = new Outer_Edge(id_edges, intersected_boundary, split_edge->outer_face, split_edge->init_boundary_location, current_event->intersected_index,
					vertex, v2, split_edge->width, alpha_1, alpha_2, split_edge->v1_2, split_edge->dir_v1v2, abscissa, t2);
				edges_for_this_corner.insert(it1, v_v2);
			}
			
			delete split_edge;
		}
	} else if (current_event->is_colliding_ray) {

		// On the other hand, if the intersected ray is a real ray
		if (vertex_is_new && !vertex_is_corner) inner_vertices.push_back(vertex);

		abscissa = current_event->t_intersected;
		list<Inner_Edge *> & edges_for_this_intersected_ray = inner_edges[current_event->intersected];
		if (edges_for_this_intersected_ray.size() == 0) {
			Vertex* v0 = initial_vertices_of_rays[current_event->intersected];
			// In case no edge modelizes the path followed by the intersected ray, we create one
			if (v0 != vertex) {
				double alpha_1, alpha_2;
				if (vertex_is_new) {
					alpha_1 = r_j->opposite()->incidence_angle; 
					alpha_2 = r_j->incidence_angle;
				}
				else {
					alpha_1 = atan2(vertex->pt.y - v0->pt.y, vertex->pt.x - v0->pt.x);
					alpha_2 = (alpha_1 <= 0 ? alpha_1 + PI : alpha_1 - PI);
				}
				Inner_Edge* v1_v = new Inner_Edge(id_edges, r_j, 0, v0, vertex, r_j->parent->width, alpha_1, alpha_2, -r_j->initial_length, abscissa);
				edges_for_this_intersected_ray.push_back(v1_v);
			}

		} else {

			// Loops on the list of edges associated to the intersected ray
			// Exhibits the one that contains the intersection point, if it does exist
			Vertex* v1 = NULL;
			Vertex* v2 = NULL;
			double last_abscissa = -FLT_MAX;
			list<Inner_Edge *>::iterator it2 = edges_for_this_intersected_ray.begin();
			while (it2 != edges_for_this_intersected_ray.end()) {
				double t1, t2;
				(*it2)->time_range(r_j, t1, t2);
				if (t1 <= abscissa && abscissa < t2) break;
				v2 = (*it2)->v2;
				last_abscissa = t2; // (*it2)->t2;
				it2++;
			}

			if (it2 == edges_for_this_intersected_ray.end()) {
				if (vertex == initial_vertices_of_rays[current_event->intersected]) {
					// in case of rounding errors
					v2 = vertex;
				}
				// Adds a new edge at the end of the list if it does not overlap with the polygon boundary, otherwise tag it as a ghost edge
				if (v2 != vertex) {

					// Check if v2-->vertex overlaps with an outer edge, or an inner edge, if vertex is a rounded vertex
					for (int i = 0; i < vertex->directions.size(); ++i) {
						Edge* e_i = vertex->directions[i].second->e;
						if (e_i->type == OUTER_EDGE && (e_i->v1 == v2 || e_i->v2 == v2)) {
							overlap_with_outer_edge = true;
						}
						else if (e_i->type == INNER_EDGE && (e_i->v1 == v2 || e_i->v2 == v2)) {
							overlap_with_inner_edge = true;
						}
					}
					int tag = overlap_with_outer_edge ? -1 : (overlap_with_inner_edge? -2 : 0);

					double alpha_1 = atan2(vertex->pt.y - v2->pt.y, vertex->pt.x - v2->pt.x);
					double alpha_2 = (alpha_1 <= 0 ? alpha_1 + PI : alpha_1 - PI);

					Inner_Edge* v2_v = new Inner_Edge(id_edges, r_j, tag, v2, vertex, r_j->parent->width, alpha_1, alpha_2, last_abscissa, abscissa, tag>=0);
					edges_for_this_intersected_ray.push_back(v2_v);
				}
			} else {
				// Splits this edge and tag child edges if it overlaps with the polygon boundary or existing inner edges
				Inner_Edge* split_edge = *it2;
				Vertex* v1 = split_edge->v1, *v2 = split_edge->v2;
				double t1, t2;
				split_edge->time_range(r_j, t1, t2);
				assert(v1 != v2);
				if (v1 != vertex && v2 != vertex) {
					Inner_Edge* v1_v, *v_v2;
					if (vertex_is_new) {
						////double alpha_1 = v1->incidence_angle(split_edge), alpha_2 = v2->incidence_angle(split_edge);
						double alpha_v1 = atan2(vertex->pt.y - v1->pt.y, vertex->pt.x - v1->pt.x);
						double alpha_v_v1 = (alpha_v1 <= 0 ? alpha_v1 + PI : alpha_v1 - PI);
						double alpha_v_v2 = atan2(v2->pt.y - vertex->pt.y, v2->pt.x - vertex->pt.x);
						double alpha_v2 = (alpha_v_v2 <= 0 ? alpha_v_v2 + PI : alpha_v_v2 - PI);

						overlap_with_outer_edge = (split_edge->tag == -1);
						overlap_with_inner_edge = (split_edge->tag == -2);

						v1_v = new Inner_Edge(id_edges, r_j, split_edge->tag, v1, vertex, r_j->parent->width, alpha_v1, alpha_v_v1, t1, abscissa);
						v_v2 = new Inner_Edge(id_edges, r_j, split_edge->tag, vertex, v2, r_j->parent->width, alpha_v_v2, alpha_v2, abscissa, t2);
					}
					else {
						// The split_edge could have a rounded vertex, therefore we recompute alpha
						double alpha_v1 = atan2(vertex->pt.y - v1->pt.y, vertex->pt.x - v1->pt.x);
						double alpha_v_v1 = (alpha_v1 <= 0 ? alpha_v1 + PI : alpha_v1 - PI);
						double alpha_v_v2 = atan2(v2->pt.y - vertex->pt.y, v2->pt.x - vertex->pt.x);
						double alpha_v2 = (alpha_v_v2 <= 0 ? alpha_v_v2 + PI : alpha_v_v2 - PI);
						// Check if v1-->vertex and/or vertex-->v2 overlaps with an outer edge
						bool overlap_with_outer_edge_v1_v = false, overlap_with_outer_edge_v_v2 = false;
						bool overlap_with_inner_edge_v1_v = false, overlap_with_inner_edge_v_v2 = false;
						for (int i = 0; i < vertex->directions.size(); ++i) {
							Edge* e_i = vertex->directions[i].second->e;
							if (e_i->type == OUTER_EDGE && (e_i->v1 == v1 || e_i->v2 == v1)) {
								overlap_with_outer_edge_v1_v = true; break;
							} else if (e_i->type == INNER_EDGE && (e_i->v1 == v1 || e_i->v2 == v1)) {
								overlap_with_inner_edge_v1_v = true; break;
							}
						}
						for (int i = 0; i < vertex->directions.size(); ++i) {
							Edge* e_i = vertex->directions[i].second->e;
							if (e_i->type == OUTER_EDGE && (e_i->v1 == v2 || e_i->v2 == v2)) {
								overlap_with_outer_edge_v_v2 = true; break;
							}
							else if (e_i->type == INNER_EDGE && (e_i->v1 == v2 || e_i->v2 == v2)) {
								overlap_with_inner_edge_v_v2 = true; break;
							}
						}
						int tag_v1_v = overlap_with_outer_edge_v1_v ? -1:  (overlap_with_inner_edge_v1_v ? -2 : 0);
						int tag_v_v2 = overlap_with_outer_edge_v_v2 ? -1: (overlap_with_inner_edge_v_v2 ? -2 : 0);
						
						overlap_with_outer_edge = (overlap_with_outer_edge_v1_v || overlap_with_outer_edge_v_v2);
						overlap_with_inner_edge = (overlap_with_inner_edge_v1_v || overlap_with_inner_edge_v_v2);

						v1_v = new Inner_Edge(id_edges, r_j, tag_v1_v, v1, vertex, r_j->parent->width, alpha_v1, alpha_v_v1, t1, abscissa, tag_v1_v >= 0);
						v_v2 = new Inner_Edge(id_edges, r_j, tag_v_v2, vertex, v2, r_j->parent->width, alpha_v_v2, alpha_v2, abscissa, t2, tag_v_v2 >= 0);
					}
					it2 = edges_for_this_intersected_ray.erase(it2);
					edges_for_this_intersected_ray.insert(it2, v1_v);
					edges_for_this_intersected_ray.insert(it2, v_v2);
					delete split_edge;
				}

			}
		}
	}

	// Second, the intersectant ray : we create an edge that comes from the last intersection point
	// in which it was actively involved, to the interesection point we just created
	// In the special case where the intersectant ray overlaps a parallel outer edge on polygon inner contours, we create a ghost edge.
	// In the special case where the intersected ray overlaps an outer edge or an inner edge, we create nothing.
	if (overlap_with_outer_edge || overlap_with_inner_edge) r_i_continue = true;

	Vertex* v3 = NULL;
	Inner_Edge* last_e = NULL;
	Inner_Edge* v3_v = NULL;
	abscissa = current_event->t_intersectant;
	double alpha_v3 = r_i->opposite()->incidence_angle;
	double alpha_v = r_i->incidence_angle;

	list<Inner_Edge *> & edges_for_this_intersectant_ray = inner_edges[current_event->intersectant];
	list<Inner_Edge *>::iterator it2 = edges_for_this_intersectant_ray.begin();

	if (edges_for_this_intersectant_ray.size() == 0) {
		v3 = initial_vertices_of_rays[current_event->intersectant];

		if (v3 != vertex) {
			if (!vertex_is_new) {
				alpha_v3 = atan2(vertex->pt.y - v3->pt.y, vertex->pt.x - v3->pt.x);
				alpha_v = (alpha_v3 <= 0 ? alpha_v3 + PI : alpha_v3 - PI);
			}
			v3_v = new Inner_Edge(id_edges, r_i, 0, v3, vertex, r_i->parent->width, alpha_v3, alpha_v, -r_i->initial_length, abscissa);
		}
	}
	else {
		double last_abscissa = -FLT_MAX;
		while (it2 != edges_for_this_intersectant_ray.end()) {
			last_e = *it2;
			v3 = last_e->v2;
			double t1, t2;
			(*it2)->time_range(r_i, t1, t2);
			last_abscissa = t2;
			it2++;
		}
		if (v3 != vertex) {
			// Check if v3-->vertex overlaps with an outer edge (i.e. check if r_i is passing through a collinear boundary)
			bool intersectant_overlap_with_outer_edge = false, intersectant_overlap_with_inner_edge = false;
			for (int k = 0; k < vertex->directions.size(); ++k) {
				Edge* e = vertex->directions[k].second->e;
				if (e->type == OUTER_EDGE) {
					Outer_Edge* oe = static_cast<Outer_Edge *>(e);
					list<Outer_Edge *> & edges_for_this_corner = outer_edges[oe->support_boundary_index];
					for (list<Outer_Edge *>::iterator it_oe = edges_for_this_corner.begin(); it_oe != edges_for_this_corner.end(); ++it_oe) {
						if ((*it_oe)->v1 == v3 || (*it_oe)->v2 == v3) {
							intersectant_overlap_with_outer_edge = true; break;
						}
					}
				}
				else if (e->type == INNER_EDGE) {
					Inner_Edge* ie = static_cast<Inner_Edge *>(e);
					if (ie->tag == -1) continue;
					if (e->v1 == v3 || e->v2 == v3) {
						intersectant_overlap_with_inner_edge = true; break;
					}
					list<Inner_Edge *> edges_for_this_direction;
					for (set<Ray *>::iterator it_r = ie->rays.begin(); it_r != ie->rays.end(); ++it_r) {
						list<Inner_Edge *> & edges_for_this_ray = inner_edges[(*it_r)->index];
						list<Inner_Edge *> & edges_for_opp_ray = inner_edges[(*it_r)->opposite()->index];
						edges_for_this_direction.insert(edges_for_this_direction.begin(), edges_for_this_ray.begin(), edges_for_this_ray.end());
						edges_for_this_direction.insert(edges_for_this_direction.begin(), edges_for_opp_ray.begin(), edges_for_opp_ray.end());
					}
					for (list<Inner_Edge *>::iterator it_ie = edges_for_this_direction.begin(); it_ie != edges_for_this_direction.end(); ++it_ie) {
						if ((*it_ie)->v1 == v3 || (*it_ie)->v2 == v3) {
							intersectant_overlap_with_inner_edge = true; break;
						}
					}
					if (intersectant_overlap_with_inner_edge) break;
				}
			}

			if (!vertex_is_new || fabs(last_e->get_alpha_1() - alpha_v3) > 1e-7 ) {
				alpha_v3 = atan2(vertex->pt.y - v3->pt.y, vertex->pt.x - v3->pt.x);
				alpha_v = (alpha_v3 <= 0 ? alpha_v3 + PI : alpha_v3 - PI);
			}

			// If v3-->vertex overlaps with an outer edge, we create a "ghost" inner edge and tag it to be deleted in the end.
			int tag = intersectant_overlap_with_outer_edge ? -1 : (intersectant_overlap_with_inner_edge ? -2 : 0);
			v3_v = new Inner_Edge(id_edges, r_i, tag, v3, vertex, r_i->parent->width, alpha_v3, alpha_v, last_abscissa, abscissa, tag >= 0);
		}
	}
	if (v3_v != NULL) edges_for_this_intersectant_ray.push_back(v3_v);
}

void Partition::find_pixels_inside_all_facets()
{
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		(*it_f)->find_pixels_inside_facet(offset);
	}
}

void Partition::compute_feature_all_facets(Matrix<double> & I_prob, Matrix<uchar> & I)
{
	for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
		(*it_f)->compute_feature(I_prob, I, offset, params->max_label);
	}
}

void Partition::find_active_gradients_all_facets(Matrix<double> & I_m, Matrix<double> & I_t, bool is_color_image)
{
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		(*it_f)->find_active_gradients(I_m, I_t, offset, params->split_params, is_color_image);
	}
}

void Partition::weight_all_edges(Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t)
{
	assert(cdf_grad_m.size() != 0);
	Edge* e = edges_head;
	while (e != NULL) {
		e->weight(used, I_m, I_t, offset, params->split_params, cdf_grad_m);
		e = e->e_next;
	}
}

void Partition::estimate_cdf_grad_m(Matrix<double> & I_m, bool is_color_image)
{
	double bin_width = params->split_params.bin_width;
	int n_bins = ceil(0.5 / bin_width);

	vector<double> & cdf = is_color_image ? cdf_grad_m : cdf_lab_grad_m;
	cdf.clear();
	cdf.resize(n_bins, 0);

	for (uint i = 0; i < I_m.rows; i++) {
		for (uint j = 0; j < I_m.cols; j++) {
			double value_p = I_m(i, j);
			int id_bin = floor(value_p / bin_width);
			id_bin = jclamp(0, id_bin, n_bins - 1);
			cdf[id_bin] += 1;
		}
	}
	int n_pixels = I_m.rows*I_m.cols;
	for (int i = 1; i < n_bins; ++i) cdf[i] += cdf[i-1];
	for (int i = 0; i < n_bins; ++i) cdf[i] /= cdf.back();

	double & m_value_thresh = is_color_image ? params->split_params.m_value_thresh : params->split_params.lab_m_value_thresh;
	m_value_thresh = 1*bin_width;
	for (int i = 0; i < n_bins; ++i) {
		if (cdf[i] < 1 - params->split_params.gradient_percentage) m_value_thresh = (i+1)*bin_width;
		else break;
	}
}

void Partition::extend_detected_segments_all_facets(Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, 
	Matrix<double> & I_m, Matrix<double> & I_t, map<int, double> & angles, vector<int> & bin_to_angle_index)
{
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		extend_detected_segments((*it_f), I, I_prob, used, I_m, I_t, angles, bin_to_angle_index, (*it_f)->split_gain);
		(*it_f)->delete_invalid_segments();
	}
}

void Partition::merge_containers(list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices, vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges)
{
	e_size = 0;
	// We are going to insert elements contained in outer_vertices, inner_vertices in a double-linked list
	for (list<Vertex *>::iterator it_v = outer_vertices.begin(); it_v != outer_vertices.end(); it_v++) push_back(*it_v);
	outer_vertices.clear();
	for (list<Vertex *>::iterator it_v = inner_vertices.begin(); it_v != inner_vertices.end(); it_v++) push_back(*it_v);
	inner_vertices.clear();

	// We also merge the different lists of edges
	// Step 1 : inner_edges
	// In particular, we eliminate the "shadow" inner edges that overlaps with polygon inner contour
	for (int i = 0 ; i < inner_edges.size() ; i++) {
		for (list<Inner_Edge *>::iterator it_e = inner_edges[i].begin(); it_e != inner_edges[i].end(); it_e++) {
			if ((*it_e)->tag < 0) {
				delete *it_e;
			}
			else {
				push_back(*it_e);
			}
		}
		inner_edges[i].clear();
	}
	inner_edges.clear();

	// Step 2 : outer_edges
	for (int i = 0; i < outer_edges.size(); i++) {
		for (list<Outer_Edge *>::iterator it_e = outer_edges[i].begin(); it_e != outer_edges[i].end(); it_e++) {
			push_back(*it_e);
		}
	}

	// cleanup repeated connection
	Vertex* v = vertices_head;
	while (v != NULL) {
		set<Vertex *> connected_vertices;
		for (int i = 0; i < v->directions.size(); ++i) {
			Edge* e_i = v->directions[i].second->e;
			Vertex* connect_v = e_i->v1 == v ? e_i->v2 : e_i->v1;
			if (connected_vertices.count(connect_v) != 0) {
				if (e_i->type == INNER_EDGE) erase(e_i, true);
				else {
					Edge* e_to_remove = v->directions[i-1].second->e;
					Vertex* v_to_remove = e_to_remove->v1 == v ? e_to_remove->v2 : e_to_remove->v1;
					assert(v_to_remove == connect_v);
					erase(e_to_remove, true);
				}
			}
			connected_vertices.insert(connect_v);
		}
		v = v->v_next;
	}
}


void Partition::build_faces()
{// build non-nested faces
	id_edges = 0;
	int nb_edges = int(e_size);
	bool** queued = new bool*[nb_edges];
	map<int, HalfEdge*> contours;

	Edge* e = edges_head;
	while (e != NULL) {
		queued[id_edges] = new bool[2];
		// The algorithm of face retrieval consists in crossing over the edges of the
		// graph in clockwise order to find cycles. This leads us to consider virtual
		// half-edges : for every e = (v1 v2), there "exists" the half-edge (v1 v2)
		// and the half-edge (v2 v1), except for the outer edges for which only one
		// half-edge exists. We eliminate such half-edges from the algorithm by
		// considering that they have been already processed.

		bool top_image = false, bottom_image = false, right_image = false, left_image = false, polygon_outer_edge = false, polygon_inner_edge = false;
		if (e->type == OUTER_EDGE) {
			Outer_Edge* oe = static_cast<Outer_Edge *>(e);
			top_image = (oe->boundary == TOP_IMAGE);
			bottom_image = (oe->boundary == BOTTOM_IMAGE);
			right_image = (oe->boundary == RIGHT_IMAGE);
			left_image = (oe->boundary == LEFT_IMAGE);
			polygon_outer_edge = (oe->boundary == POLYGON_OUTER);
			polygon_inner_edge = (oe->boundary == POLYGON_INNER);
			int boundary_id = oe->init_boundary_location.first;
			if (!contours.count(boundary_id)) {
				if (bottom_image || right_image) {
					contours[boundary_id] = oe->v2_v1;
				}
				else {
					contours[boundary_id] = oe->v1_v2;
				}
			}
		}

		if (top_image || left_image) {
			// Entry 1 correspond to the direction (v1 v2), entry 0 to (v2 v1)
			// Values are determined according to the way outer edges are created
			queued[id_edges][1] = false;
			queued[id_edges][0] = true;
		} else if (bottom_image || right_image) {
			queued[id_edges][1] = true;
			queued[id_edges][0] = false;
		} else if (polygon_outer_edge || polygon_inner_edge) {
			queued[id_edges][1] = false;
			queued[id_edges][0] = true;
		}
		else {
			// General case : inner edge
			queued[id_edges][1] = false;
			queued[id_edges][0] = false;
		}
		e->id_edge = id_edges++;
		e = e->e_next;
	}

	// Defines a queue and initializes with the top left half-edge of the figure, or a vertex on the outer boundary of the polygon
	std::queue<HalfEdge *> queue;
	Vertex* v_start = vertices_head;
	HalfEdge* h = v_start->directions[v_start->directions.size() - 1].second;
	// If the chosen halfedge is not an unvisited halfedge on the outer boundary, we need to choose another direction
	int i_dir = 0;
	while (queued[h->e->id_edge][h->v1_v2]) {
		h = v_start->directions[i_dir].second;
		++i_dir;
	}
	queue.push(h);

	list<vector<list<HalfEdge*>>> faces_edges;
	set<int> contours_visited;
	// While there are elements left in the queue
	while (queue.size() != 0) {
		// We define a reference to the first element of the queue,
		// which is the first element of a new face
		HalfEdge* h_e = queue.front();
		HalfEdge* h_f = h_e;
		vector<list<HalfEdge *>> current_face(1);

		// If the element hasn't been already processed in the meanwhile
		if (!queued[h_f->e->id_edge][h_f->v1_v2]) {
			queued[h_f->e->id_edge][h_f->v1_v2] = true;

			do {
				// Finds the next half-edge in clockwise order
				// Adds it to the face we are currently defining
				Vertex* v = (h_f->v1_v2 ? h_f->e->v2 : h_f->e->v1);
				int path = 0;
				for (int i = 0; i < v->directions.size(); i++) {
					if (v->directions[i].second->e == h_f->e) {
						path = (i + 1) % v->directions.size();
						break;
					}
				}
				h_f = v->directions[path].second;
				assert((h_f->v1_v2 && h_f->e->v1 == v) || (!h_f->v1_v2 && h_f->e->v2 == v));
				current_face[0].push_back(h_f);
				if (h_f->e->type == OUTER_EDGE) {
					Outer_Edge* oe = static_cast<Outer_Edge *>(h_f->e);
					contours_visited.insert(oe->init_boundary_location.first);
				}

				// The current half-edge shouldn't be queued
				queued[h_f->e->id_edge][h_f->v1_v2] = true;

				// Inserts the opposite half-edge into the queue
				if (!queued[h_f->e->id_edge][!h_f->v1_v2]) {
					if (h_f->e->v1_v2 == h_f) {
						queue.push(h_f->e->v2_v1);
					} else {
						queue.push(h_f->e->v1_v2);
					}
				}

			} while (h_f != h_e);
			faces_edges.push_back(current_face);
		}

		// Removes the first element of the queue, since we're done with it
		queue.pop();
	}

	// Add unvisited inner contours to the corresponding facet
	for (const auto& id_he_pair: contours) {
		if (id_he_pair.first==0 || contours_visited.count(id_he_pair.first) != 0) continue;

		Point2d pt = 0.5*(id_he_pair.second->e->v1->pt + id_he_pair.second->e->v2->pt);
		// search for the facet that contains this inner contour
		for (vector<list<HalfEdge*>>& f_edges : faces_edges) {
			if (Geometry::is_inside_contour(pt, f_edges[0])) {
				// form the inner contour in anti-clockwise order
				list<HalfEdge *> inner_contour;
				HalfEdge* h_e = id_he_pair.second;
				HalfEdge* h_f = h_e;
				do {
					inner_contour.push_back(h_f);
					Vertex* v = (h_f->v1_v2 ? h_f->e->v2 : h_f->e->v1);
					int path = 0;
					for (int i = 0; i < v->directions.size(); i++) {
						if (v->directions[i].second->e == h_f->e) {
							path = (i - 1) % v->directions.size();
							break;
						}
					}
					h_f = v->directions[path].second;
					assert(h_f->e->type == OUTER_EDGE);
				} while (h_f != h_e);
				f_edges.push_back(inner_contour);
				break;
			}
		}
	}

	for (vector<list<HalfEdge*>>& f_edges : faces_edges) {
		faces.push_back(new Face(id_faces, f_edges));
	}

	for (int i = 0; i < int(e_size); i++) {
		delete[] queued[i];
	}
	delete[] queued;

	// debug
    //bool simple = (count_nonsimple_faces() == 0);
}


void Partition::build_split_faces(const set<Edge*>& split_edges, Face* parent_face, vector<Face *> & split_facets) {
	
	map<int, bool*> queued;
	set<Edge*> edges_to_visit(split_edges);

	for (Edge* e : split_edges) {
		queued[e->id_edge] = new bool[2];
		// The algorithm of face retrieval consists in crossing over the edges of the
		// graph in clockwise order to find cycles, as well as the hierachy of cycles
		// Entry 1 correspond to the direction (v1 v2), entry 0 to (v2 v1)
		queued[e->id_edge][1] = false;
		queued[e->id_edge][0] = false;
		e = e->e_next;
	}

	map<HalfEdge*, int> boundary_edge_contour_id;
	for (size_t i=0; i<parent_face->edges.size(); ++i) {
		for (HalfEdge* he : parent_face->edges[i]) {
			edges_to_visit.insert(he->e);
			boundary_edge_contour_id[he] = i;
			queued[he->e->id_edge] = new bool[2];
			if (he->v1_v2) {
				queued[he->e->id_edge][1] = false;
				queued[he->e->id_edge][0] = true;
			} else {
				queued[he->e->id_edge][1] = true;
				queued[he->e->id_edge][0] = false;
			}
		}
	}

	vector<list<HalfEdge*>> inner_contours; 
	vector<vector<list<HalfEdge*>>> new_faces_edges;
	while (edges_to_visit.size() != 0) {
		std::queue<HalfEdge *> queue;
		HalfEdge* h = (*edges_to_visit.begin())->v1_v2;
		if (queued[h->e->id_edge][h->v1_v2]) h = h->opposite(); 
		queue.push(h);
		edges_to_visit.erase(h->e);

		while (queue.size() != 0) {
			// We define a reference to the first element of the queue,
			// which is the first element of a new face
			HalfEdge* h_e = queue.front();
			HalfEdge* h_f = h_e;
			list<HalfEdge *> current_contour;

			// If the element hasn't been already processed in the meanwhile
			if (!queued[h_f->e->id_edge][h_f->v1_v2]) {
				queued[h_f->e->id_edge][h_f->v1_v2] = true;

				do {
					// Attempt to find the next half-edge in clockwise order
					// Adds it to the face we are currently defining
					Vertex* v = (h_f->v1_v2 ? h_f->e->v2 : h_f->e->v1);
					int path = 0;
					for (int i = 0; i < v->directions.size(); i++) {
						if (v->directions[i].second->e == h_f->e) {
							path = (i + 1) % v->directions.size();
							break;
						}
					}
					h_f = v->directions[path].second;
					assert((h_f->v1_v2 && h_f->e->v1 == v) || (!h_f->v1_v2 && h_f->e->v2 == v));
					current_contour.push_back(h_f);
					edges_to_visit.erase(h_f->e);

					// The current half-edge shouldn't be queued
					queued[h_f->e->id_edge][h_f->v1_v2] = true;

					// Inserts the opposite half-edge into the queue
					if (!queued[h_f->e->id_edge][!h_f->v1_v2]) {
						if (h_f->e->v1_v2 == h_f) {
							queue.push(h_f->e->v2_v1);
						} else {
							queue.push(h_f->e->v1_v2);
						}
					}
				} while (h_f != h_e);
				if (HalfEdge::is_clockwise(current_contour)) {
					new_faces_edges.emplace_back(1, current_contour);
				} else { // store this cycle
					inner_contours.emplace_back(current_contour);
				}
			}
			// Removes the first element of the queue, since we're done with it
			queue.pop();
		}
	}

	// delete parent facet before constructing new facets
	faces.remove(parent_face);
	delete parent_face;

	for (vector<list<HalfEdge*>> face_edges : new_faces_edges) {
		Face* f_new = new Face(id_faces, face_edges);
		faces.push_back(f_new);
		split_facets.push_back(f_new);
	}

	// Build hierarchy if inner contours are formed. Each contour points to its inner contours that are one-layer lower
	map<Face*, std::unordered_set<Face*>> lower_layer;
	map<Face*, int> inner_face_to_inner_contour_id;
	if (!inner_contours.empty()) {
		for (size_t i=0; i<inner_contours.size(); ++i) {
			std::unordered_set<Face*> inner_faces;
			for (HalfEdge* he : inner_contours[i]) {
				assert(he->opposite()->f != NULL);
				inner_faces.insert(he->opposite()->f);
				inner_face_to_inner_contour_id[he->opposite()->f] = i;
			}
			HalfEdge* he = inner_contours[i].front();
			Point2d pt = 0.5*(he->e->v1->pt + he->e->v2->pt);
			for (Face* f : split_facets) {
				if (inner_faces.count(f)) continue;
				if (Geometry::is_inside_contour(pt, f->edges[0])) {
					if (lower_layer.count(f)) {
						lower_layer.at(f).insert(inner_faces.begin(), inner_faces.end());
					} else {
						lower_layer[f] = inner_faces;
					}
				}
			}
		}
		// clean: remove nested faces inside inner contour
		for (auto& layer_pair : lower_layer) {
			std::unordered_set<Face*>& inner_faces = layer_pair.second;
			auto it_inner = inner_faces.begin();
			while (it_inner != inner_faces.end()) {
				Face* f_inner = *it_inner;
				if (lower_layer.count(f_inner)) {
					for (Face* f_nested : lower_layer.at(f_inner)) {
						assert(f_nested!=f_inner);
						inner_faces.find(f_nested);
					}
				}
				++it_inner;
			}
		}
	}
	// Add inner contours to corresponding faces
	for (auto& layer_pair : lower_layer) {
		Face* f = layer_pair.first;
		std::unordered_set<int> inner_contours_id;
		std::unordered_set<Face*>& inner_faces = layer_pair.second;
		for (Face* f_inner : inner_faces) {
			inner_contours_id.insert(inner_face_to_inner_contour_id.at(f_inner));
		}
		if (!inner_contours_id.empty()) {
			vector<list<HalfEdge*>> f_edges_new(f->edges);
			f_edges_new.reserve(1+inner_contours_id.size());
			for (int i : inner_contours_id) {
				f_edges_new.push_back(inner_contours[i]);
			}
			f->edges = f_edges_new;
			for (size_t id = 0; id < f->edges.size(); ++id) {
				for (list<HalfEdge *>::iterator it = f->edges[id].begin(); it != f->edges[id].end(); it++) {
					(*it)->set(f);
				}
			}
			f->list_vertices();
			f->classify();
		}
	}

	for (auto& qpair: queued) {
		delete[] qpair.second;
	}

	// debug
    //bool valid = (count_nonsimple_faces() == 0);
}


Edge* Partition::merge_collinear_edges(Vertex* v, list<Vertex*>& vertices_erased, bool destroy) {

	Edge* e_12 = nullptr;
	if (v == nullptr) return nullptr;

	// The 4 first vertices are image corners
	if (v->connectivity() == 2 && v->id_vertex > 3) {
		Point2d p = v->pt;
		list<Edge *> edges_to_delete;
		set<Ray *> rays;
		Image_Boundary boundary = INVALID_BORDER;
		bool v1_2 = true;

		double alpha_eps = 1e-4;
		double dist_eps = 1e-3;

		// The vertex v is bivalent : that's why we find the vertices to which v is connected
		// However, maybe the connected vertices are themselves bivalent : we need to iterate on v1 and v2 until
		// finding trivalent vertices. Until then, we define vertices and edges to delete
		Edge* e1 = v->directions[0].second->e;
		Edge* e2 = v->directions[1].second->e;

		Edge_Type e_type = e1->type;
		Inner_Edge *ie = nullptr;
		Outer_Edge* oe = nullptr;

		Vertex* v1 = (e1->v1 == v ? e1->v2 : e1->v1);
		Vertex* v2 = (e2->v1 == v ? e2->v2 : e2->v1);

		// We define a sequence of vertices (v1 .. v .. v2)
		// Obtains the two facets adjacent to the big edge (v1 v2) : to this end, we can use e2
		Face *f_v1_v2 = nullptr, *f_v2_v1 = nullptr;
		HalfEdge* h = e2->v1_v2;
		if (h->e->v1 == v && h->e->v2 == v2) {
			f_v1_v2 = h->f;
			f_v2_v1 = h->opposite()->f;
		}
		else {
			f_v1_v2 = h->opposite()->f;
			f_v2_v1 = h->f;
		}

		if (e_type == INNER_EDGE) {
			ie = static_cast<Inner_Edge*>(e1);
			ie->get_supporting_rays(rays);
			ie = static_cast<Inner_Edge*>(e2);
			ie->get_supporting_rays(rays);
		}
		else if (e_type == OUTER_EDGE) {
			oe = static_cast<Outer_Edge*>(e1);
			boundary = oe->boundary;
			v1_2 = oe->v1_2;
		}
		else {
			// ARTIFICIAL_EDGE
			Artificial_Edge* ae = static_cast<Artificial_Edge *>(e1);
		}

		// test if v is collinear with v1--v2
		bool collinear_v, collinear_1, collinear_2;
		double alpha_v = e1->get_alpha_1();

		Point2d line_dir = v2->pt - v1->pt;
		line_dir /= cv::norm(line_dir);
		Point2d disp = v->pt - v1->pt;
		double transverse_distance = fabs(disp.x*(-line_dir.y) + disp.y*line_dir.x);
		collinear_v = (transverse_distance < dist_eps && (fabs(alpha_v - e2->get_alpha_1()) < alpha_eps || fabs(alpha_v - e2->get_alpha_2()) < alpha_eps));

		collinear_1 = collinear_2 = collinear_v;

		if (collinear_v)
		{
			// Finds the edge we have not explored so far on v1's side
			edges_to_delete.push_back(e1);
			Vertex* v1_next = v1;
			while (v1->connectivity() == 2 && v1->id_vertex > 3 && collinear_1) {
				if (v1->directions[0].second->e == e1) {
					e1 = v1->directions[1].second->e;
				}
				else {
					e1 = v1->directions[0].second->e;
				}
				v1_next = (e1->v1 == v1 ? e1->v2 : e1->v1);
				// test if v1 and v is collinear with v1_next--v2
				line_dir = v1_next->pt - v2->pt;
				line_dir /= cv::norm(line_dir);
				disp = v->pt - v1_next->pt;
				transverse_distance = fabs(disp.x*(-line_dir.y) + disp.y*line_dir.x);
				collinear_1 = (transverse_distance < 1e-4 && (fabs(alpha_v - e1->get_alpha_1()) < alpha_eps || fabs(alpha_v - e1->get_alpha_2()) < alpha_eps));

				if (collinear_1) {
					vertices_erased.push_back(v1);
					edges_to_delete.push_back(e1);
					if (e_type == INNER_EDGE) {
						ie = static_cast<Inner_Edge*>(e1);
						ie->get_supporting_rays(rays);
					}
					v1 = v1_next;
				}
			}
			e1 = edges_to_delete.back();

			// Same operation for v2
			edges_to_delete.push_back(e2);
			Vertex* v2_next = v2;
			while (v2->connectivity() == 2 && v2->id_vertex > 3 && collinear_2) {
				if (v2->directions[0].second->e == e2) {
					e2 = v2->directions[1].second->e;
				}
				else {
					e2 = v2->directions[0].second->e;
				}
				v2_next = (e2->v1 == v2 ? e2->v2 : e2->v1);
				// test if v2 and v is collinear with v1--v2_next
				line_dir = v2_next->pt - v1->pt;
				line_dir /= cv::norm(line_dir);
				transverse_distance = fabs(disp.x*(-line_dir.y) + disp.y*line_dir.x);
				collinear_2 = (transverse_distance < dist_eps && (fabs(alpha_v - e2->get_alpha_1()) < alpha_eps || fabs(alpha_v - e2->get_alpha_2()) < alpha_eps));

				if (collinear_2) {
					vertices_erased.push_back(v2);
					edges_to_delete.push_back(e2);
					if (e_type == INNER_EDGE) {
						ie = static_cast<Inner_Edge*>(e2);
						ie->get_supporting_rays(rays);
					}
					v2 = v2_next;
				}
			}
			e2 = edges_to_delete.back();

			// Creates the merged edge
			for (auto dir : v1->directions) {
				if (dir.second->e->v1 == v2 || dir.second->e->v2 == v2) return e_12;
			}
			/* Get the smallest width among collinear edges */
			double width = FLT_MAX;
			for (auto it = edges_to_delete.begin(); it != edges_to_delete.end(); ++it) { 
				if ((*it)->width < width) width = (*it)->width;
			}
			double alpha_v1 = atan2(v2->pt.y - v1->pt.y, v2->pt.x - v1->pt.x);
			double alpha_v2 = (alpha_v1 <= 0 ? alpha_v1 + PI : alpha_v1 - PI);
			if (e_type == INNER_EDGE) {
				e_12 = new Inner_Edge(id_edges, rays, 0, v1, v2, width, alpha_v1, alpha_v2);
			}
			else if (e_type == OUTER_EDGE) {
				e_12 = new Outer_Edge(id_edges, boundary, v1, v2, width, alpha_v1, alpha_v2, v1_2);
			}
			else {
				double ae_a = -sin(alpha_v1), ae_b = cos(alpha_v1), ae_c = -ae_a * v1->pt.y - ae_b * v1->pt.y;
				e_12 = new Artificial_Edge(id_edges, ae_a, ae_b, ae_c, v1, v2, width, alpha_v1, alpha_v2);
			}
			/*int n_pixels = 0;
			for (list<Edge*>::iterator it_e = edges_to_delete.begin(); it_e != edges_to_delete.end(); ++it_e) n_pixels += (*it_e)->region.size();
			e_12->region.reserve(n_pixels);
			for (list<Edge*>::iterator it_e = edges_to_delete.begin(); it_e != edges_to_delete.end(); ++it_e) {
				for (auto it_p = (*it_e)->region.begin(); it_p != (*it_e)->region.end(); ++it_p)
					e_12->region.push_back(std::move(*it_p));
			}*/

			if (f_v1_v2 != nullptr) {
				f_v1_v2->remove_non_corners(v1, e1, v2, e_12->v1_v2);
			}
			if (f_v2_v1 != nullptr) {
				f_v2_v1->remove_non_corners(v2, e2, v1, e_12->v2_v1);			
			}

			// Deletes the references to the splitted edges and the bivalent edges
			erase(edges_to_delete, true);
			for (list<Vertex *>::iterator it_v = vertices_erased.begin(); it_v != vertices_erased.end(); it_v++) {
				queue_bivalent.erase(pair<Vertex*, double>(*it_v, (*it_v)->energy_gain));
				bboxes_by_x.erase(pair<Vertex*, double>(*it_v, (*it_v)->bbox.first));
				bboxes_by_y.erase(pair<Vertex*, double>(*it_v, (*it_v)->bbox.second));
				erase(*it_v);
				quadtree->remove((*it_v), destroy);
			}

			// Deletes the last bivalent edge
			queue_bivalent.erase(pair<Vertex*, double>(v,v->energy_gain));
			bboxes_by_x.erase(pair<Vertex*, double>(v, v->bbox.first));
			bboxes_by_y.erase(pair<Vertex*, double>(v, v->bbox.second));
			erase(v);
			quadtree->remove(v, destroy);

			// Adds the merged edge to the list of vertices
			push_back(e_12);
		}
	}
	return e_12;
}

void Partition::seal_and_remove_bivalent_vertices(int n_boundary_vertices)
{
	for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
		Face* f = (*it_f);
		for (list<HalfEdge *>::iterator it_h = f->edges[0].begin() ; it_h != f->edges[0].end() ; it_h++) {
			HalfEdge* h = (*it_h);
			int i = h->e->id_edge;
		}
	}

	Vertex* v = vertices_head;
	while (v != NULL) {
		assert(v->connectivity() > 1);

		// The first few vertices are image corners
		if (v->connectivity() == 2 && v->id_vertex > n_boundary_vertices-1) {

			list<Edge *> edges_to_delete;
			list<Vertex *> vertices_to_delete;
			set<Ray *> rays;
			Image_Boundary boundary = INVALID_BORDER;
			pair<int, int> init_boundary_location;
			uint support_boundary_index;
			Face* outer_face = nullptr;
			bool v1_2 = true;

			// The vertex v is bivalent : that's why we find the vertices to which v is connected
			// However, maybe the connected vertices are themselves bivalent : we need to iterate on v1 and v2 until
			// finding trivalent vertices. Until then, we define vertices and edges to delete
			Edge* e1 = v->directions[0].second->e;
			Edge* e2 = v->directions[1].second->e;

			Edge_Type e_type = e1->type;
			Inner_Edge *ie = nullptr;
			Outer_Edge* oe = nullptr;

			Vertex* v1 = (e1->v1 == v ? e1->v2 : e1->v1);
			Vertex* v2 = (e2->v1 == v ? e2->v2 : e2->v1);
			
			// We define a sequence of vertices (v1 .. v .. v2)
			// Obtains the two facets adjacent to the big edge (v1 v2) : to this end, we can use e2
			Face *f_v1_v2 = nullptr, *f_v2_v1 = nullptr;
			HalfEdge* h = e2->v1_v2;
			if (h->e->v1 == v && h->e->v2 == v2) {
				f_v1_v2 = h->f;
				f_v2_v1 = h->opposite()->f;
			} else {
				f_v1_v2 = h->opposite()->f;
				f_v2_v1 = h->f;
			}

			if (e_type == INNER_EDGE) {
				ie = static_cast<Inner_Edge*>(e1);
				ie->get_supporting_rays(rays);
				ie = static_cast<Inner_Edge*>(e2);
				ie->get_supporting_rays(rays);
			} else {
				oe = static_cast<Outer_Edge*>(e1);
				boundary = oe->boundary;
				init_boundary_location = oe->init_boundary_location;
				support_boundary_index = oe->support_boundary_index;
				outer_face = oe->outer_face;
				v1_2 = oe->v1_2;
			}

			
			// Finds the edge we have not explored so far on v1's side
			edges_to_delete.push_back(e1);
			while (v1->connectivity() == 2 && v1->id_vertex > n_boundary_vertices - 1) {
				vertices_to_delete.push_back(v1);
				if (v1->directions[0].second->e == e1) {
					e1 = v1->directions[1].second->e;
				} else {
					e1 = v1->directions[0].second->e;
				}
				edges_to_delete.push_back(e1);
				if (e_type == INNER_EDGE) {
					ie = static_cast<Inner_Edge*>(e1);
					ie->get_supporting_rays(rays);
				}
				v1 = (e1->v1 == v1 ? e1->v2 : e1->v1);
			}
			e1 = edges_to_delete.back();

			// Same operation for v2
			edges_to_delete.push_back(e2);
			while (v2->connectivity() == 2 && v2->id_vertex > n_boundary_vertices - 1) {
				vertices_to_delete.push_back(v2);
				if (v2->directions[0].second->e == e2) {
					e2 = v2->directions[1].second->e;
				} else {
					e2 = v2->directions[0].second->e;
				}
				edges_to_delete.push_back(e2);

				if (e_type == INNER_EDGE) {
					ie = static_cast<Inner_Edge*>(e2);
					ie->get_supporting_rays(rays);
				}
				v2 = (e2->v1 == v2 ? e2->v2 : e2->v1);
			}
			e2 = edges_to_delete.back();

			// Creates the merged edge
			Edge* e_12 = nullptr;
			/* Get the smallest width among collinear edges */
			double width = FLT_MAX;
			for (auto it = edges_to_delete.begin(); it != edges_to_delete.end(); ++it) {
				if ((*it)->width < width) width = (*it)->width;
			}
			double alpha_v1 = atan2(v2->pt.y - v1->pt.y, v2->pt.x - v1->pt.x);
			double alpha_v2 = (alpha_v1 <= 0 ? alpha_v1 + PI : alpha_v1 - PI);
			if (e_type == INNER_EDGE) {
				e_12 = new Inner_Edge(id_edges, rays, 0, v1, v2, width, alpha_v1, alpha_v2);
			} else {
				e_12 = new Outer_Edge(id_edges, boundary, outer_face, init_boundary_location, support_boundary_index, 
					v1, v2, width, alpha_v1, alpha_v2, v1_2, 0, 0);
			}

			if (f_v1_v2 != nullptr) {
				f_v1_v2->remove_non_corners(v1, e1, v2, e_12->v1_v2);

				for (list<HalfEdge *>::iterator it_h = f_v1_v2->edges[0].begin() ; it_h != f_v1_v2->edges[0].end() ; it_h++) {
					int i = (*it_h)->e->id_edge;
				}
			}

			if (f_v2_v1 != nullptr) {
				f_v2_v1->remove_non_corners(v2, e2, v1, e_12->v2_v1);

				for (list<HalfEdge *>::iterator it_h = f_v2_v1->edges[0].begin() ; it_h != f_v2_v1->edges[0].end() ; it_h++) {
					int i = (*it_h)->e->id_edge;
				}
			}

			// Deletes the references to the splitted edges and the bivalent edges
			erase(edges_to_delete, true);
			for (list<Vertex *>::iterator it_v = vertices_to_delete.begin(); it_v != vertices_to_delete.end(); it_v++) {
				erase(*it_v);
				quadtree->remove((*it_v), true);
			}

			// Deletes the last bivalent edge
			Vertex* v_next = erase(v);
			quadtree->remove(v, true);
			v = v_next;

			// Adds the merged edge to the list of vertices
			push_back(e_12);

		} else {
			v = v->v_next;
		}
	}

	for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
		Face* f = (*it_f);
		for (list<HalfEdge *>::iterator it_h = f->edges[0].begin() ; it_h != f->edges[0].end() ; it_h++) {
			HalfEdge* h = (*it_h);
			int i = h->e->id_edge;
		}
	}
}


bool Partition::check_removal_criteria(Vertex* v, Matrix<double> & I_m, const Matrix<double> & I_t, Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used)
{
	bool valid = true;
	if (v->connectivity() != 2 || (v->id_vertex >= 0 && v->id_vertex < 4) ) {
		valid = false;	
	}

	HalfEdge* h1 = v->directions[0].second;
	HalfEdge* h2 = v->directions[1].second;
	Vertex *v1 = h1->e->v1 != v ? h1->e->v1 : h1->e->v2;
	Vertex *v2 = h2->e->v1 != v ? h2->e->v1 : h2->e->v2;
	// v is not considered a valid bivalent vertex to remove if removing v results in deleting a face
	for (const auto& dir: v1->directions) {
		Vertex* v_dir = dir.second->e->v1 != v1 ? dir.second->e->v1 : dir.second->e->v2;
		if (v_dir == v2) {
			valid = false;
			break;
		}
	}

	if (!valid) {
		queue_bivalent.erase(pair<Vertex*, double>(v,v->energy_gain));
		bboxes_by_x.erase(pair<Vertex*, double>(v, v->bbox.first));
		bboxes_by_y.erase(pair<Vertex*, double>(v, v->bbox.second));
		v->energy_gain = -10000;
		return false;
	}
	
	double alpha_v1v = h1->e->v1 != v ? h1->e->get_alpha_1() : h1->e->get_alpha_2();
	double alpha_v2v = h2->e->v1 != v ? h2->e->get_alpha_1() : h2->e->get_alpha_2();

	// set a bounding box
	double x_min, x_max, y_min, y_max;
	x_min = min(min(v1->pt.x, v->pt.x), v2->pt.x);
	x_max = max(max(v1->pt.x, v->pt.x), v2->pt.x);
	y_min = min(min(v1->pt.y, v->pt.y), v2->pt.y);
	y_max = max(max(v1->pt.y, v->pt.y), v2->pt.y);
	bboxes_by_x.erase(pair<Vertex*, double>(v, v->bbox.first));
	bboxes_by_y.erase(pair<Vertex*, double>(v, v->bbox.second));
	v->bbox = pair<double, double>(x_max - x_min, y_max - y_min);
	bboxes_by_x.emplace(pair<Vertex*, double>(v, v->bbox.first));
	bboxes_by_y.emplace(pair<Vertex*, double>(v, v->bbox.second));

	
	list<Point2d*> triangle;
	triangle.push_back(&(v1->pt));
	triangle.push_back(&(v->pt));
	triangle.push_back(&(v2->pt));

	// check if neighbor vertices fall in the triangle
	std::list<Vertex *> neighbors;
	quadtree->search(x_min, x_max, y_min, y_max, neighbors);
	for (Vertex* v_neighbor : neighbors) {
		if (v_neighbor == v || v_neighbor == v1 || v_neighbor == v2) continue;
		if (Geometry::is_inside_contour(v_neighbor->pt, triangle)) {
			queue_bivalent.erase(pair<Vertex*, double>(v, v->energy_gain));
			v->energy_gain = -10000;
			return false;
		}
	}

	double fidelity_p_before(0), fidelity_e_before(0), fidelity_p_after(0), fidelity_e_after(0);
	double alpha_1 = atan2(v2->pt.y - v1->pt.y, v2->pt.x - v1->pt.x);
	double alpha_2 = alpha_1 <= 0 ? alpha_1 + PI : alpha_1 - PI;
	for (auto dir : v1->directions) {
		if (dir.second->e == h1->e) continue;
		if (Geometry::is_between(dir.first, alpha_1, alpha_v1v) || std::abs(dir.first - alpha_1) < 1e-4 || std::abs(std::abs(dir.first - alpha_1) - 2 * PI) < 1e-4) {
			queue_bivalent.erase(pair<Vertex*, double>(v, v->energy_gain));
			v->energy_gain = -10000;
			return false;
		}
	}
	for (auto dir : v2->directions) {
		if (dir.second->e == h2->e) continue;
		if (Geometry::is_between(dir.first, alpha_2, alpha_v2v) || std::abs(dir.first - alpha_2) < 1e-4 || std::abs(std::abs(dir.first - alpha_2) - 2 * PI) < 1e-4) {
			queue_bivalent.erase(pair<Vertex*, double>(v, v->energy_gain));
			v->energy_gain = -10000;
			return false;
		}
	}

	const Face* f_a = h1->f, *f_b = h1->opposite()->f;
	bool consistent_label = true;
	if (f_a != nullptr && f_b != nullptr) consistent_label = (f_a->label == f_b->label);

	int id = -3;
	Vertex v_copy(id, v->pt.x, v->pt.y), v1_copy(id, v1->pt.x, v1->pt.y), v2_copy(id, v2->pt.x, v2->pt.y);
	id = -1;
	Edge* e_v1v2 = new Inner_Edge(id, set<Ray*>(), 0, &v1_copy, &v2_copy, 1.0, alpha_1, alpha_2);
	e_v1v2->weight(used, I_m, I_t, offset, params->split_params, cdf_grad_m);

	if (consistent_label) {
		double fidelity_e_before = (fidelity_edge(consistent_label, h1->e->unit_grad_cost(), h1->e->length) + 
									fidelity_edge(consistent_label, h2->e->unit_grad_cost(), h2->e->length))/ diag_length;
		double fidelity_e_after = fidelity_edge(consistent_label, e_v1v2->unit_grad_cost(), e_v1v2->length) / diag_length;
		queue_bivalent.erase(pair<Vertex*, double>(v, v->energy_gain));
		v->energy_gain = params->fidelity_beta*(fidelity_e_before - fidelity_e_after) + params->prior_lambda;
		queue_bivalent.emplace(v, v->energy_gain);
		delete e_v1v2;
		return true;
	}

	Edge* e_vv1, *e_v2v;
	// copy existing edges
	double alpha_vv1 = h1->e->v1 == v ? h1->e->get_alpha_1() : h1->e->get_alpha_2();
	double alpha_vv2 = h2->e->v1 == v ? h2->e->get_alpha_1() : h2->e->get_alpha_2();
	id = -1;
	e_vv1 = new Inner_Edge(id, set<Ray*>(), 0, &v_copy, &v1_copy, 1.0, alpha_vv1, alpha_v1v);
	id = -1;
	e_v2v = new Inner_Edge(id, set<Ray*>(), 0, &v2_copy, &v_copy, 1.0, alpha_v2v, alpha_vv2);
	// triangle: v->v1->v2
	list<HalfEdge*> triangle_edges = { e_vv1->v1_v2, e_v1v2->v1_v2, e_v2v->v1_v2 };
	Label_ID label_before = f_a->label, label_after = f_b->label;
	if (HalfEdge::is_clockwise(triangle_edges) != 1) {
		list<HalfEdge*> triangle_copy = triangle_edges;
		triangle_edges.clear();
		for (HalfEdge* h: triangle_copy) triangle_edges.push_front(h->opposite());
		Label_ID temp = label_before;
		label_before = label_after;
		label_after = temp;
	}
	id = -1;
	Face* triangle_face = new Face(id, triangle_edges);
	triangle_face->find_pixels_inside_facet(offset);
	if (triangle_face->pixels.size() >= f_a->pixels.size() || triangle_face->pixels.size() >= f_b->pixels.size()) {
		delete triangle_face;
		delete e_v1v2;
		delete e_vv1;
		delete e_v2v;
		queue_bivalent.erase(pair<Vertex*, double>(v, v->energy_gain));
		v->energy_gain = -10000;
		return false;
	}
	triangle_face->compute_feature(I_prob, I, offset, params->max_label);
	fidelity_p_before += triangle_face->pixels.size() - triangle_face->semantic_probabilities[label_before];
	fidelity_p_before /= double(area);
	fidelity_p_after += triangle_face->pixels.size() - triangle_face->semantic_probabilities[label_after];
	fidelity_p_after /= double(area);

	fidelity_e_before += fidelity_edge(consistent_label, h1->e->unit_grad_cost(), h1->e->length);
	fidelity_e_before += fidelity_edge(consistent_label, h2->e->unit_grad_cost(), h2->e->length);
	fidelity_e_before /= diag_length;
	fidelity_e_after = fidelity_edge(consistent_label, e_v1v2->unit_grad_cost(), e_v1v2->length) / diag_length;
	double fidelity_before = (1 - params->fidelity_beta)*fidelity_p_before + params->fidelity_beta*fidelity_e_before;
	double fidelity_after = (1 - params->fidelity_beta)*fidelity_p_after + params->fidelity_beta*fidelity_e_after;

	queue_bivalent.erase(pair<Vertex*, double>(v, v->energy_gain));
	v->energy_gain = fidelity_before - fidelity_after + params->prior_lambda;
	queue_bivalent.emplace(v, v->energy_gain);


	// free memory
	delete triangle_face;
	delete e_v1v2;
	delete e_vv1; 
	delete e_v2v;

	return (v->energy_gain>0);
}

void Partition::remove_bivalent_vertiex(Vertex* v, Matrix<double> & I_m, Matrix<double> & I_t, Matrix<double> & I_lab_m, Matrix<double> & I_lab_t, 
	Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, map<int, double> & angles, vector<int> & bin_to_angle_index)
{
	Edge* e1 = v->directions[0].second->e;
	Edge* e2 = v->directions[1].second->e;
	Vertex* v1 = (e1->v1 == v ? e1->v2 : e1->v1);
	Vertex* v2 = (e2->v1 == v ? e2->v2 : e2->v1);

	if (!check_removal_criteria(v, I_m, I_t, I, I_prob, used)) {
		trace(params->verbose_level, 5, "Cannot remove vertex " + std::to_string(v->id_vertex) + " between vertices " +
			std::to_string(v1->id_vertex) + " and " + std::to_string(v2->id_vertex));
		return;
	}
	
	Face *f_v1_v2 = nullptr, *f_v2_v1 = nullptr;
	HalfEdge* h = e2->v1_v2;
	if (h->e->v1 == v) {
		f_v1_v2 = h->f;
		f_v2_v1 = h->opposite()->f;
	}
	else {
		f_v1_v2 = h->opposite()->f;
		f_v2_v1 = h->f;
	}
	

	for (auto dir : v1->directions) {
		Vertex* v_dir = dir.second->e->v1 != v1 ? dir.second->e->v1 : dir.second->e->v2;
		if (v_dir == v2) {
			if (f_v1_v2 != nullptr && f_v2_v1 != nullptr) {
				vector<Dual_Graph::vertex_descriptor> vds;
				vds.push_back(face_to_graph.at(f_v1_v2));
				vds.push_back(face_to_graph.at(f_v2_v1));
				merge_operator(I, I_prob, used, I_m, I_t, I_lab_m, I_lab_t, vds, angles, bin_to_angle_index);
			}
			else {
				assert(false);
			}
			return;
		}
	}

	const vector<Point2i> & e1_region = e1->region;
	for (auto it_p = e1_region.begin(); it_p != e1_region.end(); ++it_p) {
		int i = rows - 1 - it_p->y; // rows
		int j = it_p->x; // columns
		used(i, j) = false;
	}
	const vector<Point2i> & e2_region = e2->region;
	for (auto it_p = e2_region.begin(); it_p != e2_region.end(); ++it_p) {
		int i = rows - 1 - it_p->y; // rows
		int j = it_p->x; // columns
		used(i, j) = false;
	}

	int id_v1_v2 = f_v1_v2 != nullptr ? f_v1_v2->id_face : -1;
	int id_v2_v1 = f_v2_v1 != nullptr ? f_v2_v1->id_face : -1;
	if (f_v1_v2 != nullptr && f_v2_v1 != nullptr) {
		trace(params->verbose_level, 5, "remove bivalent vertex " + std::to_string((*queue_bivalent.begin()).first->id_vertex) + " between vertices " +
			std::to_string(v1->id_vertex) + " and " + std::to_string(v2->id_vertex) + " from face " + std::to_string(id_v1_v2) + " (lab=" + std::to_string(f_v1_v2->label)
			+ ") and face " + std::to_string(id_v2_v1) + " (lab=" + std::to_string(f_v2_v1->label) + ")");
	}
	else {
		trace(params->verbose_level, 5, "remove bivalent vertex " + std::to_string((*queue_bivalent.begin()).first->id_vertex) + " between vertices " +
			std::to_string(v1->id_vertex) + " and " + std::to_string(v2->id_vertex) );
	}
	

	set<Vertex*> vertices_to_check;
	vertices_to_check.emplace(v1);
	vertices_to_check.emplace(v2);
	list<Vertex*> neighbors;
	double x_min = v->pt.x - bboxes_by_x.begin()->second;
	double x_max = v->pt.x + bboxes_by_x.begin()->second;
	double y_min = v->pt.y - bboxes_by_y.begin()->second;
	double y_max = v->pt.y + bboxes_by_y.begin()->second;
	quadtree->search(x_min, x_max, y_min, y_max, neighbors);
	vertices_to_check.insert(neighbors.begin(), neighbors.end());
	vertices_to_check.erase(v);
	
	double width = min(e1->width, e2->width);
	double alpha_1 = atan2(v2->pt.y - v1->pt.y, v2->pt.x - v1->pt.x);
	double alpha_2 = (alpha_1 <= 0 ? alpha_1 + PI : alpha_1 - PI);
	double alpha_v1v = e1->v1 != v ? e1->get_alpha_1() : e1->get_alpha_2();
	double alpha_v2v = e2->v1 != v ? e2->get_alpha_1() : e2->get_alpha_2();
	bool resplit = true;
	if (std::abs(alpha_v1v - alpha_1) < 0.09 || std::abs(std::abs(alpha_v1v - alpha_1) - 2 * PI) < 0.09) {
		if (std::abs(alpha_v2v - alpha_2) < 0.09 || std::abs(std::abs(alpha_v2v - alpha_2) - 2 * PI) < 0.09) {
			resplit = false;
		}
	}

	Edge* e_v1v2;
	if (e1->type == OUTER_EDGE) {
		Outer_Edge* oe = static_cast<Outer_Edge*>(e1);
		e_v1v2 = new Outer_Edge(id_edges, oe->boundary, v1, v2, width, alpha_1, alpha_2, oe->v1_2);
	}
	else {
		e_v1v2 = new Inner_Edge(id_edges, set<Ray*>(), 0, v1, v2, width, alpha_1, alpha_2);
	}

	e_v1v2->weight(used, I_m, I_t, offset, params->split_params, cdf_grad_m);

	const vector<Point2i> & ev1v2_region = e_v1v2->region;
	for (auto it_p = ev1v2_region.begin(); it_p != ev1v2_region.end(); ++it_p) {
		int i = rows - 1 - it_p->y; // rows
		int j = it_p->x; // columns
		used(i, j) = true;
	}

	if (params->verbose_level >= 6 && f_v1_v2 != nullptr && f_v2_v1 != nullptr) {
		Label_ID label_before = v->directions[0].second->f->label;
		Label_ID label_after = v->directions[0].second->opposite()->f->label;
		// triangle: v->v1->v2
		list<HalfEdge*> triangle_edges = { v->directions[0].second, e_v1v2->v1_v2, v->directions[0].second->opposite() }; 
		if (HalfEdge::is_clockwise(triangle_edges) != 1) {
			Label_ID temp = label_before;
			label_before = label_after;
			label_after = temp;
		}
		trace(params->verbose_level, 6, "triangle label before = " + std::to_string(label_before) +
			", label after = " + std::to_string(label_after));
	}
	
	vector<Segment*> segments_v1_v2, segments_v2_v1;
	size_t n_v1_v2_old(0), n_v2_v1_old(0);
	if (f_v1_v2 != nullptr) {
		n_v1_v2_old = f_v1_v2->pixels.size();
		if (resplit) {
			segments_v1_v2.reserve(f_v1_v2->detected_segments.size());
			std::copy(f_v1_v2->detected_segments.begin(), f_v1_v2->detected_segments.end(), segments_v1_v2.begin());
			f_v1_v2->split_vertex_chains.clear();
			f_v1_v2->split_edges_property.clear();
			f_v1_v2->split_vertex_outer_edge_location.clear();
		}
		f_v1_v2->remove_non_corners(v1, e1, v2, e_v1v2->v1_v2);
	}
	
	if (f_v2_v1 != nullptr) {
		n_v2_v1_old = f_v2_v1->pixels.size();
		if (resplit) {
			segments_v2_v1.reserve(f_v2_v1->detected_segments.size());
			std::copy(f_v2_v1->detected_segments.begin(), f_v2_v1->detected_segments.end(), segments_v2_v1.begin());
			f_v2_v1->split_vertex_chains.clear();
			f_v2_v1->split_edges_property.clear();
			f_v2_v1->split_vertex_outer_edge_location.clear();
		}

		f_v2_v1->remove_non_corners(v2, e2, v1, e_v1v2->v2_v1);
	}
	;
	erase(e1, true);
	erase(e2, true);
	push_back(e_v1v2);

	

	if (f_v1_v2 != nullptr) {
		f_v1_v2->find_pixels_inside_facet(offset);

		f_v1_v2->find_active_gradients(I_m, I_t, offset, params->split_params, true); // call after weighting edges
		f_v1_v2->find_active_gradients(I_lab_m, I_lab_t, offset, params->split_params, false);
		f_v1_v2->compute_feature(I_prob, I, offset, params->max_label);

		if (n_v1_v2_old < f_v1_v2->pixels.size() && resplit) {
			for (size_t i = 0; i < segments_v2_v1.size(); i++) {
				Point2d & s_i_O = segments_v2_v1[i]->finalBarycenter;
				bool is_inside_polygon = Geometry::is_inside_contour(s_i_O, f_v1_v2->edges[0]);
				if (is_inside_polygon) {
					for (int id = 1; id < f_v1_v2->edges.size(); ++id) {
						if (Geometry::is_inside_contour(s_i_O, f_v1_v2->edges[id])) {
							is_inside_polygon = false;
							break;
						}
					}
				}
				if (is_inside_polygon) {
					f_v1_v2->detected_segments.insert(segments_v2_v1[i]);
					if (f_v2_v1 != nullptr) f_v2_v1->detected_segments.erase(segments_v2_v1[i]);
				}
			}
		}
	}
	if (f_v2_v1 != nullptr) {
		f_v2_v1->find_pixels_inside_facet(offset);

		f_v2_v1->find_active_gradients(I_m, I_t, offset, params->split_params, true);
		f_v2_v1->find_active_gradients(I_lab_m, I_lab_t, offset, params->split_params, false);
		f_v2_v1->compute_feature(I_prob, I, offset, params->max_label);

		if (n_v2_v1_old < f_v2_v1->pixels.size() && resplit) {
			for (size_t i = 0; i < segments_v1_v2.size(); i++) {
				Point2d & s_i_O = segments_v1_v2[i]->finalBarycenter;
				bool is_inside_polygon = Geometry::is_inside_contour(s_i_O, f_v2_v1->edges[0]);
				if (is_inside_polygon) {
					for (int id = 1; id < f_v2_v1->edges.size(); ++id) {
						if (Geometry::is_inside_contour(s_i_O, f_v2_v1->edges[id])) {
							is_inside_polygon = false;
							break;
						}
					}
				}
				if (is_inside_polygon) {
					f_v2_v1->detected_segments.insert(segments_v1_v2[i]);
					if (f_v1_v2 != nullptr) f_v1_v2->detected_segments.erase(segments_v1_v2[i]);
				}
			}
		}
	}


	// update bivalent queue
	queue_bivalent.erase(pair<Vertex*, double>(v, v->energy_gain));
	bboxes_by_x.erase(pair<Vertex*, double>(v, v->bbox.first));
	bboxes_by_y.erase(pair<Vertex*, double>(v, v->bbox.second));
	erase(v);
	quadtree->remove(v, true);
	for (Vertex* v_check : vertices_to_check) check_removal_criteria(v_check, I_m, I_t, I, I_prob, used);

	// update graph
	if (f_v1_v2 != nullptr) {
		Dual_Graph::vertex_descriptor v_12 = face_to_graph.at(f_v1_v2);
		// remove old vertex from the queue
		queue_split.erase(make_pair(v_12, dual_graph[v_12].t));
		// add cost on the vertices here
		extend_detected_segments(f_v1_v2, I, I_prob, used, I_m, I_t, angles, bin_to_angle_index, dual_graph[v_12].t);
		f_v1_v2->delete_invalid_segments();
		queue_split.insert(make_pair(v_12, dual_graph[v_12].t));

		Dual_Graph::adjacency_iterator v_12_it, v_12_end;
		boost::tie(v_12_it, v_12_end) = adjacent_vertices(v_12, dual_graph);

		while (v_12_it != v_12_end) {
			Dual_Graph::vertex_descriptor v_adj = *v_12_it;
			bool found;
			Dual_Graph::edge_descriptor connect_e;
			boost::tie(connect_e, found) = boost::edge(v_adj, v_12, dual_graph);
			queue_merge.erase(make_pair(connect_e, dual_graph[connect_e].t));

			merge_gain(f_v1_v2, dual_graph[v_adj].f, I, I_m, I_t, dual_graph[connect_e].t, dual_graph[connect_e].merged_label);
			queue_merge.insert(make_pair(connect_e, dual_graph[connect_e].t));
			
			++v_12_it;
		}
	}
	if (f_v2_v1 != nullptr) {
		Dual_Graph::vertex_descriptor v_21 = face_to_graph.at(f_v2_v1);
		queue_split.erase(make_pair(v_21, dual_graph[v_21].t));
		extend_detected_segments(f_v2_v1, I, I_prob, used, I_m, I_t, angles, bin_to_angle_index, dual_graph[v_21].t);
		f_v2_v1->delete_invalid_segments();
		queue_split.insert(make_pair(v_21, dual_graph[v_21].t));

		Dual_Graph::adjacency_iterator v_21_it, v_21_end;
		boost::tie(v_21_it, v_21_end) = adjacent_vertices(v_21, dual_graph);

		while (v_21_it != v_21_end) {
			Dual_Graph::vertex_descriptor v_adj = *v_21_it;
			Face* f_adj = dual_graph[v_adj].f;
			if (f_adj != f_v1_v2) {
				bool found;
				Dual_Graph::edge_descriptor connect_e;
				boost::tie(connect_e, found) = boost::edge(v_adj, v_21, dual_graph);
				queue_merge.erase(make_pair(connect_e, dual_graph[connect_e].t));

				merge_gain(f_v2_v1, dual_graph[v_adj].f, I, I_m, I_t, dual_graph[connect_e].t, dual_graph[connect_e].merged_label);
				queue_merge.insert(make_pair(connect_e, dual_graph[connect_e].t));
				
			}
			++v_21_it;
		}
	}
}

void Partition::build_merge_graph(Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t)
{
	merge_graph.clear();
	queue_merge_graph.clear();
	map<Face *, Dual_Merge_Graph::vertex_descriptor> f_to_vid;
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); ++it_f) {
		Dual_Merge_Graph::vertex_descriptor  u = boost::add_vertex(merge_graph);
		merge_graph[u].superfacet.push_back(*it_f);
		merge_graph[u].label = (*it_f)->label;
		merge_graph[u].semantic_probabilities = (*it_f)->semantic_probabilities;
		merge_graph[u].n_pixels = (*it_f)->pixels.size();
		size_t n_edges = 0;
		double perimeter = 0;
		for (int id = 0; id < (*it_f)->edges.size(); id++) {
			n_edges += (*it_f)->edges[id].size();
			for (list<HalfEdge *>::iterator it_h = (*it_f)->edges[id].begin(); it_h != (*it_f)->edges[id].end(); ++it_h) {
				if ((*it_h)->e->type == OUTER_EDGE) {
					Outer_Edge* oe = static_cast<Outer_Edge*>((*it_h)->e);
					Face* outer_face = oe->outer_face;
					Label_ID boundary_label = outer_face == nullptr ? INVALID_LABEL : outer_face->label;
					merge_graph[u].boundary_edges.push_back({ (*it_h)->e, boundary_label });
				}
				perimeter += (*it_h)->e->length;
			}
		}
		merge_graph[u].perimeter = perimeter;
		merge_graph[u].n_edges = n_edges;
		f_to_vid[*it_f] = u;
	}

	Edge* e = edges_head;
	while (e != NULL) {
		if (e->v1_v2->f != nullptr && e->v2_v1->f != nullptr) {
			Dual_Merge_Graph::vertex_descriptor  u = f_to_vid.at(e->v1_v2->f);
			Dual_Merge_Graph::vertex_descriptor  v = f_to_vid.at(e->v2_v1->f);
			Dual_Merge_Graph::edge_descriptor connect_e;
			bool found;
			boost::tie(connect_e, found) = boost::edge(u, v, merge_graph);
			if (!found) {
				bool inserted;
				boost::tie(connect_e, inserted) = boost::add_edge(u, v, merge_graph);
				merge_graph[connect_e].shared_edges.push_back(e);
			}
			else {
				list<Edge *> & shared = merge_graph[connect_e].shared_edges;
				merge_graph[connect_e].shared_edges.push_back(e);
			}
		}
		e = e->e_next;
	}

	Dual_Merge_Graph::edge_iterator it, it_end;
	boost::tie(it, it_end) = boost::edges(merge_graph);
	while (it != it_end) {
		Dual_Merge_Graph::vertex_descriptor  u = boost::source(*it, merge_graph), v = boost::target(*it, merge_graph);
		double &t = merge_graph[*it].t;
		Label_ID & merged_label = merge_graph[*it].merged_label;
		map<list<Edge*>, bool> shared_edges;
		shared_edges[merge_graph[*it].shared_edges] = (merge_graph[u].label==merge_graph[v].label);
		merge_gain({u, v}, shared_edges, I, I_m, I_t, t, merged_label, false);
		queue_merge_graph.insert(make_pair(*it, t));
		++it;
	}
	assert(boost::num_vertices(merge_graph) == faces.size());
}


void Partition::build_dual_graph(Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t)
{	
	queue_split.clear();
	queue_merge.clear();
	dual_graph.clear();
	face_to_graph.clear();

	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); ++it_f) {
		Dual_Graph::vertex_descriptor  u = boost::add_vertex(dual_graph);
		dual_graph[u].f = *it_f;
		face_to_graph[*it_f] = u;
		dual_graph[u].t = (*it_f)->split_gain;

		auto inserted = queue_split.insert(make_pair(u, dual_graph[u].t));
	}

	Edge* e = edges_head;
	while (e != NULL) {
		if ((e->v1_v2->f != nullptr) && (e->v2_v1->f != nullptr)) {
			assert(e->v1_v2->f != e->v2_v1->f);
			Face* f_u = e->v1_v2->f, *f_v = e->v2_v1->f;
			Dual_Graph::vertex_descriptor  u = face_to_graph.at(f_u);
			Dual_Graph::vertex_descriptor  v = face_to_graph.at(f_v);
			Dual_Graph::edge_descriptor connect_e;
			bool found;
			boost::tie(connect_e, found) = boost::edge(u, v, dual_graph);
			if (!found) {
				bool inserted;
				boost::tie(connect_e, inserted) = (f_u->pixels.size() > f_v->pixels.size()) ? boost::add_edge(u, v, dual_graph) : boost::add_edge(v, u, dual_graph);
				double &t = dual_graph[connect_e].t;
				Label_ID &merged_label = dual_graph[connect_e].merged_label;
				merge_gain(e->v1_v2->f, e->v2_v1->f, I, I_m, I_t, t, merged_label);
				queue_merge.insert(make_pair(connect_e, t));
			}
		}
		e = e->e_next;
	}

	assert(queue_split.size() == boost::num_vertices(dual_graph) && queue_merge.size() == boost::num_edges(dual_graph));
}

void Partition::compute_energy()
{
	float fidelity_p = 0, fidelity_e = 0;
	int num_e = 0;
	
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); ++it_f) {
		fidelity_p += ((*it_f)->pixels.size() - (*it_f)->semantic_probabilities[(*it_f)->label]);
	}
	fidelity_p /= area;

	Edge* e = edges_head;
	while (e != NULL) {
		++num_e;
		assert(e->unit_edge_cost >= 0);
		bool consistent_label = true;
		if (e->v1_v2->f != nullptr && e->v2_v1->f != nullptr) {
			consistent_label = (e->v1_v2->f->label == e->v2_v1->f->label);
		}
		fidelity_e += fidelity_edge(consistent_label, e->unit_edge_cost, e->length);
		e = e->e_next;
	}
	fidelity_e /= diag_length;
	energy = (1 - params->fidelity_beta) * fidelity_p + params->fidelity_beta * fidelity_e + params->prior_lambda*double(num_e);
}

void Partition::propose_next_operation()
{
	// update vs: a vector of size 2 for merging, 1 for splitting
	vs.clear();

	double merge_gain = (!queue_merge.empty()) ? queue_merge.begin()->second : -10000;
	double split_gain = queue_split.begin()->second;
	double remove_gain = (!queue_bivalent.empty()) ? (*queue_bivalent.begin()).first->energy_gain : -10000;

	if (remove_gain > 0) {
		next_op = REMOVE;
	}
	else if (merge_gain > 0 || split_gain > 0) {
		if (split_gain > merge_gain) {
			vs.push_back(queue_split.begin()->first);
			next_op = SPLIT;
		}
		else {
			Dual_Graph::vertex_descriptor v_i = boost::source(queue_merge.begin()->first, dual_graph), v_j = boost::target(queue_merge.begin()->first, dual_graph);
			vs.push_back(v_i);
			vs.push_back(v_j);
			next_op = MERGE;
		}
	}
	else {
		next_op = STOP;
	}
}

Operator_Type Partition::realize_next_operation(Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, 
	Matrix<double> & I_m, Matrix<double> & I_t, Matrix<double> & I_lab_m, Matrix<double> & I_lab_t, map<int, double> & angles, vector<int> & bin_to_angle_index) {

	++n_iterations;
	propose_next_operation();

	if (next_op == MERGE) {
		merge_operator(I, I_prob, used, I_m, I_t, I_lab_m, I_lab_t, vs, angles, bin_to_angle_index);
	}
	else if (next_op == SPLIT) {
		split_operator(I, I_prob, used, I_m, I_t, I_lab_m, I_lab_t, vs, angles, bin_to_angle_index);
	}
	else if (next_op == REMOVE) {
		if (params->verbose_level > 5) {
			compute_energy();
			trace(params->verbose_level, 5, "energy before : " + std::to_string(energy) + ", estimated energy gain : " + std::to_string((*queue_bivalent.begin()).second));
		}
		remove_bivalent_vertiex((*queue_bivalent.begin()).first, I_m, I_t, I_lab_m, I_lab_t, I, I_prob, used, angles, bin_to_angle_index);
		if (params->verbose_level > 5) {
			compute_energy();
			trace(params->verbose_level, 5, "energy after : " + std::to_string(energy));
		}
	}
	else if (next_op == STOP) {
		std::cout << "** STOP!" << std::endl;
	}

	vs.clear();

	return next_op;
}

void Partition::naive_labeling()
{
	// Label each facet according to the distribution of label histogram by choosing the index that corresponds to the majority
	for (list<Face *>::iterator it = faces.begin(); it != faces.end(); ++it) {
		(*it)->estimate_label(0, area, diag_length);
	}
	// Perform a few iterations to reach a local minimum
	bool stop = true;
	int i = 0;
	do {
		for (list<Face *>::iterator it = faces.begin(); it != faces.end(); ++it) {
			int old_label = (*it)->label;
			(*it)->estimate_label(params->fidelity_beta, area, diag_length);
			if (old_label != (*it)->label) stop = false;
		}
		++i;
	} while ((!stop) && i < 10);
}

pair<double, double> Partition::highest_merge_split_gains() {
	return make_pair(!queue_merge.empty() ? queue_merge.begin()->second : 0, queue_split.begin()->second);
}

void Partition::mark_edges_next_operation()
{
	for (int i = 0; i < edges_next_operation.size(); ++i) edges_next_operation[i].clear();
	edges_next_operation.clear();
	faces_next_operation.clear();
	edges_next_operation.resize(3, vector<pair<Point2d, Point2d>>());
	
	if (vs.empty()) propose_next_operation();

	// save edges to be operated next
	if (next_op == MERGE) {
		bool found;
		Dual_Graph::edge_descriptor connect_e;
		boost::tie(connect_e, found) = boost::edge(vs[0], vs[1], dual_graph);
		assert(found);
		Face* & f_i = dual_graph[vs[0]].f;
		Face* & f_j = dual_graph[vs[1]].f;
		faces_next_operation.push_back(f_i);
		faces_next_operation.push_back(f_j);
		
		for (int id = 0; id < f_i->edges.size(); ++id) {
			Vertex* A, *B;
			for (list<HalfEdge *>::iterator it_h = f_i->edges[id].begin(); it_h != f_i->edges[id].end(); ++it_h) {
				A = (*it_h)->e->v1;
				B = (*it_h)->e->v2;
				edges_next_operation[0].push_back(make_pair(A->pt, B->pt));
			}
		}
		for (int id = 0; id < f_j->edges.size(); ++id) {
			Vertex* A, *B;
			for (list<HalfEdge *>::iterator it_h = f_j->edges[id].begin(); it_h != f_j->edges[id].end(); ++it_h) {
				A = (*it_h)->e->v1;
				B = (*it_h)->e->v2;
				edges_next_operation[1].push_back(make_pair(A->pt, B->pt));
			}
		}
	}
	else if (next_op == SPLIT) {
		Face* f = dual_graph[vs[0]].f;
		faces_next_operation.push_back(f);
		
		// save the best e among candidates
		for (map<pair<Vertex *, Vertex *>, Split_Edge_Property>::iterator it_s = f->split_edges_property.begin(); it_s != f->split_edges_property.end(); ++it_s) {
			Vertex* v1 = it_s->first.first, *v2 = it_s->first.second;
			edges_next_operation[2].push_back(make_pair(v1->pt, v2->pt));
		}
		// save facet edges
		for (int id = 0; id < f->edges.size(); ++id) {
			Vertex* A, *B;
			for (list<HalfEdge *>::iterator it_h = f->edges[id].begin(); it_h != f->edges[id].end(); ++it_h) {
				A = (*it_h)->e->v1;
				B = (*it_h)->e->v2;
				edges_next_operation[0].push_back(make_pair(A->pt, B->pt));
			}
		}
	} 
	else if (next_op == REMOVE) {
		Vertex* v = (*queue_bivalent.begin()).first;
		Edge* e1 = v->directions[0].second->e;
		Edge* e2 = v->directions[1].second->e;
		edges_next_operation[2].push_back(make_pair(e1->v1->pt, e1->v2->pt));
		edges_next_operation[2].push_back(make_pair(e2->v1->pt, e2->v2->pt));
	}
}

void Partition::merge_gain(const list<Dual_Merge_Graph::vertex_descriptor> & face_vertices_list, const map<list<Edge *>, bool>& shared_edges, 
	Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t, double & gain, Label_ID & merged_label, bool ignore_boundary_conditions, bool verbose)
{
	list<vector<double>> faces_prob_hists;
	list<Label_ID> labels;
	list<int> num_pixels;
	list<vector<Boundary_Edge_Info>> boundary_edges_list;

	for (auto graph_vertex: face_vertices_list) {
		Merge_Graph_Vertex_Property vertex_property = merge_graph[graph_vertex];
		faces_prob_hists.push_back(vertex_property.semantic_probabilities);
		labels.push_back(vertex_property.label);
		num_pixels.push_back(vertex_property.n_pixels);

		vector<Boundary_Edge_Info> temp;
		for (auto it = vertex_property.boundary_edges.begin(); it != vertex_property.boundary_edges.end(); ++it) {
			if (it ->label == INVALID_LABEL) continue;
			temp.emplace_back(*it);
		}

		Dual_Merge_Graph::adjacency_iterator adj_vit, adj_vend;
		boost::tie(adj_vit, adj_vend) = adjacent_vertices(graph_vertex, merge_graph);
		for (auto it = adj_vit; it != adj_vend; ++it) {
			if (find(face_vertices_list.begin(), face_vertices_list.end(), *it) != face_vertices_list.end()) {
				continue;
			}
			auto found = boost::edge(*it, graph_vertex, merge_graph);
			const list<Edge *> & common_edges = merge_graph[found.first].shared_edges;
			for (auto it_e = common_edges.begin(); it_e != common_edges.end(); ++it_e) {
				temp.emplace_back(*it_e, merge_graph[*it].label);
			}
		}
		boundary_edges_list.push_back(temp);
	}

	merge_gain(faces_prob_hists, boundary_edges_list, labels, num_pixels, shared_edges, gain, merged_label, ignore_boundary_conditions, verbose);
}


void Partition::merge_gain(Face* f_1, Face* f_2, Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t, double & gain, Label_ID & merged_label, bool verbose)
{
	vector<double> & feature_1 = f_1->semantic_probabilities, &feature_2 = f_2->semantic_probabilities;
	int n_1 = f_1->pixels.size(), n_2 = f_2->pixels.size();
	int label_1 = f_1->label, label_2 = f_2->label;
	int n_edges_1 = 0, n_edges_2 = 0;
	for (int id = 0; id < f_1->edges.size(); ++id) n_edges_1 += f_1->edges[id].size();
	for (int id = 0; id < f_2->edges.size(); ++id) n_edges_2 += f_2->edges[id].size();

	map<list<Edge *>, bool> shared_edges; 
	list<Edge *> temp;
	vector<Boundary_Edge_Info> other_adj_edges_1, other_adj_edges_2;
	other_adj_edges_1.reserve(n_edges_1);
	other_adj_edges_2.reserve(n_edges_2);
	for (int id = 0; id < f_1->edges.size(); ++id) {
		for (list<HalfEdge *>::iterator it_h = f_1->edges[id].begin(); it_h != f_1->edges[id].end(); ++it_h) {
			Face* f_adj = (*it_h)->opposite()->f;
			if (f_adj != nullptr && f_adj->id_face == f_2->id_face) temp.emplace_back((*it_h)->e);
			else if (f_adj != nullptr) {
				other_adj_edges_1.emplace_back((*it_h)->e, f_adj->label);
			}
		}
	}
	for (int id = 0; id < f_2->edges.size(); ++id) {
		for (list<HalfEdge *>::iterator it_h = f_2->edges[id].begin(); it_h != f_2->edges[id].end(); ++it_h) {
			Face* f_adj = (*it_h)->opposite()->f;
			if (f_adj != f_1 && f_adj != nullptr) {
				other_adj_edges_2.emplace_back((*it_h)->e, f_adj->label);
			}
		}
	}
	if (label_1 == label_2) shared_edges[temp] = true;
	else shared_edges[temp] = false;
	merge_gain({feature_1, feature_2}, {other_adj_edges_1, other_adj_edges_2}, {label_1, label_2}, {n_1, n_2}, shared_edges, gain, merged_label, false, verbose);
}

void Partition::merge_gain(const list<vector<double>>& faces_prob_hists, const list<vector<Boundary_Edge_Info>> & adj_face_edges_1ist, 
	const list<Label_ID>& labels, const list<int>& num_pixels, const map<list<Edge *>, bool>& shared_edges, 
	double & gain, Label_ID & merged_label, bool ignore_boundary_conditions, bool verbose)
{
	/*
	Merge gain = energy before - energy after (we minimize the energy)
	The energy is expressed as the sum of two terms: fidelity and prior. 
	The fidelity term should be the sum of two parts: agreement with probability map, and variance multiplied by the number of pixels inside facet.
	The prior term penalizes number of edges.
	*/

	/* edge information */
	double grad_cost = 0, shared_length = 0;
	list<double> shared_edges_costs;
	for (auto& edges: shared_edges) {
		boost::tie(grad_cost, shared_length) = joint_grad_cost_and_length(edges.first);
		shared_edges_costs.push_back(fidelity_edge(edges.second, grad_cost, shared_length));
	}
	
	/* compute merged facet probability histogram */
	int n_total = accumulate(num_pixels.begin(), num_pixels.end(), 0);
	vector<double> merged_hist;
	Face::compute_merged_histogram(faces_prob_hists, merged_hist);

	vector<vector<Adj_Edge_Info>*> edges_to_neighbor_labels;
	vector<Adj_Edge_Info> temp;
	if (!ignore_boundary_conditions) {
		int num_edges = 0;
		for (auto& adj_face_edges : adj_face_edges_1ist) num_edges += adj_face_edges.size();
		temp.reserve(num_edges);
		for (auto& adj_face_edges : adj_face_edges_1ist) {
			for (int i = 0; i < adj_face_edges.size(); ++i) {
				temp.push_back({ adj_face_edges[i].e, adj_face_edges[i].label });
			}
		}
		edges_to_neighbor_labels.push_back(&temp);
	}
	
	Face::estimate_label(merged_hist, n_total, merged_label, params->fidelity_beta, area, diag_length, edges_to_neighbor_labels, ignore_boundary_conditions);
	
	/* compute local energy for facets before merging */
	double fidelity_before(0), fidelity_p_before(0), fidelity_e_before(0);
	int num_e_reduce = shared_edges.size();

	auto it_f = faces_prob_hists.begin();
	auto it_p = num_pixels.begin();
	auto it_l = labels.begin();
	while (it_f!= faces_prob_hists.end()) {
		fidelity_p_before += (*it_p) - (*it_f)[*it_l];
		++it_f;
		++it_p;
		++it_l;
	}
	fidelity_p_before /= double(area);

	fidelity_e_before = accumulate(shared_edges_costs.begin(), shared_edges_costs.end(), 0);
	it_l = labels.begin();
	for (auto& adj_face_edges : adj_face_edges_1ist) {
		for (auto it = adj_face_edges.begin(); it != adj_face_edges.end(); ++it) {
			bool same_label = (it->label == *it_l);
			Edge* adj_e = it->e;
			fidelity_e_before += fidelity_edge(same_label, adj_e->unit_grad_cost(), adj_e->length);
		}
		++it_l;
	}
	fidelity_e_before /= diag_length;
	
	fidelity_before = (1- params->fidelity_beta)*fidelity_p_before + params->fidelity_beta*fidelity_e_before;


	/* compute local energy after merging  */
	double fidelity_after(0), fidelity_p_after(0), fidelity_e_after(0);
	fidelity_p_after += n_total - merged_hist[merged_label];
	fidelity_p_after /= double(area);

	for (int i = 0; i < edges_to_neighbor_labels.size(); ++i) {
		for (auto it = edges_to_neighbor_labels[i]->begin(); it != edges_to_neighbor_labels[i]->end(); ++it) {
			bool same_label = (it->label == *it_l);
			Edge* adj_e = it->e;
			fidelity_e_after += fidelity_edge(same_label, adj_e->unit_grad_cost(), adj_e->length);
		}
	}
	fidelity_e_after /= diag_length;

	fidelity_after = (1- params->fidelity_beta)*fidelity_p_after + params->fidelity_beta*fidelity_e_after;


	/* compute change in energy  */
	gain = fidelity_before - fidelity_after + params->prior_lambda*double(num_e_reduce);

	if (verbose) {
		std::cout << "fidelity_p_before, fidelity_e_before, num_e_reduce; fidelity_p_after, fidelity_e_after = " 
			<< fidelity_p_before << ", " << fidelity_e_before << ", " << num_e_reduce << "; " << fidelity_p_after << fidelity_e_after << std::endl;
	}
}


pair<double, double> Partition::joint_grad_cost_and_length(const list<Edge *> & edges) {
	double mean_grad_cost = 0;
	double shared_length = 0;
	for (auto it = edges.begin(); it != edges.end(); ++it) {
		Edge * e = (*it);
		double unit_grad_cost = e->unit_grad_cost();
		mean_grad_cost += unit_grad_cost*e->length;
		shared_length += e->length;
	}
	mean_grad_cost /= max(shared_length, 1e-5);
	return pair<double,double>(mean_grad_cost, shared_length);
}


void Partition::detect_line_segments_all_facets(Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t, bool is_color_image)
{
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		detect_line_segments(I, I_m, I_t, *it_f, is_color_image);
	}
}

int Partition::detect_line_segments(Matrix<uchar> & I, Matrix<double> & I_m, Matrix<double> & I_t, Face* f, bool is_color_image)
{
	/*
	Detect line segments by region growing. We keep two vectors, alpha (weighted average gradient direction) and nv (normal vector to the line segment).
	As the line segment grows longer, alpha and nv should get closer because the influence of discretization decreases.
	Due to the discrete nature of image grid, the region growing should look 2 steps ahead instead of 1, that is,
	the neighbor pixel is accepted if:
	1) the neighbor pixel has similar gradient direction to alpha.
	*/
	if (f->check_split == false) return 0;

	int n_lines = 0;
	// set parameters
	double length_thresh = params->split_params.length_thresh;
	double spatial_thresh = 0.71;
	double value_thresh = is_color_image ? params->split_params.m_value_thresh : params->split_params.lab_m_value_thresh;
	double bin_width = params->split_params.bin_width;
	double distance_to_edge_thresh = params->split_params.distance_to_edge;
	double angle_thresh = params->split_params.angle;
	double p = angle_thresh / PI;
	double min_dist_angle_score = params->split_params.min_dist_angle_score;

	vector<set<Point2i, pointcomp>> bins = is_color_image ? f->active_gradients : f->active_lab_gradients;

	int n_active_pixels = 0;
	for (int i = 0; i < bins.size(); i++) n_active_pixels += bins[i].size();
	
	int m_rows = static_cast<int>(I_m.rows);
	int m_cols = static_cast<int>(I_m.cols);

	bool exist(false);
	Point2i seed_pixel;
	while (true) {
		// select a seed from high active gradient pixels
		exist = false;
		for (auto rit = bins.rbegin(); rit != bins.rend() ; ++rit) {
			if (!rit->empty()) {
				auto it = rit->begin();
				seed_pixel = *it;
				rit->erase(it); // remove from bins
				exist = true;
				break;
			}
		}
		if (!exist) break;

		// create container
		vector<Point2i> region;
		double grad_angle;
		Geometry::region_grow(seed_pixel, I_m, I_t, region, grad_angle, bins, value_thresh, bin_width, angle_thresh);
		/* if the region is too small, reject */
		int reg_size = region.size();
		if (reg_size < 3) continue;
		/* construct rectangular approximation for the region */
		rect rec;
		double region_angle = grad_angle + PI/2; // level line angle
		if (region_angle > PI) region_angle -= M_2__PI;
		Geometry::region2rect(region, I_m, region_angle, angle_thresh, p, &rec);

		bool accept_region_density = Geometry::refine(region, &reg_size, I_m, bins, value_thresh, bin_width, region_angle, grad_angle,
			angle_thresh, p, &rec, I_t, params->lsd_density);

		Vec2d nv(rec.dy, -rec.dx);
		if (nv[0] * cos(grad_angle) + nv[1] * sin(grad_angle) < 0) nv = -nv;
		// add vertex coordinate offset if I_m is computed on pixel corners
		Vec2d line_dir(rec.dx, rec.dy);
		Point2d end1 = Point2d(rec.x1, rec.y1) + offset;
		Point2d end2 = Point2d(rec.x2, rec.y2) + offset;
		if (line_dir.dot(end2 - end1) < 0) line_dir = -line_dir;
		end1 = end1 - 0.5*Point2d(line_dir[0], line_dir[1]);
		end2 = end2 + 0.5*Point2d(line_dir[0], line_dir[1]);
		Point2d center = 0.5*(end1 + end2);
		double l_length = cv::norm(end2 - end1);

		// compute score and return unrelated pixels back to bins
		double unit_grad_cost = 0, unit_score = 1;
		if (is_color_image && accept_region_density) {
			double score_sum = 0, sum_overlap = 0;
			double mean_angle_diff = 0;
			for (int k = 0; k < reg_size; ++k) {
				Point2i p_k = region[k];
				Point2d disp = (Point2d)p_k + offset - end1;
				double transverse_dist = fabs(disp.ddot(nv));
				int i = m_rows - 1 - p_k.y; // rows
				int j = p_k.x; // columns
				double dbin = I_m(i, j) / bin_width;
				int bin = floor(dbin);
				if (bin >= cdf_grad_m.size()) bin = cdf_grad_m.size() - 1;
				double d_angle = Geometry::angle_diff_abs(I_t(i, j), grad_angle);
				if (d_angle < angle_thresh) { 
					if (transverse_dist > spatial_thresh) continue;
					double cdf_prev_bin = bin > 0 ? cdf_grad_m[bin - 1] : 0;
					double raw_pixel_score = cdf_prev_bin + (dbin - bin)*(cdf_grad_m[bin] - cdf_prev_bin);  // Give a raw pixel score in [0, 1]
					double axial_dist = disp.ddot(line_dir);
					double d_min_pos = axial_dist - 0.5, d_max_pos = axial_dist + 0.5;
					double overlap = 1;
					if (d_max_pos > l_length) overlap -= d_max_pos - l_length;
					if (d_min_pos < 0) overlap -= 0 - d_min_pos;
					overlap = max(overlap, 0.);
					double pixel_score = compute_pixel_score(raw_pixel_score, d_angle);
					sum_overlap += overlap;
					mean_angle_diff += overlap*d_angle;
					score_sum += overlap*pixel_score;
				}
				else if (p_k != seed_pixel) {
					int id_bin = jclamp(0, floor((I_m(i, j) - value_thresh) / bin_width), bins.size() - 1);
					bins[id_bin].emplace(p_k);
					/* remove point from the region by replacing with the last element*/
					region[k].x = region[reg_size - 1].x; /* if k==*reg_size-1 copy itself */
					region[k].y = region[reg_size - 1].y;
					--reg_size;
					--k; /* to avoid skipping one point */
				}
			}
			region.resize(reg_size); // resize to remove repeated points
			mean_angle_diff /= l_length;
			/* estimate an average score */
			unit_score = score_sum / max(sum_overlap, l_length);
			unit_grad_cost = compute_unit_grad_cost(1 - unit_score, mean_angle_diff, 0);
		}

		double t_0 = -1e5, t_1 = 1e5;
		double min_dist_to_edge = FLT_MAX;
		bool overlap_with_edge_region = false;
		double covered_length = 0;
		bool align_with_edge = false;
		int ii = 0;
		while (ii < f->edges.size() && !align_with_edge) {
			for (list<HalfEdge *>::iterator it_h = f->edges[ii].begin(); it_h != f->edges[ii].end(); ++it_h) {
				Vec2d & e_normal = (*it_h)->e->normal;
				bool aligned_normal = is_color_image ? (nv[0]*e_normal[0] + nv[1]*e_normal[1] > 0.9848) : (fabs(nv[0]*e_normal[0] + nv[1]*e_normal[1]) > 0.9848); // within +-10 degrees
				if (aligned_normal) {
					Point2d & pB = (*it_h)->e->v2->pt;
					double axial_disp = e_normal[1] * (center.x - pB.x) + (-e_normal[0]) * (center.y - pB.y);
					double d_1 = fabs(e_normal.dot(pB - end1)), d_2 = fabs(e_normal.dot(pB - end2));
					bool between_edge = (axial_disp <= 0 && axial_disp >= -(*it_h)->e->length);
					if (between_edge) {
						double d = min(d_1, d_2) - 0.5*rec.width*fabs(e_normal.dot(nv));
						if (d < min_dist_to_edge) min_dist_to_edge = d;
						if (d < 0.5*(*it_h)->e->width) overlap_with_edge_region = true;
					}
					double dist_to_AB = max(d_1, d_2);
					if (dist_to_AB <= distance_to_edge_thresh) {
						if (between_edge) {
							covered_length += min(-axial_disp, 0.5*l_length) + min(axial_disp + (*it_h)->e->length, 0.5*l_length);
							if (is_color_image) align_with_edge = (covered_length > 0.95*l_length) || (dist_to_AB < 1.0 || (1-(*it_h)->e->unit_grad_cost())*(*it_h)->e->length > (1-unit_grad_cost)*l_length);
							else align_with_edge = (covered_length > 0.9*l_length);
						}
						else {
							if (axial_disp - 0.5*l_length < 0) covered_length += -(axial_disp - 0.5*l_length);
							else if (axial_disp + 0.5*l_length > -(*it_h)->e->length) covered_length += axial_disp + 0.5*l_length + (*it_h)->e->length;
						}
					}
				}
				if (align_with_edge) break;

				double determinant = line_dir[0] * (-e_normal[0]) - line_dir[1] * e_normal[1];
				if (fabs(determinant) > 1e-4) {// If not parallel, compute intersection.
					Point2d & pB = (*it_h)->e->v2->pt;
					double y_i = -(line_dir[0] * (pB.y - center.y) - line_dir[1] * (pB.x - center.x)) / determinant;
					if (y_i <= 0 && y_i >= -(*it_h)->e->length) {
						double x_i = ((pB.x - center.x) * (-e_normal[0]) - (pB.y - center.y) * e_normal[1]) / determinant;
						t_1 = (x_i > 0 && x_i < t_1) ? x_i : t_1;
						t_0 = (x_i <= 0 && x_i > t_0) ? x_i : t_0;
					}
				}
			}
			ii++;
		}

		/* Segments are accepted according to the following criteria:
		Region density: region should contain dense enough points
		Location: should not align with a nearby edge
		Score:  unit gradient score & length score
		Lenght in the facet: percentage of segment support when extended till the facet boundary */
		double support_percent = l_length / (t_1 - t_0);
		double length_cost = 1 / (1 + exp(5 * (2 * l_length / length_thresh - 1))); // = 0.5 when it's half length_thresh
		bool accept = (accept_region_density && (!align_with_edge) && ((1- unit_grad_cost)*(1 - length_cost) > 0.5 || support_percent >= 0.3));
		/* remove pixels if the segment is accepted or the region does not overlap an edge */
		bool remove_pixels_from_candidates = (accept || (!align_with_edge && !overlap_with_edge_region));
		if (accept) {
			/* Crop line segment*/
			if (t_0 <= -1e4 || t_1 >= 1e4) { // segment falls outside the facet
				if (t_0 <= -1e4 && t_1 < 1e4) end1 = center + t_1*Point2d(line_dir[0], line_dir[1]);
				if (t_1 >= 1e4 && t_0 > -1e4) end2 = center + t_0*Point2d(line_dir[0], line_dir[1]);
			}
			else {
				if (t_0 > -0.5*l_length) end1 = center + t_0*Point2d(line_dir[0], line_dir[1]);
				if (t_1 < 0.5*l_length) end2 = center + t_1*Point2d(line_dir[0], line_dir[1]);
			}
			double m_mean = 0;
			Segment* s = new Segment(id_segments++, end1.x, end1.y, end2.x, end2.y, rec.width, unit_grad_cost, Point2d(nv[0], nv[1]), false, is_color_image);
			f->detected_segments.insert(s);
			s->pixels.reserve(reg_size);
			/*save pixels*/
			for (int k = 0; k < reg_size; ++k) s->pixels.push_back(region[k]);
			++n_lines;
		}
		if (remove_pixels_from_candidates) { // remove pixels
			n_active_pixels -= reg_size;
			for (auto it = region.begin(); it != region.end(); ++it) {
				int i = m_rows - 1 - it->y; // rows
				int j = it->x; // columns
				int id_bin = jclamp(0, floor((I_m(i, j) - value_thresh) / bin_width), bins.size() - 1);
				if (is_color_image) f->active_gradients[id_bin].erase(*it);
				else f->active_lab_gradients[id_bin].erase(*it);
			}
		}
	}
	return n_lines;
}


void Partition::extend_detected_segments(Face* f, Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, 
	Matrix<double> & I_m, Matrix<double> & I_t, map<int, double> & angles, vector<int> & bin_to_angle_index, double & gain)
{
	vector<Segment*> selected_split_segments;
	f->split_vertex_chains.clear();
	f->split_edges_property.clear();
	f->split_vertex_outer_edge_location.clear();

	// Step 1: Select candidate line segments
	int n_segs = 0;
	for (auto it = f->detected_segments.begin(); it != f->detected_segments.end(); ++it) {
		selected_split_segments.emplace_back(*it);
		++n_segs;
		if (n_segs >= 10) break;
	}

	if (f->pixels.size() < 2 || selected_split_segments.empty() || f->check_split == false) {
		f->split_gain = gain = -10000;
		return;
	}

	// Step 2: Propagate segments and record all subfacets
	vector<Ray *> rays = vector<Ray *>();
	Partition subgraph(rows, cols, params);
	Segment_Regularization_Tree subtree;
	int n_init_vertices = 0;
	Propagation::propagate_inside_polygon(f->edges, &subgraph, &subtree, angles, bin_to_angle_index, rays, selected_split_segments, n_init_vertices);
	
	if (subgraph.faces.size() < 2) {
		f->split_gain = gain = -10000;
		return;
	}

	// inherit cdf from parent partition
	subgraph.cdf_grad_m = this->cdf_grad_m;
	subgraph.cdf_lab_grad_m = this->cdf_lab_grad_m;
	subgraph.find_pixels_inside_all_facets();
	subgraph.weight_all_edges(used, I_m, I_t);
	subgraph.compute_feature_all_facets(I_prob, I);
	subgraph.naive_labeling();
	subgraph.build_merge_graph(I, I_m, I_t);

	/* recover used gradients for inner edges of subgraph*/
	Edge* sub_e = subgraph.edges_head;
	while (sub_e != nullptr) {
		if (sub_e->type == INNER_EDGE) {
			const vector<Point2i> & grad_pixels = sub_e->region;
			for (int k = 0; k < grad_pixels.size(); ++k) {
				int i = rows - 1 - grad_pixels[k].y; // rows
				int j = grad_pixels[k].x; // columns
				used(i, j) = false;
			}
		}
		sub_e = sub_e->e_next;
	}

	// Step 3: Group subfacets
	set<pair<Dual_Merge_Graph::edge_descriptor, double>, queuecomp> & sub_queue = subgraph.queue_merge_graph;
	Dual_Merge_Graph & sub_merge_graph = subgraph.merge_graph;
	if (sub_queue.rbegin()->second > 0) {
		f->split_gain = gain = -10000;
		return;
	}
	while (!sub_queue.empty() && sub_queue.begin()->second > 0) {

		Dual_Merge_Graph::edge_descriptor e_top = sub_queue.begin()->first;
		Dual_Merge_Graph::vertex_descriptor v_i = boost::source(e_top, sub_merge_graph), v_j = boost::target(e_top, sub_merge_graph);
		assert(v_i != v_j);
		Label_ID label_i = sub_merge_graph[v_i].label, label_j = sub_merge_graph[v_j].label;
		Label_ID label_ij = sub_merge_graph[e_top].merged_label;

		/* record adjacency iterators of v_j */
		Dual_Merge_Graph::adjacency_iterator vit_j, vend_j;
		boost::tie(vit_j, vend_j) = adjacent_vertices(v_j, sub_merge_graph);

		/* add new edges from adjacent vertices of v_j to v_i */
		while (vit_j != vend_j) {
			Dual_Merge_Graph::vertex_descriptor v = *vit_j;
			bool found_vi_v, found_vj_v;
			Dual_Merge_Graph::edge_descriptor e_vi_v, e_vj_v;
			boost::tie(e_vj_v, found_vj_v) = boost::edge(v_j, v, sub_merge_graph);
			// remove queue entry related to v_j
			sub_queue.erase(make_pair(e_vj_v, sub_merge_graph[e_vj_v].t));
			if (v != v_i) {
				boost::tie(e_vi_v, found_vi_v) = boost::edge(v_i, v, sub_merge_graph);
				if (found_vi_v) {
					list<Edge *> & shared_e_vi_v = sub_merge_graph[e_vi_v].shared_edges, &shared_e_vj_v = sub_merge_graph[e_vj_v].shared_edges;
					shared_e_vi_v.insert(shared_e_vi_v.begin(), shared_e_vj_v.begin(), shared_e_vj_v.end());
				}
				else {
					bool inserted;
					boost::tie(e_vi_v, inserted) = boost::add_edge(v_i, v, sub_merge_graph);
					sub_merge_graph[e_vi_v].shared_edges = sub_merge_graph[e_vj_v].shared_edges;
				}
			}
			++vit_j;
		}
		
		/* update the graph vertex by merging v_j into v_i */
		list<Face *> & superfacet_i = sub_merge_graph[v_i].superfacet;
		list<Face *> & superfacet_j = sub_merge_graph[v_j].superfacet;
		superfacet_i.insert(superfacet_i.begin(), superfacet_j.begin(), superfacet_j.end());
		list<Boundary_Edge_Info> & boundary_edges_i = sub_merge_graph[v_i].boundary_edges;
		list<Boundary_Edge_Info> & boundary_edges_j = sub_merge_graph[v_j].boundary_edges;
		boundary_edges_i.insert(boundary_edges_i.begin(), boundary_edges_j.begin(), boundary_edges_j.end());
		sub_merge_graph[v_i].n_edges = sub_merge_graph[v_i].n_edges + sub_merge_graph[v_j].n_edges - 2*sub_merge_graph[e_top].shared_edges.size();

		/* compute semantic_probabilities of f_ij */
		vector<double> feature_ij;
		vector<int> color_hist_ij;
		int n_total = sub_merge_graph[v_i].n_pixels + sub_merge_graph[v_j].n_pixels;
		Face::compute_merged_histogram({sub_merge_graph[v_i].semantic_probabilities, sub_merge_graph[v_j].semantic_probabilities}, feature_ij);

		sub_merge_graph[v_i].semantic_probabilities = feature_ij;
		sub_merge_graph[v_i].label = label_ij;
		sub_merge_graph[v_i].n_pixels = n_total;
		sub_merge_graph[v_i].perimeter += sub_merge_graph[v_j].perimeter;
		for (list<Edge *>::iterator it_se = sub_merge_graph[e_top].shared_edges.begin(); it_se != sub_merge_graph[e_top].shared_edges.end(); ++it_se) {
			sub_merge_graph[v_i].perimeter -= 2 * (*it_se)->length;
			(*it_se)->disable();
		}

		boost::clear_vertex(v_j, sub_merge_graph);
		// update edge properties of edges of v_i (remove edges connecting to v_j before this step)
		Dual_Merge_Graph::adjacency_iterator vit_i, vend_i;
		boost::tie(vit_i, vend_i) = adjacent_vertices(v_i, sub_merge_graph);
		while (vit_i != vend_i) {
			Dual_Merge_Graph::vertex_descriptor v = *vit_i;
			if (v != v_j) {
				Dual_Merge_Graph::edge_descriptor e_vi_v = boost::edge(v_i, v, sub_merge_graph).first;
				// remove corresponding entry in the queue
				sub_queue.erase(make_pair(e_vi_v, sub_merge_graph[e_vi_v].t));

				map<list<Edge *>, bool> shared_e_vi_v;
				shared_e_vi_v[sub_merge_graph[e_vi_v].shared_edges] = (label_i==sub_merge_graph[v].label);
				Label_ID & merged_label = sub_merge_graph[e_vi_v].merged_label;
				subgraph.merge_gain({v_i, v}, shared_e_vi_v, I, I_m, I_t, sub_merge_graph[e_vi_v].t, merged_label, false);
				// insert updated queue entry
				sub_queue.insert(make_pair(e_vi_v, sub_merge_graph[e_vi_v].t));

				// Update other merge queue entries associated with this neighbor face
				if (label_ij != label_i || label_ij != label_j) {
					Dual_Merge_Graph::adjacency_iterator vit, vend;
					boost::tie(vit, vend) = adjacent_vertices(v, sub_merge_graph);
					while (vit != vend) {
						Dual_Merge_Graph::vertex_descriptor v_adj = *vit;
						if (v_adj != v_i && v_adj != v_j) {
							Dual_Merge_Graph::edge_descriptor e_v_adj;
							bool found_old;
							boost::tie(e_v_adj, found_old) = boost::edge(*vit, v, sub_merge_graph);
							sub_queue.erase(make_pair(e_v_adj, sub_merge_graph[e_v_adj].t));
							list<Edge *> & common_edges = sub_merge_graph[e_v_adj].shared_edges;
							map<list<Edge *>, bool> shared_e_vit_v;
							shared_e_vit_v[common_edges] = (sub_merge_graph[v].label==sub_merge_graph[v_adj].label);
							subgraph.merge_gain({v, v_adj}, shared_e_vit_v, I, I_m, I_t, sub_merge_graph[e_v_adj].t, sub_merge_graph[e_v_adj].merged_label, false);
							sub_queue.insert(make_pair(e_v_adj, sub_merge_graph[e_v_adj].t));
						}
						++vit;
					}
				}
			}
			++vit_i;
		}
		boost::remove_vertex(v_j, sub_merge_graph);
	}
	if (boost::num_edges(sub_merge_graph) == 0) {
		f->split_gain = gain = -10000;
		return;
	}
	
	// Finally we collect shared edges between resulting superfacets
	set<Edge *> split_edges;
	map<list<Edge *>, bool> split_edges_map;
	for (auto& merge_op : sub_queue) {
		list<Edge *> & shared_edges = sub_merge_graph[merge_op.first].shared_edges;
		split_edges.insert(shared_edges.begin(), shared_edges.end());
		split_edges_map[shared_edges] = false;
	}
	
	vector<list<HalfEdge *>> ordered_edge_chains;
	map<Vertex*, pair<int, int>> parent_outer_edge_location_subgraph; //record location of vertices that touch the polygon boundary, including those coincide with polygon vertices

	// form chains that stops whenever touching the polygon boundary, or reaching a vertex whose other edges are all visited
	while (split_edges.size() > 0) {
		ordered_edge_chains.push_back(list<HalfEdge *>());
		HalfEdge* h = (*split_edges.begin())->v1_v2;
		Vertex* v1 = h->v1_v2 ? h->e->v1 : h->e->v2;
		Vertex* v2 = h->v1_v2 ? h->e->v2 : h->e->v1;
		split_edges.erase(h->e);
		ordered_edge_chains.back().emplace_back(h);
		for (int i = 0; i < v1->directions.size(); i++) {
			HalfEdge* h_v1 = v1->directions[i].second;
			if (h_v1->e->type == OUTER_EDGE && h_v1->f != nullptr) {
				Outer_Edge* oe = static_cast<Outer_Edge *>(h_v1->e);
				parent_outer_edge_location_subgraph.insert(make_pair(v1, oe->init_boundary_location));
				break;
			}
		}

		bool stop = false;
		Vertex* v_next = v2;
		do {
			for (int i = 0; i < v_next->directions.size(); i++) {
				HalfEdge* h_next = v_next->directions[i].second;
				if (h_next->e->type == OUTER_EDGE && h_next->f != nullptr) {
					Outer_Edge* oe = static_cast<Outer_Edge *>(h_next->e);
					parent_outer_edge_location_subgraph.insert(make_pair(v_next, oe->init_boundary_location));
					stop = true; 
					break;
				}
			}
			if (!stop) {
				bool added=false;
				for (int i = 0; i < v_next->directions.size(); i++) {
					Edge* e = v_next->directions[i].second->e;
					if (split_edges.count(e)) {
						if (v_next == e->v1) {
							h = e->v1_v2;
							v_next = e->v2;
						}
						else {
							h = e->v2_v1;
							v_next = e->v1;
						}
						split_edges.erase(e);
						ordered_edge_chains.back().emplace_back(h);
						added = true;
						break;
					}
				}
				if (v_next == v1 || !added) {
					stop = true;
				}
			}
		} while (!stop);
		// now check along the other side if necessary
		stop = false;
		if (v_next != v1) {
			v_next = v1;
			do {
				for (int i = 0; i < v_next->directions.size(); i++) {
					HalfEdge* h_next = v_next->directions[i].second;
					if (h_next->e->type == OUTER_EDGE && h_next->f != nullptr) {
						Outer_Edge* oe = static_cast<Outer_Edge *>(h_next->e);
						parent_outer_edge_location_subgraph.insert(make_pair(v_next, oe->init_boundary_location));
						stop = true;
						break;
					}
				}
				if (!stop) {
					bool added=false;
					for (int i = 0; i < v_next->directions.size(); i++) {
						Edge* e = v_next->directions[i].second->e;
						if (split_edges.count(e)) {
							if (v_next == e->v1) {
								h = e->v2_v1;
								v_next = e->v2;
							}
							else {
								h = e->v1_v2;
								v_next = e->v1;
							}
							split_edges.erase(e);
							ordered_edge_chains.back().push_front(h);
							added = true;
							break;
						}
					}
					if (!added) stop = true;
				}
			} while (!stop);
		}
	}

	// Step 4: Create new vertices and make ordered lists while skipping unneccessary vertices
	map<Vertex*, pair<int, int>> & split_vertex_outer_edge_location = f->split_vertex_outer_edge_location; // location of vertices that touch the polygon boundary, including those coincide with polygon vertices
	map<pair<Vertex *, Vertex *>, Split_Edge_Property> & split_edges_property = f->split_edges_property;
	vector<vector<Vertex *>> & ordered_vertex_chains = f->split_vertex_chains;
	ordered_vertex_chains.resize(ordered_edge_chains.size(), vector<Vertex *>());

	map<Point2d, Vertex*, pointcomp> added_vertices;
	for (size_t i = 0; i < ordered_edge_chains.size(); ++i) {
		list<HalfEdge *>::iterator it_h = ordered_edge_chains[i].begin();
		Vertex* v_start = (*it_h)->v1_v2 ? (*it_h)->e->v1 : (*it_h)->e->v2;
		Vertex* v_end = ordered_edge_chains[i].back()->v1_v2 ? ordered_edge_chains[i].back()->e->v2 : ordered_edge_chains[i].back()->e->v1;
		if (added_vertices.count(v_start->pt)) {
			ordered_vertex_chains[i].emplace_back(added_vertices.at(v_start->pt));
		}
		else if (parent_outer_edge_location_subgraph.count(v_start) > 0) { // When the first vertex touches the polygon boundary
			pair<int, int> contour_edge_id = parent_outer_edge_location_subgraph.at(v_start);
			if (v_start->id_vertex < n_init_vertices) { // if it coincide with a polygon vertex, we point to the existing vertex
				list<pair<Vertex *, Vertex_Type> >::iterator it_v = f->vertices[contour_edge_id.first].begin();
				std::advance(it_v, contour_edge_id.second);
				assert(it_v->first->pt == v_start->pt);
				ordered_vertex_chains[i].emplace_back(it_v->first);
			}
			else { // create a new vertex when it does not coincide with an existing polygon vertex
				Vertex* v_outer = new Vertex(id_vertices, v_start->pt.x, v_start->pt.y);
				added_vertices[v_start->pt] = v_outer;
				ordered_vertex_chains[i].emplace_back(v_outer);
				split_vertex_outer_edge_location.insert(make_pair(v_outer, contour_edge_id));
			}
		}
		else { // When the first vertex does not touch the polygon boundary
			Vertex* v_inner = new Vertex(id_vertices, v_start->pt.x, v_start->pt.y);
			added_vertices[v_start->pt] = v_inner;
			ordered_vertex_chains[i].emplace_back(v_inner);
		}

		while (it_h != ordered_edge_chains[i].end()) {
			if (*it_h == ordered_edge_chains[i].back()) {
				Vertex* v_prev = ordered_vertex_chains[i].back();
				Vertex* v_curr = nullptr;
				if (v_start == v_end) { // closed loop
					v_curr = ordered_vertex_chains[i].front();
				}
				else if (added_vertices.count(v_end->pt)) {
					v_curr = added_vertices.at(v_end->pt);
					ordered_vertex_chains[i].emplace_back(v_curr);
				}
				else if (parent_outer_edge_location_subgraph.count(v_end)) {
					pair<int, int> contour_edge_id = parent_outer_edge_location_subgraph.at(v_end);
					if (v_end->id_vertex < n_init_vertices) { // if it coincide with a polygon vertex, we point to the corresponding polygon vertex
						list<pair<Vertex *, Vertex_Type> >::iterator it_v = f->vertices[contour_edge_id.first].begin();
						std::advance(it_v, contour_edge_id.second);
						if (it_v->first->pt != v_end->pt) {
							std::cout << contour_edge_id.second << std::endl;
							for (auto v_pair : f->vertices[contour_edge_id.first]) {
								std::cout << v_pair.first->pt << "  ";
							}
							std::cout << std::endl;
							std::cout << v_end->pt << std::endl;
						}
						assert(it_v->first->pt == v_end->pt);
						v_curr = it_v->first;
						ordered_vertex_chains[i].emplace_back(v_curr);
					}
					else { // create a new vertex when it does not coincide with an existing polygon vertex
						v_curr = new Vertex(id_vertices, v_end->pt.x, v_end->pt.y);
						ordered_vertex_chains[i].emplace_back(v_curr);
						split_vertex_outer_edge_location.insert(make_pair(v_curr, contour_edge_id));
					}
				}
				else {
					v_curr = new Vertex(id_vertices, v_end->pt.x, v_end->pt.y);
					ordered_vertex_chains[i].emplace_back(v_curr);
				}

				Split_Edge_Property edge_property;
				edge_property.width = (*it_h)->e->width;
				edge_property.length = cv::norm(v_curr->pt - v_prev->pt);
				Edge::weight(v_prev->pt, v_curr->pt, edge_property.width, edge_property.region, edge_property.normal, edge_property.grad_center_offset,
					used, I_m, I_t, offset, params->split_params, cdf_grad_m, edge_property.unit_edge_cost, edge_property.mean_angle_diff, edge_property.grad_weight, false);
				split_edges_property[make_pair(v_prev, v_curr)] = edge_property;
			}
			else { // check if it's parallel with the next edge
				assert((*it_h)->e->type == INNER_EDGE && (*std::next(it_h))->e->type == INNER_EDGE);
				Inner_Edge * ie_next = static_cast<Inner_Edge *>((*std::next(it_h))->e);
				Inner_Edge * ie = static_cast<Inner_Edge *>((*it_h)->e);
				Vertex* v_curr = (*it_h)->v1_v2 ? (*it_h)->e->v2 : (*it_h)->e->v1;
				if (v_curr->connectivity()>2 || (ie->get_front_supporting_segment()->node_colinear != ie_next->get_front_supporting_segment()->node_colinear) ) {
					Vertex* v_prev = ordered_vertex_chains[i].back();
					Vertex* v_inner = new Vertex(id_vertices, v_curr->pt.x, v_curr->pt.y);
					added_vertices[v_curr->pt] = v_inner;
					ordered_vertex_chains[i].emplace_back(v_inner);
					Split_Edge_Property edge_property;
					edge_property.width = ie->width;
					edge_property.length = cv::norm(v_inner->pt - v_prev->pt);
					Edge::weight(v_prev->pt, v_inner->pt, edge_property.width, edge_property.region, edge_property.normal, edge_property.grad_center_offset,
						used, I_m, I_t, offset, params->split_params, cdf_grad_m, edge_property.unit_edge_cost, edge_property.mean_angle_diff, edge_property.grad_weight, false);
					split_edges_property[make_pair(v_prev, v_inner)] = edge_property;
				}
			}
			++it_h;
		}
	}

	/* recover gradients used by candidate splitting edges */
	for (auto it = split_edges_property.begin(); it != split_edges_property.end(); ++it) {
		vector<Point2i> & grad_pixels = it->second.region;
		for (size_t k = 0; k < grad_pixels.size(); ++k) used(rows - 1 - grad_pixels[k].y, grad_pixels[k].x) = false;
	}

	// compute energy decrease gained by this split
	double gain_if_merge = 0;
	Label_ID merged_label;
	list<Dual_Merge_Graph::vertex_descriptor> face_vertices_list;
	Dual_Merge_Graph::vertex_iterator vit, vend;
	for (boost::tie(vit, vend) = boost::vertices(sub_merge_graph); vit!=vend; ++vit) {
		face_vertices_list.push_back(*vit);
	}
	subgraph.merge_gain(face_vertices_list, split_edges_map, I, I_m, I_t, gain_if_merge, merged_label, false);

	f->split_gain = gain = -gain_if_merge;

	size_t num_split_edges = 0;
	for (size_t i = 0; i < ordered_vertex_chains.size(); ++i) num_split_edges += ordered_vertex_chains[i].size() - 1;
	assert(num_split_edges = split_edges_property.size());

	// free memory
	Ray::clear_rays(rays);
}


void Partition::split_operator(Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t, Matrix<double> & I_lab_m, Matrix<double> & I_lab_t, 
	vector<Dual_Graph::vertex_descriptor> & graph_vertices, map<int, double> & angles, vector<int> & bin_to_angle_index)
{
	assert(graph_vertices.size() == 1);
	Dual_Graph::vertex_descriptor v = graph_vertices[0];
	// get facet info
	Face* f = dual_graph[v].f;
	int label_f = f->label;
	vector<Point2i>& pixels = f->pixels;

	trace(params->verbose_level, 5, "Split face " + std::to_string(f->id_face) + " (size = " + std::to_string(f->pixels.size()) + ", lab = " + std::to_string(f->label) + ") ");

	// prepare segments in f to be deleted
	for (set<Segment *, costcomp>::iterator it_s = f->detected_segments.begin(); it_s != f->detected_segments.end(); ++it_s) {
		(*it_s)->to_delete = true;
	}

	vector<Face *> f_list;
	f_list = split_facet(f, I, I_prob, used, I_m, I_t, I_lab_m, I_lab_t);
	face_to_graph.erase(f);
	// debug
	for (Face* f_new : f_list) {
		assert(f_new->detected_segments.empty());
	}

	// update the connect graph and priority queues
	Dual_Graph::adjacency_iterator vit, vend;
	boost::tie(vit, vend) = adjacent_vertices(v, dual_graph);

	// add new vertices corresponding to new facets
	for (size_t i=0; i<f_list.size(); ++i) {
		Dual_Graph::vertex_descriptor v_new = boost::add_vertex(dual_graph);
		face_to_graph[f_list[i]] = v_new;
		dual_graph[v_new].f = f_list[i];
		detect_line_segments(I, I_m, I_t, f_list[i], true);
		detect_line_segments(I, I_lab_m, I_lab_t, f_list[i], false);
		extend_detected_segments(f_list[i], I, I_prob, used, I_m, I_t, angles, bin_to_angle_index, dual_graph[v_new].t);
		f_list[i]->delete_invalid_segments();
		queue_split.insert(make_pair(v_new, dual_graph[v_new].t));
	}
	// remove old vertex from the queue
	queue_split.erase(make_pair(v, dual_graph[v].t));

	// add graph edges
	for (size_t i=0; i<f_list.size(); ++i) {
		Face* f_new = f_list[i];
		Dual_Graph::vertex_descriptor v_new = face_to_graph.at(f_new);
		set<Face *> adjacent_faces;
		f_new->get_neighboring_faces(adjacent_faces);

		for (Face* f_adj : adjacent_faces) {
			Dual_Graph::vertex_descriptor v_adj = face_to_graph.at(f_adj);
			bool found_old, found;
			Dual_Graph::edge_descriptor old_connect_e, connect_e;
			boost::tie(old_connect_e, found_old) = boost::edge(v_adj, v, dual_graph);
			if (found_old) queue_merge.erase(make_pair(old_connect_e, dual_graph[old_connect_e].t));

			boost::tie(connect_e, found) = boost::edge(v_adj, v_new, dual_graph);
			if (!found) {
				bool added;
				boost::tie(connect_e, added) = boost::add_edge(v_new, v_adj, dual_graph);
				if (!added) continue;

				merge_gain(f_new, f_adj, I, I_m, I_t, dual_graph[connect_e].t, dual_graph[connect_e].merged_label);
				queue_merge.insert(make_pair(connect_e, dual_graph[connect_e].t));
				// Update other merge queue entries associated with this neighbor face
				if (f_new->label != label_f) {
					Dual_Graph::adjacency_iterator vit_adj, vend_adj;
					boost::tie(vit_adj, vend_adj) = adjacent_vertices(v_adj, dual_graph);
					while (vit_adj != vend_adj) {
						Face *v_f_adj = dual_graph[*vit_adj].f;
						if (*vit_adj != v && *vit_adj != v_new && v_f_adj != nullptr) {
							boost::tie(old_connect_e, found_old) = boost::edge(*vit_adj, v_adj, dual_graph);
							auto erased = queue_merge.erase(make_pair(old_connect_e, dual_graph[old_connect_e].t));
							assert(erased > 0);
							merge_gain(f_adj, v_f_adj, I, I_m, I_t, dual_graph[old_connect_e].t, dual_graph[old_connect_e].merged_label);
							queue_merge.insert(make_pair(old_connect_e, dual_graph[old_connect_e].t));
						}
						++vit_adj;
					}
				}
			}
		}
	}

	// debug
	while (vit != vend) {
		bool connected = false;
		for (size_t i=0; i<f_list.size(); ++i) {
			Dual_Graph::vertex_descriptor v_new = face_to_graph.at(f_list[i]);
			Dual_Graph::edge_descriptor connect_e;
			boost::tie(connect_e, connected) = boost::edge(*vit, v_new, dual_graph);
			if (connected) break;
		}
		assert(connected);
		++vit;
	}

	// Display
	std::string out_msg = "New faces ";
	for (Face* f_new : f_list) {
		out_msg += std::to_string(f_new->id_face) + " (size = " + std::to_string(f_new->pixels.size()) + ", lab = " + std::to_string(f_new->label) + "), ";
	}
	trace(params->verbose_level, 5, out_msg);

	// remove vertex v in the connect graph
	boost::clear_vertex(v, dual_graph);
	boost::remove_vertex(v, dual_graph);
	graph_vertices.clear();
}


void Partition::merge_operator(Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t, Matrix<double> & I_lab_m, Matrix<double> & I_lab_t, 
	vector<Dual_Graph::vertex_descriptor> & vi_vj, map<int, double> & angles, vector<int> & bin_to_angle_index)
{
	// locate vertices for f_i and f_j
	Dual_Graph::vertex_descriptor v_i = vi_vj[0], v_j = vi_vj[1];
	Face* f_i = dual_graph[v_i].f;
	Face* f_j = dual_graph[v_j].f;
	trace(params->verbose_level, 5, "Merge faces " + std::to_string(f_i->id_face) + " (size = " + std::to_string(f_i->pixels.size()) + ", lab = " + std::to_string(f_i->label) + ") and " + std::to_string(f_j->id_face) + " (size = " + std::to_string(f_j->pixels.size()) + ", lab = " + std::to_string(f_j->label) + ")");

	auto find_e_ij = boost::edge(vi_vj[0], vi_vj[1], dual_graph);

	// record adjacency iterators
	Dual_Graph::adjacency_iterator vit_i, vend_i, vit_j, vend_j;
	boost::tie(vit_i, vend_i) = adjacent_vertices(v_i, dual_graph);
	boost::tie(vit_j, vend_j) = adjacent_vertices(v_j, dual_graph);
	// record original face labels
	int label_i = f_i->label, label_j = f_j->label;

	map<int, list<Face *>::iterator> map_facets_to_delete = map<int, list<Face *>::iterator>();
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		if (*it_f == f_i || *it_f == f_j) map_facets_to_delete[(*it_f)->id_face] = it_f;
	}

	// the new facet f_ij
	Face* f_ij = merge_two_facets(map_facets_to_delete, f_i, f_j, I, I_prob, used, I_m, I_t, dual_graph[find_e_ij.first].merged_label);
	//f_ij->find_active_gradients(I_m, I_t, offset, params->split_params, true);
	//f_ij->find_active_gradients(I_lab_m, I_lab_t, offset, params->split_params, false);
	
	// Check inherited segments
	int n_bins = ceil((0.5 - params->split_params.m_value_thresh)/ params->split_params.bin_width);
	int lab_n_bins = ceil((0.5 - params->split_params.lab_m_value_thresh)/ params->split_params.bin_width);
	vector<list<Point2i>> active_gradients_segs(n_bins, list<Point2i>()), active_gradients_lab_segs(lab_n_bins, list<Point2i>());
	
	auto it_s = f_ij->detected_segments.begin();
	while ( it_s != f_ij->detected_segments.end()) {
		Segment* s = *it_s;
		double length_cost = 1 / (1 + exp(5 * (2 * (*it_s)->length / params->split_params.length_thresh - 1))); // = 0.5 when it's half length_thresh
		if ((1 - (*it_s)->support_cost)*(1 - length_cost) >= 0.5) {
			++it_s; continue;
		} 
		else { // reject
			double value_thresh = s->is_color_edge ? params->split_params.m_value_thresh : params->split_params.lab_m_value_thresh;
			int max_id_bin = s->is_color_edge ? f_ij->active_gradients.size() - 1 : f_ij->active_lab_gradients.size() - 1;
			double grad_bin_width = params->split_params.bin_width;
			for (auto it = s->pixels.begin(); it != s->pixels.end(); ++it) {
				int i = rows - 1 - it->y; // rows
				int j = it->x; // columns
				int id_bin = jclamp(0, floor((I_m(i, j) - value_thresh) / grad_bin_width), max_id_bin);
				if (s->is_color_edge) active_gradients_segs[id_bin].push_back(*it);
				else active_gradients_lab_segs[id_bin].push_back(*it);
			}
			it_s = f_ij->detected_segments.erase(it_s);
			delete s;
		}
	}
	
    for (int i = 0; i < n_bins; ++i) {
		set<Point2i, pointcomp> grad_pixs = move(f_ij->active_gradients[i]);
		grad_pixs.insert(make_move_iterator(active_gradients_segs[i].begin()), make_move_iterator(active_gradients_segs[i].end()));
		f_ij->active_gradients[i] = move(grad_pixs);
	}
	for (int i = 0; i < lab_n_bins; ++i) {
		set<Point2i, pointcomp> grad_pixs = move(f_ij->active_lab_gradients[i]);
		grad_pixs.insert(make_move_iterator(active_gradients_lab_segs[i].begin()), make_move_iterator(active_gradients_lab_segs[i].end()));
		f_ij->active_lab_gradients[i] = move(grad_pixs);
	}
	

	// update the connect graph and priority queues
	// add new vertex corresponding to f_ij
	Dual_Graph::vertex_descriptor v_ij = boost::add_vertex(dual_graph);
	dual_graph[v_ij].f = f_ij;
	// detect extra line segments
	detect_line_segments(I, I_m, I_t, f_ij, true);
	detect_line_segments(I, I_lab_m, I_lab_t, f_ij, false);
	extend_detected_segments(f_ij, I, I_prob, used, I_m, I_t, angles, bin_to_angle_index, dual_graph[v_ij].t);
	f_ij->delete_invalid_segments();

	queue_split.insert(make_pair(v_ij, dual_graph[v_ij].t));
	// remove old vertices from the queue
	queue_split.erase(make_pair(v_i, dual_graph[v_i].t));
	queue_split.erase(make_pair(v_j, dual_graph[v_j].t));

	// add edges connecting the new vertex to adjacent vertices
	while (vit_i != vend_i) {
		Dual_Graph::vertex_descriptor v = *vit_i;
		bool found_old;
		Dual_Graph::edge_descriptor old_connect_e;
		boost::tie(old_connect_e, found_old) = boost::edge(v_i, v, dual_graph);
		queue_merge.erase(make_pair(old_connect_e, dual_graph[old_connect_e].t));

		Face *f_adj = dual_graph[v].f;
		if (v != v_j && f_adj != nullptr) {
			Dual_Graph::edge_descriptor connect_e;
			bool inserted;
			
			boost::tie(connect_e, inserted) = (f_ij->pixels.size() > f_adj->pixels.size()) ? boost::add_edge(v_ij, v, dual_graph) : boost::add_edge(v, v_ij, dual_graph);
			merge_gain(f_ij, f_adj, I, I_m, I_t, dual_graph[connect_e].t, dual_graph[connect_e].merged_label);
			if (inserted) queue_merge.insert(make_pair(connect_e, dual_graph[connect_e].t));
			// Update other merge queue entries associated with this neighbor face
			if (f_ij->label != label_i) {
				Dual_Graph::adjacency_iterator vit, vend;
				boost::tie(vit, vend) = adjacent_vertices(v, dual_graph);
				while (vit != vend) {
					Face *v_f_adj = dual_graph[*vit].f;
					if (*vit != v_ij && *vit != v_i && *vit != v_j && v_f_adj != nullptr) {
						boost::tie(old_connect_e, found_old) = boost::edge(*vit, v, dual_graph);
						queue_merge.erase(make_pair(old_connect_e, dual_graph[old_connect_e].t));
						merge_gain(f_adj, v_f_adj, I, I_m, I_t, dual_graph[old_connect_e].t, dual_graph[old_connect_e].merged_label);
						queue_merge.insert(make_pair(old_connect_e, dual_graph[old_connect_e].t));
					}
					++vit;
				}
			}
		}
		++vit_i;
	}
	while (vit_j != vend_j) {
		Dual_Graph::vertex_descriptor v = *vit_j;
		bool found_old;
		Dual_Graph::edge_descriptor old_connect_e;
		boost::tie(old_connect_e, found_old) = boost::edge(v_j, v, dual_graph);
		queue_merge.erase(make_pair(old_connect_e, dual_graph[old_connect_e].t));

		Face *f_adj = dual_graph[v].f;
		if (v != v_i && f_adj != nullptr) {
			Dual_Graph::edge_descriptor connect_e;
			bool found; 
			boost::tie(connect_e, found) = boost::edge(v_ij, v, dual_graph);
			if (!found) {
				bool inserted;
				boost::tie(connect_e, inserted) = (f_ij->pixels.size() > f_adj->pixels.size()) ? boost::add_edge(v_ij, v, dual_graph) : boost::add_edge(v, v_ij, dual_graph);
				merge_gain(f_ij, f_adj, I, I_m, I_t, dual_graph[connect_e].t, dual_graph[connect_e].merged_label);
				if (inserted) queue_merge.insert(make_pair(connect_e, dual_graph[connect_e].t));
				// Update other merge queue entries associated with this neighbor face
				if (f_ij->label != label_j) {
					Dual_Graph::adjacency_iterator vit, vend;
					boost::tie(vit, vend) = adjacent_vertices(v, dual_graph);
					while (vit != vend) {
						Face *v_f_adj = dual_graph[*vit].f;
						if (*vit != v_ij && *vit != v_i && *vit != v_j && v_f_adj != nullptr) {
							boost::tie(old_connect_e, found_old) = boost::edge(*vit, v, dual_graph);
							queue_merge.erase(make_pair(old_connect_e, dual_graph[old_connect_e].t));
							merge_gain(f_adj, v_f_adj, I, I_m, I_t, dual_graph[old_connect_e].t, dual_graph[old_connect_e].merged_label);
							queue_merge.insert(make_pair(old_connect_e, dual_graph[old_connect_e].t));
						}
						++vit;
					}
				}
			}
		}
		++vit_j;
	}
	trace(params->verbose_level, 5, "New faces " + std::to_string(f_ij->id_face) + " (size = " + std::to_string(f_ij->pixels.size()) + ", lab = " + std::to_string(f_ij->label) + ") ");
	// remove vertices in the connect graph
	boost::clear_vertex(v_i, dual_graph);
	boost::clear_vertex(v_j, dual_graph);
	boost::remove_vertex(v_i, dual_graph);
	boost::remove_vertex(v_j, dual_graph);
	vi_vj.clear();

	face_to_graph.erase(f_i);
	face_to_graph.erase(f_j);
	face_to_graph[f_ij] = v_ij;
}


void Partition::reset_indices()
{
	id_vertices = id_edges = id_faces = 0;

	Vertex* v = vertices_head;
	while (v != nullptr) {
		v->id_vertex = id_vertices++;
		v = v->v_next;
	}

	Edge* e = edges_head;
	while (e != nullptr) {
		e->id_edge = id_edges++;
		e = e->e_next;
	}

	for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
		(*it_f)->id_face = id_faces++;
	}
}


Face* Partition::third_facet(Vertex* v_i, HalfEdge* h_ij, Face* f_i, Face* f_j)
{
	for (uint i = 0; i < v_i->directions.size(); ++i) {
		HalfEdge* h_i = v_i->directions[i].second;
		if (h_i == h_ij) {
			continue;
		} else {
			if (h_i->f == f_i || h_i->f == f_j) {
				return h_i->opposite()->f;
			} else {
				return h_i->f;
			}
		}
	}
	return nullptr;
}


Face* Partition::merge_two_facets(map<int, list<Face *>::iterator> & map_facets, Face* f_i, Face* f_j, 
	Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t, Label_ID& label)
{
	// make a list that contains pairs of starting points and lengths in f_i
	list<Edge *> shared_edges;
	vector<HalfEdges_Location> l_ij, l_ji;
	list<pair<int, int>> h_ij, h_ji;
	double shared_length = 0;
	for (int id = 0; id < f_i->edges.size(); ++id) {
		h_ij.clear();
		vector<HalfEdge *> vi_edges(f_i->edges[id].begin(), f_i->edges[id].end());
		uint n_i = f_i->vertices[id].size();
		uint j = 0;
		uint start_i = n_i;
		for (int i = start_i - 1; i >= 0; --i) {
			if ((vi_edges[i]->opposite()->f != nullptr) && (vi_edges[i]->opposite()->f == f_j)) {
				++j;
				shared_edges.emplace_back(vi_edges[i]->e);
				shared_length += vi_edges[i]->e->length;
				if (start_i > i) start_i = i;
			}
			else {
				if (j > 0) h_ij.emplace_back(make_pair(start_i, j));
				j = 0;
			}
		}
		if (j > 0) h_ij.emplace_back(make_pair(0, j));
		// deal with the case where 0 and n-1 are included in different pairs
		if ((h_ij.size() > 1) && (h_ij.front().first + h_ij.front().second == n_i) && (h_ij.back().first == 0)) {
			h_ij.front().second += h_ij.back().second;
			h_ij.pop_back();
		}
		if (h_ij.size() > 0) l_ij.emplace_back(make_pair(id, h_ij));
	}

	// make a list that contains pairs of starting points and lengths in f_j
	for (int id = 0; id < f_j->edges.size(); ++id) {
		h_ji.clear();
		vector<HalfEdge *> vj_edges(f_j->edges[id].begin(), f_j->edges[id].end());
		uint n_j = f_j->vertices[id].size();
		uint j = 0;
		uint start_i = n_j;
		for (int i = start_i - 1; i >= 0; --i) {
			if ((vj_edges[i]->opposite()->f != nullptr) && (vj_edges[i]->opposite()->f == f_i)) {
				++j;
				if (start_i > i) start_i = i;
			}
			else {
				if (j > 0) h_ji.emplace_back(make_pair(start_i, j));
				j = 0;
			}
		}
		if (j > 0) h_ji.emplace_back(make_pair(0, j));
		if ((h_ji.size() > 1) && (h_ji.front().first + h_ji.front().second == n_j) && (h_ji.back().first == 0)) {
			h_ji.front().second += h_ji.back().second;
			h_ji.pop_back();
		}
		if (h_ji.size() > 0) l_ji.emplace_back(make_pair(id, h_ji));
	}
	assert(l_ji.size() > 0 && l_ij.size() > 0);

	typedef pair<Vertex*, Vertex_Type> Vertex_Element;

	int id_i = f_i->id_face;
	int id_j = f_j->id_face;

	// copy features of f_i and f_j
	vector<double> rf_hist_i = f_i->semantic_probabilities, rf_hist_j = f_j->semantic_probabilities; 
	int label_i = f_i->label, label_j = f_j->label;
	set<Segment *, costcomp> split_support_ij(std::move(f_i->detected_segments));
	for (auto it = f_j->detected_segments.begin(); it != f_j->detected_segments.end(); ++it) split_support_ij.insert(*it);

	// copy edges and vertices of f_i and f_j
	vector<list<HalfEdge *>> f_i_edges(f_i->edges.begin(), f_i->edges.end());
	vector<list<HalfEdge *>> f_j_edges(f_j->edges.begin(), f_j->edges.end());
	vector<list<Vertex_Element >> f_i_vertices(f_i->vertices.begin(), f_i->vertices.end());
	vector<list<Vertex_Element >> f_j_vertices(f_j->vertices.begin(), f_j->vertices.end());
	int n_pixels_i = f_i->pixels.size(), n_pixels_j = f_j->pixels.size();
	vector<Point2i> merged_pixels;
	merged_pixels.reserve(n_pixels_i + n_pixels_j);
	merged_pixels.insert(merged_pixels.end(), make_move_iterator(f_i->pixels.begin()), make_move_iterator(f_i->pixels.end()));
	merged_pixels.insert(merged_pixels.end(), make_move_iterator(f_j->pixels.begin()), make_move_iterator(f_j->pixels.end()));
	f_i->pixels = vector<Point2i>();
	f_j->pixels = vector<Point2i>();
	// record aggregation of active gradients and add back pixels around edges to be removed
	double value_thresh = params->split_params.m_value_thresh;
	double lab_value_thresh = params->split_params.lab_m_value_thresh;
	int n_bins = ceil((0.5 - value_thresh)/ params->split_params.bin_width);
	int lab_n_bins = ceil((0.5 - lab_value_thresh)/ params->split_params.bin_width);
	
	vector<set<Point2i, pointcomp>> active_gradients_ij(n_bins, set<Point2i, pointcomp>()), active_lab_gradients_ij(lab_n_bins, set<Point2i, pointcomp>());
	vector<list<Point2i>> active_gradients_edges(n_bins, list<Point2i>());
	for (auto edge_it = shared_edges.begin(); edge_it != shared_edges.end(); ++edge_it) {
		const vector<Point2i> & edge_pixels = (*edge_it)->region;
		for (auto it_p = edge_pixels.begin(); it_p != edge_pixels.end(); ++it_p) {
			int i = rows - 1 - it_p->y; // rows
			int j = it_p->x; // columns
			used(i, j) = false;
			int id_bin = jclamp(0, floor((I_m(i, j) - value_thresh) / params->split_params.bin_width), n_bins - 1);
			active_gradients_edges[id_bin].push_back(*it_p);
		}
	}
	for (int i = 0; i < n_bins; ++i) {
		set<Point2i, pointcomp> grad_pixs;
		if (i < f_i->active_gradients.size()) grad_pixs.insert(make_move_iterator(f_i->active_gradients[i].begin()), make_move_iterator(f_i->active_gradients[i].end()));
		if (i < f_j->active_gradients.size()) grad_pixs.insert(make_move_iterator(f_j->active_gradients[i].begin()), make_move_iterator(f_j->active_gradients[i].end()));
		grad_pixs.insert(make_move_iterator(active_gradients_edges[i].begin()), make_move_iterator(active_gradients_edges[i].end()));
		active_gradients_ij[i] = move(grad_pixs);
	}
	for (int i = 0; i < lab_n_bins; ++i) {
		set<Point2i, pointcomp> grad_pixs;
		if (i < f_i->active_lab_gradients.size()) grad_pixs.insert(make_move_iterator(f_i->active_lab_gradients[i].begin()), make_move_iterator(f_i->active_lab_gradients[i].end()));
		if (i < f_j->active_lab_gradients.size()) grad_pixs.insert(make_move_iterator(f_j->active_lab_gradients[i].begin()), make_move_iterator(f_j->active_lab_gradients[i].end()));
		active_lab_gradients_ij[i] = move(grad_pixs);
	}
	
	vector<list<HalfEdge *>> halfedges_ij;
	// Anticipates the merge of f_i and f_j, by listing the edges of the new facet (f_i U f_j)
	Face::merged_halfedges(f_i, f_j, l_ij, l_ji, halfedges_ij);
	

	// deal with bivalent vertices
	// and record endpoints of overlapping edges of f_i and f_j
	list<Vertex*> endpoints;
	for (int i_ij = 0; i_ij < l_ij.size(); ++i_ij) {
		int id = l_ij[i_ij].first;
		list<pair<int, int>> & h_ij = l_ij[i_ij].second;
		vector<HalfEdge *> edges_i(f_i_edges[id].begin(), f_i_edges[id].end());
		int n_i = edges_i.size();
		for (list<pair<int, int> >::iterator it_h = h_ij.begin(); it_h != h_ij.end(); ++it_h) {
			int seq_start = it_h->first;
			int seq_length = it_h->second;
			HalfEdge* h_1 = edges_i[seq_start], *h_2 = edges_i[(seq_start + seq_length - 1) % n_i];
			Vertex* v_1 = h_1->v1_v2 ? h_1->e->v1 : h_1->e->v2;
			Vertex* v_2 = h_2->v1_v2 ? h_2->e->v2 : h_2->e->v1;
			if (v_1 != v_2) { // if not a closed path
				endpoints.push_back(v_1);
				endpoints.push_back(v_2);
			}
			if (v_1->connectivity() == 3) {
				Face* f_1 = third_facet(v_1, h_1, f_i, f_j);
				if (f_1 != nullptr) f_1->set_as_bivalent(v_1);
			}
			if (v_2->connectivity() == 3) {
				Face* f_2 = third_facet(v_2, h_2, f_i, f_j);
				if (f_2 != nullptr) f_2->set_as_bivalent(v_2);
			}
		}
	}
	
	// Erases the old facets f_i and f_j 
	faces.erase(map_facets[id_i]);
	faces.erase(map_facets[id_j]);
	delete f_j;
	delete f_i;
	for (Edge* e : shared_edges) erase(e, true);

	list<Point2d> pts_to_check;
	for (int i_ij = 0; i_ij < l_ij.size(); ++i_ij) {
		int id = l_ij[i_ij].first;
		list<pair<int, int>> & h_ij = l_ij[i_ij].second;
		vector<Vertex_Element> vertices_i(f_i_vertices[id].begin(), f_i_vertices[id].end());
		int n_i = vertices_i.size();
		for (list<pair<int, int> >::iterator it_h = h_ij.begin(); it_h != h_ij.end(); ++it_h) {
			// Erases unconnected vertices along the deleted shared edges between v_1 and v_2
			int seq_start = it_h->first;
			int seq_length = it_h->second;
			if (seq_length <= 1) continue;
			for (uint i = 0; i < seq_length; ++i) {
				Vertex*& v_i = vertices_i[(seq_start + i) % n_i].first;
				if (v_i != nullptr && v_i->connectivity()==0) {
					queue_bivalent.erase(pair<Vertex*, double>(v_i,v_i->energy_gain));
					bboxes_by_x.erase(pair<Vertex*, double>(v_i, v_i->bbox.first));
					bboxes_by_y.erase(pair<Vertex*, double>(v_i, v_i->bbox.second));
					
					pts_to_check.push_back(v_i->pt);
					erase(v_i);
					quadtree->remove(v_i, true);
					v_i = nullptr;
				}
			}
		}
	}

	// Builds the new facet
	Face* f_ij = new Face(id_faces, halfedges_ij);
	f_ij->label = label;

	set<Vertex*> vertices_checked;
	// Merge collinear edges that passes through endpoints of overlapping edges of f_i and f_j, THEN pixels inside adjacent facets need to be updated, also the intersection points of detected line segments
	bool area_changed = false;
	list<Vertex*> vertices_erased;
	for (list<Vertex*>::iterator it_v = endpoints.begin(); it_v != endpoints.end(); ++it_v) {
		Vertex* v_1 = *it_v;
		Edge* me = nullptr;
		me = merge_collinear_edges(v_1, vertices_erased, false);
		if (me != nullptr) {
			trace(params->verbose_level, 5, "remove bivalent vertices between vertices " 
			+ std::to_string(me->v1->id_vertex) + " and " + std::to_string(me->v2->id_vertex));
			me->weight(used, I_m, I_t, offset, params->split_params, cdf_grad_m);
		}
	}
	for (Vertex* v: vertices_erased) delete v;

	//// check proximity of removed vertices
	/*for (Point2d pt : pts_to_check) {
		list<Vertex *> neighbors;
		if (!bboxes_by_x.empty()) {
			double x_min = pt.x - bboxes_by_x.begin()->second;
			double x_max = pt.x + bboxes_by_x.begin()->second;
			double y_min = pt.y - bboxes_by_y.begin()->second;
			double y_max = pt.y + bboxes_by_y.begin()->second;
			quadtree->search(x_min, x_max, y_min, y_max, neighbors);
			for (Vertex* v_neighbor : neighbors) {
				if (vertices_checked.find(v_neighbor) != vertices_checked.end()) continue;
				check_removal_criteria(v_neighbor, I_m, I_t, I, I_prob, used);
				vertices_checked.emplace(v_neighbor);
			}
		}
	}*/
	for (int id = 0; id < f_ij->vertices.size(); ++id) {
		for (pair<Vertex *, Vertex_Type> v_ij : f_ij->vertices[id]) {
			if (vertices_checked.find(v_ij.first) != vertices_checked.end()) continue;
			check_removal_criteria(v_ij.first, I_m, I_t, I, I_prob, used);
		}
	}

	// assign pixels to f_ij
	if (area_changed) f_ij->find_pixels_inside_facet();
	else f_ij->pixels = std::move(merged_pixels);
	f_ij->active_gradients = std::move(active_gradients_ij);
	f_ij->active_lab_gradients = std::move(active_lab_gradients_ij);
	f_ij->detected_segments = std::move(split_support_ij);
	Face::compute_merged_histogram({rf_hist_i, rf_hist_j}, f_ij->semantic_probabilities);

	/* check if new facet is a candidate for splitting: 
	if the amount of label inconsistency has increased
	*/
	//if (f_ij->pixels.size() > 0.15*area) {
	//	for (auto s : f_ij->detected_segments) delete s;
	//	f_ij->detected_segments.clear();
	//	f_ij->check_split = false;
	//}

	faces.push_back(f_ij);
	return f_ij;
}


bool Partition::loop_inside_sub_facet(HalfEdge * h_e, vector<list<HalfEdge*>> & sub_facet, vector<list<HalfEdge*>> & parent_facet)
{
	// There are two cases: 
	// Case 1: h_e touches the polygon contour.
	// Case 2: h_e does not touch any polygon contour. In this case we will trace a self-closed loop (either clockwise or anti-clockwise)

	Vertex* v1 = (h_e->v1_v2 ? h_e->e->v1 : h_e->e->v2);
	Vertex* v2 = (h_e->v1_v2 ? h_e->e->v2 : h_e->e->v1);

	sub_facet.resize(1, list<HalfEdge*>());

	Face* f = parent_facet[0].front()->f;

	// first list all inner facets of f
	map<Face *, int> f_inner;
	for (int i = 1; i < parent_facet.size(); ++i) {
		// loop through the inner contour to find all inner facets inside this contour
		for (list<HalfEdge *>::iterator it_ie = parent_facet[i].begin(); it_ie != parent_facet[i].end(); ++it_ie) {
			Face *inner_facet = (*it_ie)->opposite()->f;
			f_inner[inner_facet] = i;
		}
	}

	// Construct the contour that includes h_e, record indices of inner contours of f visited
	set<int> contours_visited;
	list<HalfEdge*> contour_edges;
	
	HalfEdge* h_f = h_e;
	do {
		contour_edges.push_back(h_f);
		// Finds the next half-edge by turning right
		// Adds it to the contour we are currently defining
		Vertex* v = (h_f->v1_v2 ? h_f->e->v2 : h_f->e->v1);
		uint v_n = uint(v->directions.size());
		int path = 0;
		for (uint i = 0; i < v_n; i++) {
			if (v->directions[i].second->e == h_f->e) {
				path = (i + 1) % v_n;
				break;
			}
		}
		h_f = v->directions[path].second;
		
		// Check if h_f belongs to an inner contour of f.
		if (f_inner.count(h_f->opposite()->f) != 0) {
			contours_visited.insert(f_inner.at(h_f->opposite()->f));
		}
		else if (h_f->f == f) { 
			// Continue to check if h_f belongs to the outer contour of f, provided that h_f is not on an inner contour.
			contours_visited.insert(0);
		}	

	} while (h_f != h_e);

	bool is_first_loop_clockwise = true;
	if (contours_visited.count(0)) {
		// If the outer contour is visited, the contour must be clockwise and form the outer contour of the resulting facet.
		sub_facet[0] = contour_edges;
	}
	else {
		// check if the first contour is anti-clockwise.
		int clockwise = HalfEdge::is_clockwise(contour_edges);
		if (clockwise == -1) { // area=0
			return false;
		}
		if (clockwise == 1) {
			sub_facet[0] = contour_edges;
		}
		else {
			is_first_loop_clockwise = false;
			sub_facet.reserve(2);
			sub_facet[0] = parent_facet[0];
			sub_facet.push_back(contour_edges);
		}
	}

	// Now add unvisited inner contours that belongs to this subfacet
	for (int id = 1; id < parent_facet.size(); ++id) {
		if (contours_visited.count(id) != 0) continue;
		// select any inner point of the inner contour (but not a vertex because a vertex can lie also on the outer boundary)
		Point2d pt = 0.5*(parent_facet[id].back()->e->v1->pt + parent_facet[id].back()->e->v2->pt);
		// check if the selected vertex lies inside the outer contour
		bool is_inside;
		if (is_first_loop_clockwise) {
			is_inside = Geometry::is_inside_contour(pt, sub_facet[0]);
		}
		else {
			is_inside = (Geometry::is_inside_contour(pt, sub_facet[0]) && (!Geometry::is_inside_contour(pt, sub_facet[1])));
		}
		if (is_inside) {
			sub_facet.push_back(parent_facet[id]);
		}
	}

	return true;
}



vector<Face *> Partition::split_facet(Face* f, Matrix<uchar> & I, Matrix<double> & I_prob, Matrix<bool> & used, Matrix<double> & I_m, Matrix<double> & I_t, Matrix<double> & I_lab_m, Matrix<double> & I_lab_t)
{
	vector<Face *> f_list;

	vector<vector<Vertex *>>& split_vertex_chains = f->split_vertex_chains;
	int n_chains = split_vertex_chains.size();

	vector<list<pair<Vertex *, Vertex_Type> >> f_vertices = f->vertices;  // for debug only
	map<Vertex*, pair<int, int>> split_vertex_outer_edge_location = f->split_vertex_outer_edge_location;
	Label_ID init_label = f->label;
	map<pair<Vertex *, Vertex *>, Split_Edge_Property> split_edges_property = f->split_edges_property;
	int n_split_edges = split_edges_property.size();
	int m_rows = static_cast<int>(I_m.rows);
	// debug
	vector<Point2i> set_pixels_f = f->pixels;

	// Add new vertices
	set<Vertex *> vertex_added;
	set<Vertex*> vertices_to_check;
	for (int i = 0; i < n_chains; ++i) {
		int n_chain_vertices = split_vertex_chains[i].size();
		for (int j = 0; j < n_chain_vertices; ++j) {
			Vertex* v = split_vertex_chains[i][j];
			vertices_to_check.insert(v);
			if (j == 0 || j == n_chain_vertices - 1) {
				if (vertex_added.count(v) == 0 && v->directions.size() == 0) {
					quadtree->add(v);
					push_back(v); 
					vertex_added.insert(v);
				}
			}
			else {
				quadtree->add(v);
				push_back(v);
				vertex_added.insert(v);
			}
		}
	}
	
	// Create new edges
	set<Edge *> new_edges;
	for (int i = 0; i < n_chains; ++i) {
		int n_chain_vertices = split_vertex_chains[i].size();
		for (int j = 0; j < n_chain_vertices; ++j) {
			Vertex* v = split_vertex_chains[i][j];
			Vertex* v_next = split_vertex_chains[i][(j + 1) % n_chain_vertices];
			
			for (int path = 0; path < v->directions.size(); ++path) {
				HalfEdge* h_i = v->directions[path].second;
				assert(h_i->e->v1 != v_next || h_i->e->v2 != v_next);
			}

			// add edge
			map<pair<Vertex *, Vertex *>, Split_Edge_Property>::iterator it_e = f->split_edges_property.find(pair<Vertex *, Vertex *>(v, v_next));
			if (it_e != f->split_edges_property.end()) {
				Split_Edge_Property & e_property = it_e->second;
				double alpha_1 = atan2(v_next->pt.y - v->pt.y, v_next->pt.x - v->pt.x);
				double alpha_2 = (alpha_1 <= 0 ? alpha_1 + PI : alpha_1 - PI);
				double e_a = -sin(alpha_1), e_b = cos(alpha_1), e_c = -e_a * v_next->pt.x - e_b * v_next->pt.y;
				set<Ray*> support_rays;
				Edge* e = new Inner_Edge(id_edges, support_rays, 0, v, v_next, e_property.width, alpha_1, alpha_2);
				push_back(e);
				new_edges.insert(e);
				// set edge properties
				e->unit_edge_cost = e_property.unit_edge_cost;
				e->mean_angle_diff = e_property.mean_angle_diff;
				e->grad_center_offset = e_property.grad_center_offset;
				e->grad_weight = e_property.grad_weight;
				e->normal = e_property.normal;
				e->region = e_property.region;
				// remove from the candidate list
				f->split_edges_property.erase(it_e);
				/* update used gradients */
				const vector<Point2i> & edge_region = e->region;
				for (auto it_p = edge_region.begin(); it_p != edge_region.end(); ++it_p) {
					assert((it_p->y >= 0) && (it_p->y < rows) || (it_p->x >= 0) || (it_p->x < cols));
					used(m_rows - 1 - it_p->y, it_p->x) = true;
				}
			}
		}
	}
	assert(f->split_edges_property.size() == 0);

	// Splits facet edges if intersected (update split_vertex_outer_edge_location inside f after each split)
	map<Vertex*, pair<int, int>>::iterator it_v = f->split_vertex_outer_edge_location.begin();
	while (it_v != f->split_vertex_outer_edge_location.end()) {
		map<Vertex*, pair<int, int>>::iterator it_v_next = it_v;
		++it_v_next;

		int contour_id = it_v->second.first;
		int halfedge_id = it_v->second.second;
		list<HalfEdge *> & halfedges = f->edges[contour_id];
		list<HalfEdge *>::iterator it_h1 = halfedges.begin();
		std::advance(it_h1, halfedge_id);
		HalfEdge* h1 = *it_h1;
		/* update used gradients */
		const vector<Point2i> & edge_region = h1->e->region;
		for (auto it_p = edge_region.begin(); it_p != edge_region.end(); ++it_p) {
			used(m_rows - 1 - it_p->y, it_p->x) = false;
		}
		Vertex* A = (h1->v1_v2 ? h1->e->v1 : h1->e->v2);
		Vertex* B = (h1->v1_v2 ? h1->e->v2 : h1->e->v1);
		//std::cout << A->id_vertex << "--" << it_v->first->id_vertex << "--" << B->id_vertex << std::endl;
		// Splits the intersected halfedge, update adjacent facet and f
		if (it_v->first!=A && it_v->first != B) split_edge_and_update_facets(h1, it_v->first, true);  // new edges: (A v1), (v1 B)
		it_v = it_v_next;
	}

	// build new facets
	vector<list<HalfEdge *>> f_edges = f->edges;
	build_split_faces(new_edges, f, f_list);

	size_t sum_pixels = 0;
	for (Face* f_new : f_list) {
		f_new->find_pixels_inside_facet(offset);
		sum_pixels += f_new->pixels.size();
		// compute edge scores
		for (size_t id = 0; id < f_new->edges.size(); ++id) {
			for (list<HalfEdge *>::iterator it_h = f_new->edges[id].begin(); it_h != f_new->edges[id].end(); ++it_h) {
				if ((*it_h)->e->unit_edge_cost < 0) (*it_h)->e->weight(used, I_m, I_t, offset, params->split_params, cdf_grad_m);
			}
		}
		f_new->find_active_gradients(I_m, I_t, offset, params->split_params, true); // call after weighting edges
		f_new->find_active_gradients(I_lab_m, I_lab_t, offset, params->split_params, false);
		f_new->compute_feature(I_prob, I, offset, params->max_label);
		f_new->estimate_label(0, area, diag_length);

		if (f_new->label != init_label) {
			for (int id = 0; id < f_new->vertices.size(); ++id) {
				for (pair<Vertex*, Vertex_Type> v_pair : f_new->vertices[id]) vertices_to_check.insert(v_pair.first);
			}
		}
	}
	
	//for (Vertex* v : vertices_to_check) {
	//	list<Vertex*> neighbors;
	//	if (!bboxes_by_x.empty()) {
	//		double x_min = v->pt.x - bboxes_by_x.begin()->second;;
	//		double x_max = v->pt.x + bboxes_by_x.begin()->second;;
	//		double y_min = v->pt.y - bboxes_by_y.begin()->second;;
	//		double y_max = v->pt.y + bboxes_by_y.begin()->second;;
	//		quadtree->search(x_min, x_max, y_min, y_max, neighbors);
	//		vertices_to_check.insert(neighbors.begin(), neighbors.end());
	//	}
	//}
	auto it_bv = bboxes_by_x.begin();
	while (it_bv != bboxes_by_x.end()) {
		auto it_bv_next = std::next(it_bv);
		check_removal_criteria(it_bv->first, I_m, I_t, I, I_prob, used);
		it_bv = it_bv_next;
	}
	

	for (Vertex* v : vertices_to_check) check_removal_criteria(v, I_m, I_t, I, I_prob, used);

	assert(std::abs(int(sum_pixels) - int(set_pixels_f.size())) <= 20);
	
	return f_list;
}



void Partition::split_edge_and_update_facets(HalfEdge* h, Vertex* v, bool destroy)
{
	// We want to access the opposite facet of h which is going to take one more vertex in its definition : v.
	Face* f = h->f, *f_adj = h->opposite()->f;

	bool v1_v2 = h->v1_v2;
	Vertex* v1 = (v1_v2 ? h->e->v1 : h->e->v2);
	Vertex* v2 = (v1_v2 ? h->e->v2 : h->e->v1);

	double alpha_v1 = v1->incidence_angle(h->e), alpha_v2 = v2->incidence_angle(h->e);
	Edge *v1_v = nullptr, *v_v2 = nullptr;

	if (h->e->type == INNER_EDGE) {
		Inner_Edge* e = static_cast<Inner_Edge *>(h->e);
		set<Ray *> support;
		e->get_supporting_rays(support);
		v1_v = new Inner_Edge(id_edges, support, 3, v1, v, e->width, alpha_v1, alpha_v2);
		v_v2 = new Inner_Edge(id_edges, support, 3, v, v2, e->width, alpha_v1, alpha_v2);
	} else if (h->e->type == OUTER_EDGE) {
		Outer_Edge* e = static_cast<Outer_Edge *>(h->e);
		Image_Boundary boundary = e->boundary;
		v1_v = new Outer_Edge(id_edges, boundary, v1, v, e->width, alpha_v1, alpha_v2, true);
		v_v2 = new Outer_Edge(id_edges, boundary, v, v2, e->width, alpha_v1, alpha_v2, true);
	} else if (h->e->type == ARTIFICIAL_EDGE) {
		Artificial_Edge* e = static_cast<Artificial_Edge *>(h->e);
		double ae_a, ae_b, ae_c;
		e->line_equation(ae_a, ae_b, ae_c);
		v1_v = new Artificial_Edge(id_edges, ae_a, ae_b, ae_c, v1, v, e->width, alpha_v1, alpha_v2);
		v_v2 = new Artificial_Edge(id_edges, ae_a, ae_b, ae_c, v, v2, e->width, alpha_v1, alpha_v2);
	}

	if (f_adj != nullptr) {
		// The order or the vertices in the adjacent facet is (v2 v v1)
		f_adj->add_non_corner(h->opposite(), v, v_v2->v2_v1, v1_v->v2_v1, params->prop_min_edge);
	}
	if (f != nullptr) {
		// The order or the vertices in the facet is (v1 v v2)
		f->add_non_corner(h, v, v1_v->v1_v2, v_v2->v1_v2, params->prop_min_edge);
	}
	push_back(v1_v);
	push_back(v_v2);
	erase(h->e, destroy);
}


#if 0
void Partition::fill_matrix_of_iterators()
{
    for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end() ; it_f++) {
        Face* f = (*it_f);
        for (list<Point2i>::iterator it_p = f->pixels.begin() ; it_p != f->pixels.end() ; it_p++) {
            F(I.rows - 1 - it_p->y, it_p->x) = it_f;
        }
    }
}
#endif


int Partition::count_nonsimple_faces()
{
    int count = 0;
    for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
        if (!(*it_f)->is_simple(false)) {
            ++count;
        }
    }
    return count;
}


void Partition::draw_intermediate_graph(Matrix<uchar> & I, Matrix<uchar> & J, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
    vector<list<Edge *> > & outer_edges, vector<list<Edge *> > & inner_edges)
{
	Point2d shift(0.5, 0.5);
    // Copies an enlarged version of the image
    uint f = 4;
    J = Matrix<uchar>(f * I.rows, f * I.cols, 3);
    for (uint i = 0 ; i < I.rows ; i++) {
        for (uint j = 0 ; j < I.cols ; j++) {
            for (uint k = 0 ; k < f ; k++) {
                for (uint l = 0 ; l < f ; l++) {
                    for (uint c = 0 ; c < 3 ; c++) J(f * i + k, f * j + l, c) = I(i, j, c);
                }
            }
        }
    }

    uchar blue[3] = {0, 0, 255};
    uchar green[3] = {0, 255, 0};
    uchar red[3] = {255, 0, 0};
    uchar yellow[3] = {255, 255, 0};

    // Outputs edges
    for (uint i = 0; i < outer_edges.size(); i++) {
        for (list<Edge *>::iterator it = outer_edges[i].begin(); it != outer_edges[i].end(); it++) {
            Edge* e = *it;
            Point2d pt1 = e->v1->pt + shift;
            Point2d pt2 = e->v2->pt + shift;
            Point2i pi1 = Point2i(int(round(jclamp(0, f * (pt1.x), J.cols - 1))), int(round(jclamp(0, J.rows-1 - f * pt1.y, J.rows - 1))));
            Point2i pi2 = Point2i(int(round(jclamp(0, f * (pt2.x), J.cols - 1))), int(round(jclamp(0, J.rows-1 - f * pt2.y, J.rows - 1))));
            J.line(pi1.y, pi1.x, pi2.y, pi2.x, red);
        }
    }
    for (uint i = 0; i < inner_edges.size(); i++) {
        for (list<Edge *>::iterator it = inner_edges[i].begin(); it != inner_edges[i].end(); it++) {
            Edge* e = *it;
            Point2d pt1 = e->v1->pt + shift;
            Point2d pt2 = e->v2->pt + shift;
			Point2i pi1 = Point2i(int(round(jclamp(0, f * pt1.x, J.cols - 1))), int(round(jclamp(0, J.rows-1 - f * pt1.y, J.rows - 1))));
			Point2i pi2 = Point2i(int(round(jclamp(0, f * pt2.x, J.cols - 1))), int(round(jclamp(0, J.rows-1 - f * pt2.y, J.rows - 1))));
            J.line(pi1.y, pi1.x, pi2.y, pi2.x, red);
        }
    }

    // Draws vertices
    for (int l = 0 ; l < 2 ; l++) {
        list<Vertex *> & reference_to_list = (l == 0 ? inner_vertices : outer_vertices);
        for (list<Vertex *>::iterator it = reference_to_list.begin(); it != reference_to_list.end(); it++) {
            Vertex* v = *it;
            Point2d pt = v->pt + shift;
            Point2i pi = Point2i(int(round(jclamp(0, f * pt.x, J.cols - 1))), int(round(jclamp(0, J.rows-1 - f * pt.y, J.rows - 1))));
            if (v->connectivity() == 2) {
                for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = green[c];
            } else if (v->connectivity() == 3) {
                for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = blue[c];
            } else {
                for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = yellow[c];
            }
        }
    }
}


void Partition::draw_highlighted_facets(Matrix<uchar> & J)
{
	bool is_argb = (J.channels == 4);

	uchar white[3] = { 255, 255, 255 };
	uchar red[3] = { 255, 0, 0 };
	uchar green[3] = { 0, 255, 0 };

	Matrix<int> B = Matrix<int>(rows, cols, 1);
	for (uint i = 0; i < rows; i++) {
		for (uint j = 0; j < cols; j++) {
			B(i, j) = 0;
		}
	}

	// Label pixels
	for (int i = 0; i < faces_next_operation.size(); ++i) {
		Face* f = faces_next_operation[i];
		for (auto it_p = f->pixels.begin(); it_p != f->pixels.end(); it_p++) {
			uint i = uint(int(B.rows)-1 - it_p->y);
			uint j = uint(it_p->x);
			if (B(i, j) != 0) {
				std::cout << "faces: " << B(i, j) << " " << f->id_face << std::endl;
				std::cout << "pixel: " << i << " " << j << std::endl;
				std::cout << "B(i, j) = " << B(i, j) << std::endl;
			}
			assert(B(i, j) == 0);
			B(i, j) = 1;
		}
	}	

	// Colors pixels
	for (uint i = 0; i < B.rows; i++) {
		for (uint j = 0; j < B.cols; j++) {
			if (B(i, j) == 1) {
				for (int c = 0; c < 3; c++) J(i, j, c) = uchar((int)white[c]);
				if (is_argb) J(i, j, 3) = 255;
			} 
			else if (B(i, j) == 2) {
				for (int c = 0; c < 3; c++) J(i, j, c) = uchar((int)red[c]);
				if (is_argb) J(i, j, 3) = 255;
			}
			else if (B(i, j) == 3) {
				for (int c = 0; c < 3; c++) J(i, j, c) = uchar((int)green[c]);
				if (is_argb) J(i, j, 3) = 255;
			}
		}
	}

}


void Partition::draw_active_gradients(Matrix<uchar> & J, double alpha)
{
	uchar green[3] = { 0, 255, 0 };
	uchar red[3] = { 255, 0, 0 };
	uchar cyan[3] = { 0, 255, 255 };
	uchar yellow[3] = { 255, 255, 0 };

	bool is_argb = (J.channels == 4);

	Matrix<int> B = Matrix<int>(rows, cols, 1);
	for (uint i = 0; i < rows; i++) {
		for (uint j = 0; j < cols; j++) {
			B(i, j) = 0;
		}
	}

	// Label pixels
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		Face* f = (*it_f);
		for (int id = 0; id < f->active_gradients.size(); ++id) {
			for (set<Point2i, pointcomp>::iterator it_p = f->active_gradients[id].begin(); it_p != f->active_gradients[id].end(); it_p++) {
				uint i = uint(int(B.rows)-1 - it_p->y);
				uint j = uint(it_p->x);
				if (B(i, j) != 0) {
					std::cout << "faces: " << B(i, j) << " " << f->id_face << std::endl;
					std::cout << "pixel: " << i << " " << j << std::endl;
					std::cout << "B(i, j) = " << B(i, j) << std::endl;
				}
				assert(B(i, j) == 0);
				B(i, j) = 1;
			}
		}
		for (int id = 0; id < f->active_lab_gradients.size(); ++id) {
			for (set<Point2i, pointcomp>::iterator it_p = f->active_lab_gradients[id].begin(); it_p != f->active_lab_gradients[id].end(); it_p++) {
				uint i = uint(int(B.rows) - 1 - it_p->y);
				uint j = uint(it_p->x);
				B(i, j) = 2;
			}
		}
	}

	// Colors pixels
	for (uint i = 0; i < B.rows; i++) {
		for (uint j = 0; j < B.cols; j++) {
			if (B(i, j) == 1) {
				for (int c = 0; c < 3; c++) J(i, j, c) = uchar((1-alpha)*(int)J(i, j, c) + alpha*(int)cyan[c]);
				if (is_argb) J(i, j, 3) = 255;
			}
			else if (B(i, j) == 2) {
				 for (int c = 0; c < 3; c++) J(i, j, c) = uchar((1 - alpha)*(int)J(i, j, c) + alpha*(int)red[c]);
				 if (is_argb) J(i, j, 3) = 255;
			}
			else if (B(i, j) == 3) {
				for (int c = 0; c < 3; c++) J(i, j, c) = uchar((1 - alpha)*(int)J(i, j, c) + alpha*(int)yellow[c]);
				if (is_argb) J(i, j, 3) = 255;
			}
			else if (B(i, j) == 4) {
				for (int c = 0; c < 3; c++) J(i, j, c) = uchar((1 - alpha)*(int)J(i, j, c) + alpha*(int)green[c]);
				if (is_argb) J(i, j, 3) = 255;
			}
			else if (B(i, j) == 5) {
				for (int c = 0; c < 3; c++) J(i, j, c) = uchar((1 - alpha)*(int)J(i, j, c) + alpha*(int)red[c]);
				if (is_argb) J(i, j, 3) = 255;
			}
		}
	}

}

void Partition::draw_edges(Matrix<uchar> & I, Matrix<uchar> & J, double f)
{
	Point2d shift(0.5, 0.5);
	if (f == 1.0) {
		J = I;
	} else {
		J = Matrix<uchar>(I, f);
	}

	uchar neon[3] = { 57, 255, 20 };
	uchar green[3] = { 0, 255, 0 };
	uchar blue[3] = { 0, 0, 255 };
    uchar dark_green[3] = {0, 128, 0};
    uchar red[3] = {255, 0, 0};
	uchar dark_orange[3] = { 255, 147, 0 };
    uchar yellow[3] = {255, 255, 0};

	double alpha = 0.3;
	// highlight active gradient pixels
	//draw_active_gradients(J, alpha);

	// draw edges
    Edge* e = edges_head;
    while (e != NULL) {
		Point2d pt1 = e->v1->pt + shift;
		Point2d pt2 = e->v2->pt + shift;
		Point2i pi1 = Point2i(int(round(jclamp(0, f * pt1.x, J.cols - 1))), int(round(jclamp(0, J.rows-1 - f * pt1.y, J.rows - 1))));
		Point2i pi2 = Point2i(int(round(jclamp(0, f * pt2.x, J.cols - 1))), int(round(jclamp(0, J.rows-1 - f * pt2.y, J.rows - 1))));
		J.line(pi1.y, pi1.x, pi2.y, pi2.x, red);
        e = e->e_next;
    }

    Vertex* v = vertices_head;
    while (v != NULL) {
        Point2d pt = v->pt + shift;
        Point2i pi = Point2i(int(round(jclamp(0, f * pt.x, J.cols - 1))), int(round(jclamp(0, J.rows-1 - f * pt.y, J.rows - 1))));
        if (v->connectivity() == 2) {
            for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = yellow[c];
        } else if (v->connectivity() == 3) {
            for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = blue[c];
        } else {
            for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = dark_green[c];
        }
        v = v->v_next;
    }
}

void Partition::draw_labels(Matrix<uchar> & J, const map<Label_ID, int> & label_id_to_class, bool color)
{
	Matrix<uchar> B = Matrix<uchar>(rows, cols, 1, INVALID_LABEL);
	// Assigns a label index to all pixels
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		Face* f = (*it_f);
		if (f->pixels.empty()) f->find_pixels_inside_facet(offset);
		int k = label_id_to_class.empty() ? f->label : label_id_to_class.at(f->label);
		for (auto it_p = f->pixels.begin(); it_p != f->pixels.end(); it_p++) {
			uint i = uint(int(B.rows) - 1 - it_p->y);
			uint j = uint(it_p->x);
			B(i, j) = k;
		}
	}
	if (color) {
		// Prepares palette
		map<int, Vec3b> colors;
		std::uniform_int_distribution<int> uniform_dist(100, 255);
		for (int l = 0; l <= params->max_label; l++) {
			int r = uniform_dist(generator);
			int g = uniform_dist(generator);
			int b = uniform_dist(generator);
			colors.insert(std::make_pair(l, Vec3b(r, g, b)));
		}
		// Colors facets
		J = Matrix<uchar>(rows, cols, 3);
		for (uint i = 0; i < B.rows; i++) {
			for (uint j = 0; j < B.cols; j++) {
				int k = B(i, j);
				for (int c = 0; c < 3; c++) J(i, j, c) = uchar(colors[k][c]);
			}
		}
	}
	else {
		J = B;
	}
}

void Partition::draw_faces(Matrix<uchar> & J)
{
    Matrix<int> B = Matrix<int>(rows, cols, 1);
    for (uint i = 0 ; i < rows ; i++) {
        for (uint j = 0 ; j < cols ; j++) {
            B(i, j) = -1;
        }
    }

    // Assigns a facet index to all pixels
    for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
        Face* f = (*it_f);
		if (f->pixels.empty()) f->find_pixels_inside_facet(offset);
        int k = f->id_face;
		for (auto it_p = f->pixels.begin(); it_p != f->pixels.end(); it_p++) {
            uint i = uint(int(B.rows)-1 - it_p->y);
            uint j = uint(it_p->x);
			if (B(i, j) != -1) {
				std::cout << "faces: " << B(i, j) << " " << k << std::endl;
				std::cout << "pixel: " << i << " " << j << std::endl;
				std::cout << "B(i, j) = " << B(i, j) << std::endl;
			}
            //assert(B(i, j) == -1);
            B(i, j) = k;
        }
    }

    // Prepares palette
    map<int, Vec3b> colors;
	std::uniform_int_distribution<int> uniform_dist (100, 255);
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		int r = uniform_dist(generator);
		int g = uniform_dist(generator);
		int b = uniform_dist(generator);
		colors.insert(std::make_pair((*it_f)->id_face, Vec3b(r, g, b) ));
	}
    // Colors facets
    J = Matrix<uchar>(rows, cols, 3);
    for (uint i = 0; i < B.rows; i++) {
        for (uint j = 0; j < B.cols; j++) {
            int k = B(i, j);
            for (int c = 0 ; c < 3 ; c++) J(i, j, c) = uchar(colors[k][c]);
        }
    }
}


void Partition::save_edges(Matrix<uchar> & I, const std::string &directory, const std::string &basename, double f)
{
    Matrix<uchar> J;
    draw_edges(I, J, f);
	std::string filename;

    filename += directory + separator + basename + "_edges" + (params->suffix == "" ? "" : "_" + params->suffix) + ".png";

    J.write_uchar(filename);
}


void Partition::save_labels(const std::string &directory, const std::string &basename, const map<Label_ID, int> & label_id_to_class)
{
	Matrix<uchar> J;
	draw_labels(J, label_id_to_class, false);
	std::string filename = directory + separator + basename + "_labels" + (params->suffix == "" ? "" : "_" + params->suffix) + ".tiff";
	J.write_uchar(filename);
}

void Partition::save_faces(const std::string &directory, const std::string &basename)
{
    Matrix<uchar> J_1;
	draw_faces(J_1);
	std::string filename = directory + separator + basename + "_faces" + (params->suffix == "" ? "" : "_" + params->suffix) + ".tiff";
	J_1.write_uchar(filename);
}


void Partition::draw_faces_svg(const std::string & directory, const std::string & fileNameWithoutExtension, list<Face *> f_list) {
	Point2d shift(0.5, 0.5);

	std::string filename = directory + separator + fileNameWithoutExtension + ".svg";

	// Prepares palette
	map<int, Vec3b> colors;
	std::uniform_int_distribution<int> uniform_dist(100, 255);
	for (list<Face *>::iterator it_f = f_list.begin(); it_f != f_list.end(); it_f++) {
		if (*it_f == nullptr) continue;
		int r = uniform_dist(generator);
		int g = uniform_dist(generator);
		int b = uniform_dist(generator);
		colors.insert(std::make_pair((*it_f)->id_face, Vec3b(r, g, b)));
	}

	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);
	SVG::markup_header(os, int(rows), int(cols));
	for (Face * f : f_list) {
		if (f == nullptr) continue;
		// polygons and vertices
		SVG::markup_polygon(os, int(rows), int(cols), f, colors[f->id_face][0], colors[f->id_face][1], colors[f->id_face][2], 0.1, 0.4, shift.x, shift.y);
		for (int id = 0; id < f->vertices.size(); ++id) {
			vector<pair<Vertex *, Vertex_Type>> contour_vertices(f->vertices[id].begin(), f->vertices[id].end());
			for (int i_v = 0; i_v < contour_vertices.size(); ++i_v) {
				Point2d pt = contour_vertices[i_v].first->pt + shift;
				SVG::markup_point(os, pt.x, rows - pt.y, 0.2, std::string("blue"), 0.1, std::string("blue"));
				std::string str_id = std::to_string(contour_vertices[i_v].first->id_vertex);
				SVG::markup_text(os, pt.x, rows - pt.y, str_id, 1, std::string("black"));
			}
		}
		for (auto se : f->split_edges_property) {
			Point2d pt1 = se.first.first->pt + shift;
			Point2d pt2 = se.first.second->pt + shift;
			Point2d pc1(jclamp(0, pt1.x, cols), jclamp(0, rows - pt1.y, rows));
			Point2d pc2(jclamp(0, pt2.x, cols), jclamp(0, rows - pt2.y, rows));
			LineItem seg_line(pc1.x, pc1.y, pc2.x, pc2.y, 0, 255, 0);
			SVG::markup_line(os, &seg_line, 0.2);

			std::string str_id = std::to_string(se.first.first->id_vertex);
			SVG::markup_text(os, pt1.x, rows - pt1.y, str_id, 1, std::string("black"));
			str_id = std::to_string(se.first.second->id_vertex);
			SVG::markup_text(os, pt2.x, rows - pt2.y, str_id, 1, std::string("black"));
		}
	}

	SVG::markup_footer(os);
	fb.close();
}

void Partition::subpartition_info_to_svg(const std::string & directory, const vector<Segment *> & segments, Face* f)
{
	Point2d shift(0.5, 0.5);
	int f_id = f != nullptr ? f->id_face : -1;
	std::cout << "save facet " << f_id << std::endl;
	std::string fileNameWithoutExtension = params->path_input_image.substr(0, params->path_input_image.rfind(".")).substr(params->path_input_image.rfind(separator)+1);
	if (f != nullptr) {
		std::string filename = directory + separator + fileNameWithoutExtension + "_parent_facet_" + std::to_string(f->id_face) + ".svg";
		std::filebuf fb;

		fb.open(filename, std::ios::out);
		std::ostream os(&fb);
		SVG::markup_header(os, int(rows), int(cols));
		// polygons and vertices
		SVG::markup_polygon(os, int(rows), int(cols), f, 0, 100, 255, 0.1, 0.4, shift.x, shift.y);
		for (int id = 0; id < f->vertices.size(); ++id) {
			vector<pair<Vertex *, Vertex_Type>> contour_vertices(f->vertices[id].begin(), f->vertices[id].end());
			for (int i_v = 0; i_v < contour_vertices.size(); ++i_v) {
				Point2d pt = contour_vertices[i_v].first->pt + shift;
				SVG::markup_point(os, pt.x, rows - pt.y, 0.2, std::string("blue"), 0.1, std::string("blue"));
				std::string str_id = std::to_string(contour_vertices[i_v].first->id_vertex);
				SVG::markup_text(os, pt.x, rows - pt.y, str_id, 1, std::string("black"));
			}
		}
		SVG::markup_footer(os);
		fb.close();
	}
	std::string filename = directory + separator + fileNameWithoutExtension + "_subpartition_facet_" + std::to_string(f_id) + "_into_" + std::to_string(faces.size()) + ".svg";
	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);

	SVG::markup_header(os, int(rows), int(cols));

	// Prepares palette
	map<int, Vec3b> colors;
	std::uniform_int_distribution<int> uniform_dist(100, 255);
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		int r = uniform_dist(generator);
		int g = uniform_dist(generator);
		int b = uniform_dist(generator);
		colors.insert(std::make_pair((*it_f)->id_face, Vec3b(r, g, b)));
	}
	// draw facets
	int k = 0;
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		SVG::markup_polygon(os, int(rows), int(cols), *it_f, colors[k][0], colors[k][1], colors[k][2], 0.1, 0.4, shift.x, shift.y);
		k++;
	}

	uchar color_0[3] = { 100, 240, 255 };
	uchar color_1[3] = { 0, 240, 240 };
	uchar color_2[3] = { 35, 120, 180 };
	uchar color_3[3] = { 105, 20, 105 };

	for (int i = 0; i < segments.size(); ++i) {
		Point2d pt1 = segments[i]->finalEnd1 + shift;
		Point2d pt2 = segments[i]->finalEnd2 + shift;
		Point2d pc1(jclamp(0, pt1.x, cols), jclamp(0, rows - pt1.y, rows));
		Point2d pc2(jclamp(0, pt2.x, cols), jclamp(0, rows - pt2.y, rows));
		LineItem seg_line;
		if (!segments[i]->is_disabled) {
			seg_line = LineItem(pc1.x, pc1.y, pc2.x, pc2.y, color_1[0], color_1[1], color_1[2]);
		}
		else {
			seg_line = LineItem(pc1.x, pc1.y, pc2.x, pc2.y, color_3[0], color_3[1], color_3[2]);
		}

		SVG::markup_line(os, &seg_line, 0.2);
	}

	// draw vertices
	Vertex* v = vertices_head;
	while (v != NULL) {
		Point2d pt = v->pt + shift;
		Point2d pc(jclamp(0, pt.x, cols), jclamp(0, rows - pt.y, rows));
		if (v->connectivity() == 1) {
			SVG::markup_point(os, pc.x, pc.y, 0.15, std::string("gray"), 0.05, std::string("gray"));
		}
		else if (v->connectivity() == 2) {
			SVG::markup_point(os, pc.x, pc.y, 0.15, std::string("black"), 0.05, std::string("green"));
		}
		else if (v->connectivity() == 3) {
			SVG::markup_point(os, pc.x, pc.y, 0.15, std::string("black"), 0.05, std::string("blue"));
		}
		else if (v->connectivity() == 4) {
			SVG::markup_point(os, pc.x, pc.y, 0.15, std::string("black"), 0.05, std::string("cyan"));
		}
		else {
			SVG::markup_point(os, pc.x, pc.y, 0.15, std::string("black"), 0.05, std::string("yellow"));
		}
		std::string str_id = std::to_string(v->id_vertex);
		SVG::markup_text(os, pc.x, pc.y, str_id, 1, std::string("black"));
		v = v->v_next;
	}


	SVG::markup_footer(os);

	fb.close();
}


void Partition::save_svg(const std::string & directory, const std::string &basename, bool show_details, std::string suffix)
{
	// The svg coordinates system: origin (0.0, 0.0) at the top-left corner of pixel (0,0), with y-axis pointing down.
	// We transform our coordinates system to the svg coordinates.
	Point2d shift(0.5, 0.5);

	std::string filename;
	if (show_details) {
		filename = directory + separator + basename + "_partition_" + std::to_string(faces.size()) + "_faces_" + suffix + ".svg";
	}
	else {
		filename = directory + separator + basename + "_partition_" + suffix + ".svg";
	}
	
	std::filebuf fb;

	fb.open(filename, std::ios::out);
	std::ostream os(&fb);

	SVG::markup_header(os, int(rows), int(cols));

	// Prepares palette for each label
	map<int, Vec3b> colors;
	std::uniform_int_distribution<int> uniform_dist(100, 255);
	for (int i = 0; i <= params->max_label; i++) {
		int id = i;
		int r = 0, g = 0, b = 0;
		for (int j = 0; j < 8; j++) {
			r = r | (((id >> 0) & 1) << 7 - j);
			g = g | (((id >> 1) & 1) << 7 - j);
			b = b | (((id >> 2) & 1) << 7 - j);
			id = id >> 3;
		}
		if (i > 0) colors.insert(std::make_pair(i, Vec3b(r, g, b)));
	}

	// color facets
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		int k = (*it_f)->label;
		if (k == 0) SVG::markup_path(os, int(rows), int(cols), (*it_f)->vertices[0], colors[k][0], colors[k][1], colors[k][2], 0.5, 0, shift.x, shift.y);
		else {
			int r = uniform_dist(generator);
			int g = uniform_dist(generator);
			int b = uniform_dist(generator);
			SVG::markup_path(os, int(rows), int(cols), *it_f, r, g, b, 1, 0.5, shift.x, shift.y);
		}
	}

	// draw vertices
	Vertex* v = vertices_head;
	while (v != NULL) {
		for (int j = 0; j < v->directions.size(); ++j) {
			Point2d pt = v->pt + shift;
			Point2d pc(jclamp(0, pt.x, cols), jclamp(0, rows - pt.y, rows));
			SVG::markup_point(os, pc.x, pc.y, 2, std::string("red"), 0.1, std::string("red"));
		}
		v = v->v_next;
	}
	
	// draw detected segments
	if (show_details) {
		uchar color_0[3] = { 100, 240, 255 };
		uchar color_1[3] = { 0, 240, 240 };
		uchar color_2[3] = { 35, 120, 180 };
		uchar color_3[3] = { 105, 20, 105 };
		for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); ++it_f) {
			vector<Segment *> f_split_support((*it_f)->detected_segments.begin(), (*it_f)->detected_segments.end());
			for (int i = 0; i < f_split_support.size(); ++i) {
				Point2d pt1 = f_split_support[i]->finalEnd1 + shift;
				Point2d pt2 = f_split_support[i]->finalEnd2 + shift;
				Point2d pc1(jclamp(0, pt1.x, cols), jclamp(0, rows - pt1.y, rows));
				Point2d pc2(jclamp(0, pt2.x, cols), jclamp(0, rows - pt2.y, rows));
				LineItem seg_line;
				if (i < 100) {
					if (f_split_support[i]->is_disabled) seg_line = LineItem(pc1.x, pc1.y, pc2.x, pc2.y, color_3[0], color_3[1], color_3[2]);
					else seg_line = LineItem(pc1.x, pc1.y, pc2.x, pc2.y, color_0[0], color_0[1], color_0[2]);
				}
				else {
					seg_line = LineItem(pc1.x, pc1.y, pc2.x, pc2.y, color_3[0], color_3[1], color_3[2]);
				}
				SVG::markup_line(os, &seg_line, 0.3);

				std::ostringstream str_score;
				str_score << std::fixed;
				str_score << std::setprecision(2) << f_split_support[i]->support_cost;
				Point2d pc_center = (pc1 + pc2) / 2;
				SVG::markup_text(os, pc_center.x, pc_center.y, str_score.str(), 2, std::string("blue"));
			}
		}
	}

	SVG::markup_footer(os);

	fb.close();
}


void Partition::save_ply(const std::string & directory, const std::string & basename)
{
	uint n_vertices = 0, n_edges = 0;
	Vertex* v = vertices_head;
	while (v != NULL) {
		++n_vertices;
		v = v->v_next;
	}
	Edge* e = edges_head;
	while (e != NULL) {
		++n_edges;
		e = e->e_next;
	}

	std::string filename;
	filename = directory + separator + basename + ".ply";

	std::filebuf fb;
	fb.open(filename, std::ios::out);
	std::ostream os(&fb);
	// header
	os << "ply" << std::endl;
	os << "format ascii 1.0" << std::endl;
	os << "element vertex " << n_vertices << std::endl;
	os << "property float x" << std::endl;
	os << "property float y" << std::endl;
	os << "property float z" << std::endl;
	os << "element edge " << n_edges << std::endl;
	os << "property int vertex1" << std::endl;
	os << "property int vertex2" << std::endl;
	os << "end_header" << std::endl;
	// vertices
	map<Vertex*, uint> v_vertices;
	uint v_id = 0;
	v = vertices_head;
	while (v != NULL) {
		os << v->pt.x << " " << v->pt.y << " " << 0.06 << std::endl;
		v_vertices.insert(pair<Vertex*, uint>(v, v_id++));
		v = v->v_next;
	}
	//edges
	e = edges_head;
	while (e != NULL) {
		os << v_vertices.at(e->v1) << " " << v_vertices.at(e->v2) << std::endl;
		e = e->e_next;
	}
	
	fb.close();
}

void Partition::save_graph_definition(const std::string &directory, const std::string &basename)
{
	int id_vertex = 0;
	int id_edge = 0;
	int id_face = 0;
    std::string filename = directory + separator + basename + ".txt";

	FILE* file = fopen(filename.c_str(), "w");
	if (file != NULL) {
		// If we have been able to successfully open a file, then we can start writing in it
		// First of all we set the size of the associated image
		// fprintf(file, "size\n");
        fprintf(file, "%i %i\n", rows, cols);
		fprintf(file, "%i %i %i\n", int(v_size), int(e_size), int(faces.size()));

		// We list the vertices
		fprintf(file, "vertices\n");
		Vertex* v = vertices_head;
		while (v != NULL) {
			v->id_vertex = id_vertex++;
			fprintf(file, "%i %lf %lf\n", v->id_vertex, v->pt.x, v->pt.y);
			v = v->v_next;
		}

		// We list the edges
		fprintf(file, "edges\n");
		Edge* e = edges_head;
		while (e != NULL) {
			e->id_edge = id_edge++;
			fprintf(file, "%i %i %i\n", e->id_edge, e->v1->id_vertex, e->v2->id_vertex);
			e = e->e_next;
		}

		// We list the faces
		fprintf(file, "faces\n");
		for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
			Face* f = (*it_f);
			f->id_face = id_face++;
			fprintf(file, "%i %i %i\n", f->id_face, int(f->edges.size()), f->label);
			for (int id = 0; id < f->edges.size(); ++id) {
				for (list<HalfEdge *>::iterator it_h = f->edges[id].begin(); it_h != f->edges[id].end(); it_h++) {
					HalfEdge* h = (*it_h);
					fprintf(file, "%i %i  ", h->e->id_edge, int(h->v1_v2));
				}
				fprintf(file, "\n");
			}
		}
		fclose(file);
	}
}


void Partition::save_boundaries(const std::string & directory, const std::string & basename)
{
    const uchar white[3] = {255, 255, 255};

    Matrix<uchar> B = Matrix<uchar>(rows, cols, 1);
    B.set(0);

    Edge* e = edges_head;
    while (e != NULL) {
		Point2d & pt1 = e->v1->pt;
		Point2d & pt2 = e->v2->pt;
		Point2i pi1 = Point2i(int(jclamp(0, pt1.x, int(cols) - 1)), int(jclamp(0, int(rows)-1 - pt1.y, int(rows) - 1)));
		Point2i pi2 = Point2i(int(jclamp(0, pt2.x, int(cols) - 1)), int(jclamp(0, int(rows)-1 - pt2.y, int(rows) - 1)));
		bool is_boundary = (pi1.x == 0 && pi2.x == 0) || (pi1.x == int(cols) - 1 && pi2.x == int(cols) - 1)
			|| (pi1.y == 0 && pi2.y == 0) || (pi1.y == int(rows) - 1 && pi2.y == int(rows) - 1);
		if (!is_boundary) {
			B.line(pi1.y, pi1.x, pi2.y, pi2.x, white);
		}
        e = e->e_next;
    }

    std::string filename = directory + separator + basename + "_boundaries" + (params->suffix == "" ? "" : "_" + params->suffix) + ".png";
    B.write_uchar(filename);
}
