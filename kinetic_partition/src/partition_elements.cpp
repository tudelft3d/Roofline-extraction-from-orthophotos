#include "partition_elements.h"
#include <iostream>
#include <map>
#include <opencv2/core/core.hpp>
#include "svg.h"
#include "geometry.h"

using cv::Mat;
using std::min;
using std::max;

Vertex::Vertex(int & _id_vertex, double x, double y)
{
	id_vertex = _id_vertex++;

	events = list<IndexedEvent *>();

	pt = Point2d(x, y);
	directions = vector<pair<double, HalfEdge *> >();

	v_prev = v_next = nullptr;
}


Vertex::Vertex(int & _id_vertex, IndexedEvent* event, double x, double y)
{
	id_vertex = _id_vertex++;

	events = list<IndexedEvent *>();
	if (event != nullptr) events.push_back(event);

	pt = Point2d(x, y);
	directions = vector<pair<double, HalfEdge *> >();

	v_prev = v_next = nullptr;
}


Vertex::Vertex(int & _id_vertex, IndexedEvent* event, Point2d & _pt)
{
	id_vertex = _id_vertex++;

	events = list<IndexedEvent *>();
	if (event != nullptr) events.push_back(event);
	
	pt = _pt;
	directions = vector<pair<double, HalfEdge*> >();

	v_prev = v_next = nullptr;
}


Vertex::~Vertex()
{
	for (list<IndexedEvent*>::iterator it_e = events.begin() ; it_e != events.end() ; it_e++) {
		delete (*it_e);
	}
	events.clear();

	assert(directions.size() == 0);
	//assert(faces.empty());
	faces.clear();
	directions.clear();
	v_prev = v_next = nullptr;
}


double Vertex::incidence_angle(const Edge* e)
{
	for (vector<pair<double, HalfEdge *> >::const_iterator it = directions.begin(); it != directions.end() ; ++it) {
		if (it->second->e == e) {
			return it->first;
		}
	}
	//assert(false);
	return 0;
}


void Vertex::add(double alpha, HalfEdge *h, bool check)
{
	Vertex* v = h->v1_v2 ? h->e->v2 : h->e->v1;
	vector<pair<double, HalfEdge *> >::iterator it;
	for (it = directions.begin(); it != directions.end(); it++) {
		//if (check) assert(it->second->e->v1 != v && it->second->e->v2 != v);
		if (it->first > alpha) break;
	}
	directions.insert(it, std::make_pair(alpha, h));
}


void Vertex::remove(HalfEdge *h)
{
	vector<pair<double, HalfEdge *> >::iterator it = directions.begin();
	while (it != directions.end()) {
		if (it->second == h) {
			break;
		}
		else {
			++it;
		}
	}
	directions.erase(it);
}


int Vertex::connectivity()
{
	return int(directions.size());
}


void Vertex::add(Face *f)
{
	faces.insert(f);
}


void Vertex::remove(Face *f)
{
	set<Face *>::iterator it;
	for (it = faces.begin(); it != faces.end(); it++) {
		if (*it == f) {
			faces.erase(it);
			break;
		}
	}
}


bool Vertex::outer_vertex_is_very_close(Point2d & _pt)
{
	return (pt.x == _pt.x && pt.y == _pt.y);
}


bool Vertex::inner_vertex_created_by_colinear_ray(Ray* r_intersectant, Ray* r_intersected, vector<Ray *> & rays, bool verbose)
{
	// We are going to check if this vertex has been created by another ray that belongs
	// to the group of colinear segments identified by 'node'
	Node_Colinear_Segments* node = r_intersectant->parent->node_colinear;

	for (list<IndexedEvent *>::iterator it_e = events.begin(); it_e != events.end(); it_e++) {
		IndexedEvent* e = *it_e;
		int intersectant = e->intersectant, intersected = e->intersected;
		if (intersected >= 0 && ((rays[intersected] == r_intersected && rays[intersectant]->parent->node_colinear == node) || (rays[intersectant] == r_intersected && rays[intersected]->parent->node_colinear == node))) {
			return true;
		}
	}

	return false;
}

bool Vertex::boundary_vertex_on_colinear_ray(Ray* r_intersectant, vector<Ray *> & rays)
{
	// We are going to check if this vertex lies on another ray that belongs
	// to the group of colinear segments identified by 'node'
	Node_Colinear_Segments* node = r_intersectant->parent->node_colinear;

	for (list<IndexedEvent *>::iterator it_e = events.begin(); it_e != events.end(); it_e++) {
		IndexedEvent* e = *it_e;
		int intersectant = e->intersectant, intersected = e->intersected;
		if ( rays[intersectant] != r_intersectant && rays[intersectant]->parent->node_colinear == node ) {
			return true;
		}
	}

	return false;
}


bool Vertex::same_intersection_point(Ray* r_intersectant, Ray* r_intersected)
{
	int id_intersectant = r_intersectant->index, id_intersected = r_intersected->index;

	for (list<IndexedEvent *>::iterator it_e = events.begin() ; it_e != events.end() ; it_e++) {
		IndexedEvent* e = (*it_e);
		int intersectant = e->intersectant, intersected = e->intersected;
		if ((intersectant == id_intersectant && intersected == id_intersected)
			|| (intersectant == id_intersected && intersected == id_intersectant))
			return true;
	}

	return false;
}


bool Vertex::has_no_duplicated_event()
{
	if (events.size() <= 1) return true;

	list<IndexedEvent*>::iterator it_e1, it_e2;
	for (it_e1 = events.begin() ; it_e1 != events.end() ; it_e1++) {
		int e1_intersectant = (*it_e1)->intersectant;
		int e1_intersected = (*it_e1)->intersected;
		it_e2 = it_e1;
		it_e2++;
		while (it_e2 != events.end()) {
			int e2_intersectant = (*it_e2)->intersectant;
			int e2_intersected = (*it_e2)->intersected;
			if (e1_intersectant == e2_intersectant && e1_intersected == e2_intersected) {
				if (e1_intersected == POLYGON_INNER) {
					if ((*it_e1)->intersected_index == (*it_e2)->intersected_index) return false;
				}
				else {
					return false;
				}
				
			}
			it_e2++;
		}
	}

	return true;
}


bool Vertex::right_turn(HalfEdge* h_prev, Vertex* v_curr, HalfEdge* h_next)
{
	// Handle colinearity cases
	if (h_prev->e->type == INNER_EDGE && h_next->e->type == INNER_EDGE) {
		Inner_Edge* e_prev = static_cast<Inner_Edge *>(h_prev->e);
		Inner_Edge* e_next = static_cast<Inner_Edge *>(h_next->e);
		Segment* s_prev = e_prev->get_front_supporting_segment();
		Segment* s_next = e_next->get_front_supporting_segment();
		if (s_prev->node_colinear != NULL && s_prev->node_colinear == s_next->node_colinear) {
			return true;
		}
	} else if (h_prev->e->type == OUTER_EDGE && h_next->e->type == OUTER_EDGE) {
		Outer_Edge* e_prev = static_cast<Outer_Edge *>(h_prev->e);
		Outer_Edge* e_next = static_cast<Outer_Edge *>(h_next->e);
		Image_Boundary b_prev = e_prev->boundary;
		Image_Boundary b_next = e_next->boundary;
		if (b_prev == b_next) {
			return true;
		}
	}

	double alpha_p = 0, alpha_n = 0;
	for (uint d = 0 ; d < uint(v_curr->directions.size()) ; ++d) {
		double & alpha_d = v_curr->directions[d].first;
		HalfEdge* h_d = v_curr->directions[d].second;
		if (h_d == h_prev) {
			alpha_p = alpha_d;
		} else if (h_d == h_next) {
			alpha_n = alpha_d;
		}
	}

	double eps = 1e-6;
	if (alpha_p <= 0) {
		// We should have alpha_n in [alpha_p, alpha_p + PI]
		return (fabs(alpha_p - alpha_n) < eps) || (fabs(alpha_p + PI - alpha_n) < eps) || (alpha_p < alpha_n && alpha_n < alpha_p + PI);
	} else {
		// We should have alpha_n in [-PI, alpha_p - PI] U [alpha_p, PI]
		bool first_interval = (fabs(alpha_n - (alpha_p - PI)) < eps) || (alpha_n < alpha_p - PI);
		if (first_interval) {
			return true;
		} else {
			bool second_interval = (fabs(alpha_n - alpha_p) < eps) || (alpha_p < alpha_n);
			return second_interval;
		}
	}
}


bool Vertex::right_turn(Vertex* v_prev, Vertex* v_curr, Vertex* v_next)
{
	HalfEdge *h_prev = nullptr, *h_next = nullptr;
	for (uint d = 0; d < v_curr->directions.size(); d++) {
		HalfEdge* h_d = v_curr->directions[d].second;
		if ((h_d->v1_v2 && h_d->e->v2 == v_prev) || (!h_d->v1_v2 && h_d->e->v1 == v_prev)) {
			h_prev = h_d;
		} else if ((h_d->v1_v2 && h_d->e->v2 == v_next) || (!h_d->v1_v2 && h_d->e->v1 == v_next)) {
			h_next = h_d;
		}
	}

	assert(h_prev != nullptr && h_next != nullptr);
	return right_turn(h_prev, v_curr, h_next);
}

bool Vertex::approximate_coordinates_polygon_corners(IndexedEvent* upcoming_event, Ray* r_i, Ray* r_j, double t_i, vector<list<Outer_Edge *> > & outer_edges, Point2d & pt, Vertex* & intersection_vertex, double min_edge_len, double corner_eps)
{
	// The return value is true iff pt is a corner of the polygon
	intersection_vertex = nullptr;

	double a = r_i->parent->a, b = r_i->parent->b, c = r_i->parent->c;
	double x = 0, y = 0;
	if (fabs(b) > fabs(a)) {
		x = r_i->A.x + t_i * r_i->OA[0];
		y = (-c - a * x) / b;
	}
	else {
		y = r_i->A.y + t_i * r_i->OA[1];
		x = (-c - b * y) / a;
	}
	pt = Point2d(x, y);

	if (r_j != nullptr) {
		return false;
	}
	else {
		// When a polygon boundary gets intersected. We deal with collinear and non-collinear cases separately.
		list<Outer_Edge *> & edges_for_this_corner = outer_edges[upcoming_event->intersected_index];
		Vertex* v1 = edges_for_this_corner.front()->v1, *v2 = edges_for_this_corner.back()->v2;
		
		if (r_i->collinear_boundaries.count(upcoming_event->intersected_index)>0) {
			// Collinear: output the polygon corner
			if (fabs(upcoming_event->t_intersected - 0) < min_edge_len) {
				intersection_vertex = v2;
				//assert(cv::norm(pt - v2->pt) < min_edge_len);
				pt = v2->pt; 
			}
			else {
				intersection_vertex = v1;
				//assert(cv::norm(pt - v1->pt) < min_edge_len);
				pt = v1->pt;
			}
			return true;
		}
		else {
			// Non-collinear: check if it intersects a corner of the intersected boundary: time on boundary is from -length to 0
			double dist_to_v2 = fabs(upcoming_event->t_intersected - 0), dist_to_v1 = fabs(upcoming_event->t_intersected + cv::norm(v2->pt - v1->pt));
			if (dist_to_v2 < dist_to_v1) {
				if (dist_to_v2 <= corner_eps) {
					intersection_vertex = v2;
					assert(cv::norm(pt - v2->pt) <= corner_eps);
					pt = v2->pt;
					return true;
				}
			}
			else {
				if (dist_to_v1 <= corner_eps) {
					intersection_vertex = v1;
					assert(cv::norm(pt - v1->pt) <= corner_eps);
					pt = v1->pt;
					return true;
				}
			}
		}
		return false;
	}
}

bool Vertex::approximate_coordinates(Ray* r_i, Ray* r_j, double t_i, int rows, int cols, Point2d & pt, double & corner_eps)
{
	// The return value is true iff pt is a corner of the image
	double a = r_i->parent->a, b = r_i->parent->b, c = r_i->parent->c;
	double x = 0, y = 0;
	if (fabs(b) > fabs(a)) {
		x = r_i->A.x + t_i * r_i->OA[0];
		y = (-c - a * x) / b;
	} else {
		y = r_i->A.y + t_i * r_i->OA[1];
		x = (-c - b * y) / a;
	}
	pt = Point2d(x, y);
	
	if (r_j != nullptr) {
		return false;
	} else {
		// Case when a boundary gets intersected
		return Vertex::round_coordinates(pt, rows, cols, corner_eps);
	}
}


bool Vertex::round_coordinates(Point2d & pt, int rows, int cols, double eps)
{
	// The return value is true iff pt is a corner of the image
	if (fabs(-0.5 - pt.x) < eps) {
		if (fabs(pt.y) < eps) {
			// Bottom left corner
			pt.x = -0.5;
			pt.y = -0.5;
			return true;
		} else if (fabs(double(rows)-0.5 - pt.y) < eps) {
			// Top left corner
			pt.x = -0.5;
			pt.y = double(rows) - 0.5;
			return true;
		} else {
			// Left border ?
			if (fabs(-0.5 - pt.x) < eps) pt.x = -0.5;
			return false;
		}
	} else if (fabs(double(cols)-0.5 - pt.x) < eps) {
		if (fabs(-0.5 - pt.y) < eps) {
			// Bottom right corner
			pt.x = double(cols) - 0.5;
			pt.y = -0.5;
			return true;
		} else if (fabs(double(rows)-0.5 - pt.y) < eps) {
			// Top right corner
			pt.x = double(cols) - 0.5;
			pt.y = double(rows) - 0.5;
			return true;
		} else {
			// Right border ?
			if (fabs(double(cols)-0.5 - pt.x) < eps) pt.x = double(cols) - 0.5;
			return false;
		}
	} else {
		// Bottom and top borders
		if (fabs(-0.5 - pt.y) < eps) {
			pt.y = -0.5;
		} else if (fabs(double(rows)-0.5 - pt.y) < eps) {
			pt.y = double(rows) - 0.5;
		}
		return false;
	}
}

// alpha_1: direcion angle for halfedge v1_v2, in (-pi, pi]
Edge::Edge(int & _id_edge, Vertex* _v1, Vertex* _v2, double _width, double & alpha_1, double & alpha_2, bool check)
{
	id_edge = _id_edge++;
	enabled = true;
	init_vertices_and_halfedges(_v1, alpha_1, _v2, alpha_2, check);
	Vec2d v1v2 = _v2->pt - _v1->pt;
	length = sqrt(v1v2[0] * v1v2[0] + v1v2[1] * v1v2[1]);

	width = _width;
	unit_edge_cost = -1;
	mean_angle_diff = 0;
	grad_center_offset = 0;
	grad_weight = 0;
	normal = Vec2d(-v1v2[1], v1v2[0])/length; // assume a random direction out of two

	e_prev = e_next = nullptr;
}



Artificial_Edge::Artificial_Edge(int & _id_edge, double _a, double _b, double _c, Vertex* _v1, Vertex* _v2, double _width, double _alpha_1, double _alpha_2)
	: Edge (_id_edge, _v1, _v2, _width, _alpha_1, _alpha_2)
{
	type = ARTIFICIAL_EDGE;
	a = _a;
	b = _b;
	c = _c;
}



Outer_Edge::Outer_Edge(int & _id_edge, Image_Boundary _boundary, Vertex* _v1, Vertex* _v2, double _width, double alpha_1, double alpha_2, bool _v1_2, double _t1, double _t2)
	: Edge (_id_edge, _v1, _v2, _width, alpha_1, alpha_2)
{
	type = OUTER_EDGE;
	boundary = _boundary;
	init_boundary_location = pair<int, int>(-1,-1);
	support_boundary_index = 0;
	v1_2 = _v1_2;
	dir_v1v2 = (_v2->pt - _v1->pt) / length;
	outer_face = nullptr;
}

Outer_Edge::Outer_Edge(int & _id_edge, Image_Boundary _boundary, Face* _outer_face, pair<int, int> _init_boundary_location, uint _support_boundary_index,
	Vertex* _v1, Vertex* _v2, double _width, double alpha_1, double alpha_2, bool _v1_2, Vec2d _dir_v1v2, double _t1, double _t2)
	: Edge(_id_edge, _v1, _v2, _width, alpha_1, alpha_2)
{
	type = OUTER_EDGE;
	boundary = _boundary;
	init_boundary_location = _init_boundary_location;
	support_boundary_index = _support_boundary_index;
	outer_face = _outer_face;
	v1_2 = _v1_2;
	dir_v1v2 = _dir_v1v2;
}

Outer_Edge::Outer_Edge(int & _id_edge, Image_Boundary _boundary, Face* _outer_face, pair<int, int> _init_boundary_location, uint _support_boundary_index,
	Vertex* _v1, Vertex* _v2, double _width, double alpha_1, double alpha_2, bool _v1_2, double _t1, double _t2)
	: Edge(_id_edge, _v1, _v2, _width, alpha_1, alpha_2)
{
	type = OUTER_EDGE;
	boundary = _boundary;
	init_boundary_location = _init_boundary_location;
	support_boundary_index = _support_boundary_index;
	outer_face = _outer_face;
	v1_2 = _v1_2;
	dir_v1v2 = (_v2->pt - _v1->pt) / length;
}

/*
Outer_Edge::Outer_Edge(const Outer_Edge & E)
	: Edge(E)
{
	type = OUTER_EDGE;
	boundary = E.boundary;
}
*/


Inner_Edge::Inner_Edge(int & _id_edge, Ray* _r, int _tag, Vertex* _v1, Vertex* _v2, double _width, double alpha_1, double alpha_2, double _t1, double _t2, bool check)
	: Edge (_id_edge, _v1, _v2, _width, alpha_1, alpha_2, check)
{
	type = INNER_EDGE;
	tag = _tag;
	rays.insert(_r);
}


Inner_Edge::Inner_Edge(int & _id_edge, const set<Ray*> & _support, int _tag, Vertex* _v1, Vertex* _v2, double _width, double alpha_1, double alpha_2, bool check)
	: Edge(_id_edge, _v1, _v2, _width, alpha_1, alpha_2, check)
{
	type = INNER_EDGE;
	tag = _tag;
	for (set<Ray *>::iterator it = _support.begin() ; it != _support.end() ; ++it) {
		rays.insert(*it);
	}
}


/*
Inner_Edge::Inner_Edge(const Inner_Edge & E)
	: Edge(E)
{
	type = INNER_EDGE;
	for (set<Ray *>::iterator it = E.rays.begin() ; it != E.rays.end() ; ++it) {
		rays.insert(*it);
	}
}
*/


Edge::~Edge()
{
	v1->remove(v1_v2);
	v2->remove(v2_v1);
	delete v1_v2;
	delete v2_v1;
	e_next = nullptr;
	e_prev = nullptr;
	region.clear();
}


Artificial_Edge::~Artificial_Edge()
{

}

Inner_Edge::~Inner_Edge()
{

}


Outer_Edge::~Outer_Edge()
{

}


void Edge::init_vertices_and_halfedges(Vertex* _v1, double & _alpha_1, Vertex* _v2, double & _alpha_2, bool check)
{
	v1 = _v1;
	v2 = _v2;

	v1_v2 = new HalfEdge(this, true);
	v2_v1 = new HalfEdge(this, false);
	
	alpha_1 = _alpha_1;
	alpha_2 = _alpha_2;

	v1->add(alpha_1, v1_v2, check);
	v2->add(alpha_2, v2_v1, check);
}

double Edge::get_alpha_1()
{
	return alpha_1;
}

double Edge::get_alpha_2()
{
	return alpha_2;
}

void Inner_Edge::add_supporting_ray(Ray* s)
{
	if (type == INNER_EDGE) {
		rays.insert(s);
	}
}


void Inner_Edge::add_supporting_rays(set<Ray *> & s)
{
	if (type == INNER_EDGE) {
		for (set<Ray *>::iterator it_s = s.begin(); it_s != s.end(); it_s++) {
			rays.insert(*it_s);
		}
	}
}


void Inner_Edge::get_supporting_rays(set<Ray *> & s)
{
	if (type == INNER_EDGE) {
		for (set<Ray *>::iterator it_s = rays.begin() ; it_s != rays.end() ; it_s++) {
			s.insert(*it_s);
		}
	}
}


void Inner_Edge::get_supporting_segments(set<Segment *> & s)
{
	if (type == INNER_EDGE) {
		for (set<Ray *>::iterator it_s = rays.begin() ; it_s != rays.end() ; it_s++) {
			s.insert((*it_s)->parent);
		}
	}
}


Ray* Inner_Edge::get_front_supporting_ray()
{
	if (type == INNER_EDGE) {
		Ray* r = *rays.begin();
		return r;
	} else {
		return nullptr;
	}
}


Segment* Inner_Edge::get_front_supporting_segment()
{
	if (type == INNER_EDGE) {
		Ray* r = *rays.begin();
		return r->parent;
	} else {
		return nullptr;
	}
}


void Inner_Edge::time_range(Ray* r, double & t_1, double & t_2)
{
	if (v1->events.empty()) {
		t_1 = -r->initial_length;
	} else {
		for (list<IndexedEvent *>::iterator it_e = v1->events.begin() ; it_e != v1->events.end() ; it_e++) {
			IndexedEvent* e = (*it_e);
			if (e->intersectant == r->index) {
				t_1 = e->t_intersectant;
				break;
			} else if (e->intersected == r->index) {
				t_1 = e->t_intersected;
				break;
			}
			else { // the event happened on the opposite ray at a location very close to O, or happened between some other rays
				t_1 = -r->initial_length;
			}
		}
	}

	if (v2->events.empty()) {
		t_2 = -r->initial_length;
	} else {
		for (list<IndexedEvent *>::iterator it_e = v2->events.begin() ; it_e != v2->events.end() ; it_e++) {
			IndexedEvent* e = (*it_e);
			if (e->intersectant == r->index) {
				t_2 = e->t_intersectant;
				break;
			} else if (e->intersected == r->index) {
				t_2 = e->t_intersected;
				break;
			}
			else { // the event happened on the opposite ray at a location very close to O, or happened between some other rays
				t_2 = -r->initial_length;
			}
		}
	}
}

double Edge::unit_grad_cost() {
	return compute_unit_grad_cost(unit_edge_cost, mean_angle_diff, grad_center_offset);
}

void Edge::disable() {
	enabled = false;
}

void Edge::weight(const Point2d & A, const Point2d & B, double & width, vector<Point2i> & region, Vec2d & normal, double & grad_center_offset, 
	Matrix<bool> & used, const Matrix<double> & I_m, const Matrix<double> & I_t, Point2d & offset, Split_Params & split_params, 
	const vector<double> & cdf_grad_m, double & unit_edge_cost, double & mean_angle_diff, double & grad_weight, bool check_used)
{
	// The coordinates system: origin (0.0, 0.0) is set at the center of bottom-left pixel, with y-axis pointing upwards.
	region.clear();

	// Step 1: compute edge discontinuity score
	double length = cv::norm(B - A);
	Point2d line_dir = (B - A) / length;
	int n_steps = max(int(ceil(length)),1);
	double stepsize = length/ n_steps;
	double half_stepsize = 0.5*stepsize;
	double spatial_thresh = 0.70711; // sqrt(2)/2 to include pixels covered by the edge 
	double bin_width = split_params.bin_width;
	double angle_thresh = split_params.angle;
	double min_dist_angle_score = split_params.min_dist_angle_score;
	// half window size
	double halfw = 0.5*width;
	int hwsize = ceil(halfw);

	int rows = static_cast<int>(I_m.rows);
	int cols = static_cast<int>(I_m.cols);

	list<Point2i> pixels_curr;
	list<Point2i> inactive_pixels, line_pixels_1, line_pixels_2;
	vector<list<Point2i>> line_pixels_list_1, line_pixels_list_2;
	int n_prop_1 = 0, n_prop_2 = 0;

	// assume a normal direction
	normal = Vec2d(-line_dir.y, line_dir.x);
	double normal_angle = atan2(normal[1], normal[0]);

	// move along the edge from vertex A to vertex B
	for (int l = 0; l <= n_steps; ++l) {
		list<Point2i> pixels_prev(std::move(pixels_curr));
		pixels_curr.clear();
		bool stop_1 = true, stop_2 = true;
		Point2d M(A + l*stepsize*line_dir);
		int x = M.x - offset.x, y = M.y - offset.y;
		if (M.x - offset.x - x > 0.5) x++;
		if (M.y - offset.y - y > 0.5) y++;
		double lower = (l == 0) ? -0.5 : -half_stepsize;
		double upper = (l == n_steps) ? 0.5 : half_stepsize;
		for (int dx = -hwsize; dx <= hwsize; ++dx) {
			for (int dy = -hwsize; dy <= hwsize; ++dy) {
				Point2i neighbor(x + dx, y + dy);
				Point2d disp(neighbor.x + offset.x - M.x, neighbor.y + offset.y - M.y);
				double axial_distance = line_dir.x*disp.x + line_dir.y*disp.y;
				if ((axial_distance > lower) && (axial_distance <= upper)) {
					double transverse_distance = fabs(normal[0]*disp.x + normal[1]*disp.y);
					if (transverse_distance <= halfw) {
						// check if it is outside the image boundary
						if ((neighbor.y < 0) || (neighbor.y > rows-1) || (neighbor.x < 0) || (neighbor.x > cols-1)) continue;
						pixels_curr.push_back(neighbor);
						if (std::find(pixels_prev.begin(), pixels_prev.end(), neighbor) != pixels_prev.end()) continue;
						bool active = false;
						int i = rows - 1 - neighbor.y; // rows
						int j = neighbor.x; // columns
						double neighbor_angle = I_t(i, j);
						if (neighbor_angle != NOTDEF && ( (!check_used) || (check_used && !used(i, j)) )) {
							double angle_diff = Geometry::angle_diff_abs(neighbor_angle, normal_angle);
							if (angle_diff < angle_thresh) {
								active = true;
								stop_1 = false;
								line_pixels_1.push_back(neighbor);
								used(i, j) = true;
							}
							else if (angle_diff > PI - angle_thresh) {
								active = true;
								stop_2 = false;
								line_pixels_2.push_back(neighbor);
								used(i, j) = true;
							}
						}
						if (!active) {
							inactive_pixels.push_back(neighbor);
						}
					}
				}
			}
		}
		if (!stop_1) {
			++n_prop_1;
		}
		else if (!line_pixels_1.empty()) {
			line_pixels_list_1.push_back(std::move(line_pixels_1));
			line_pixels_1.clear();
			n_prop_1 = 0;
		}
		if (!stop_2) {
			++n_prop_2;
		}
		else if (!line_pixels_2.empty()) {
			line_pixels_list_2.push_back(std::move(line_pixels_2));
			line_pixels_2.clear();
			n_prop_2 = 0;
		}
	}
	if (n_prop_1 > 0) { // add last part
		line_pixels_list_1.push_back(std::move(line_pixels_1));
		line_pixels_1.clear();
		n_prop_1 = 0;
	}
	if (n_prop_2 > 0 ) {
		line_pixels_list_2.push_back(std::move(line_pixels_2));
		line_pixels_2.clear();
		n_prop_2 = 0;
	}
	
	// encourage edge overlap with meaningful gradients
	unit_edge_cost = 0;
	mean_angle_diff = 0;
	grad_weight = 0;
	grad_center_offset = 0;
	int n_total = 0;
	for (int i = 0; i < line_pixels_list_1.size(); ++i) n_total += line_pixels_list_1[i].size();
	for (int i = 0; i < line_pixels_list_2.size(); ++i) n_total += line_pixels_list_2[i].size();

	region.reserve(max(n_total,1));
	double sum_cost = 0, sum_overlap_weight = 0;
	Vec2d grad_t_1(0, 0), grad_t_2(0, 0);
	for (int dir = 0; dir < 2; ++dir) {
		vector<list<Point2i>> & line_pixels_list = (dir == 0) ? line_pixels_list_1 : line_pixels_list_2;
		Point2d nv = (dir == 0) ? normal : -normal;
		double nv_angle = (dir == 0) ? normal_angle : normal_angle + PI;
		if (nv_angle > PI) nv_angle -= 2*PI; // send to range [-PI, PI]
		Vec2d & grad_t_temp = (dir == 0) ? grad_t_1 : grad_t_2;

		for (int k = 0; k < line_pixels_list.size(); ++k) {
			for (auto it = line_pixels_list[k].begin(); it != line_pixels_list[k].end(); ++it) {
				/* save region pixel */
				region.push_back(*it);
				/* check if it is covered by the edge */
				Point2d disp = Point2d(it->x, it->y) + offset - A;
				double transverse_dist = fabs(disp.ddot(nv));
				if (transverse_dist > spatial_thresh) continue;
				/* compute cost if covered by the edge*/
				double axial_dist = disp.ddot(line_dir);
				int i = rows - 1 - it->y; // rows
				int j = it->x; // columns
				double gx = cos(I_t(i, j)), gy = sin(I_t(i, j));
				double dbin = I_m(i, j) / bin_width;
				int bin = floor(dbin);
				if (bin >= cdf_grad_m.size()) bin = cdf_grad_m.size() - 1;
				// raw pixel score in [0, 1]
				double cdf_prev_bin = bin > 0 ? cdf_grad_m[bin - 1] : 0;
				double raw_pixel_score = cdf_prev_bin + (dbin - bin)*(cdf_grad_m[bin] - cdf_prev_bin);
				
				double d_angle = I_t(i, j) - nv_angle;
				d_angle = (d_angle > PI) ? d_angle - 2 * PI: d_angle;
				d_angle = (d_angle < -PI) ? d_angle + 2 * PI : d_angle;
				double pixel_score = compute_pixel_score(raw_pixel_score, d_angle);
				//double pixel_weight = compute_pixel_weight(transverse_dist);

				grad_t_temp += I_m(i, j)*Vec2d(gx, gy);
				double proj = I_m(i, j)*max(nv.x*gx + nv.y*gy, 0.);

				double d_min_pos = axial_dist - 0.5, d_max_pos = axial_dist + 0.5;
				double overlap = 1;
				if (d_max_pos > length) overlap -= d_max_pos - length;
				if (d_min_pos < 0) overlap -= 0 - d_min_pos;
				overlap = max(overlap, 0.);
				//double weighted_overlap = overlap*pixel_weight;
				double weighted_overlap = overlap;
				sum_cost += weighted_overlap*(1 - pixel_score);
				sum_overlap_weight += weighted_overlap;
				mean_angle_diff += d_angle*overlap;
				grad_center_offset += proj*overlap*transverse_dist;
				grad_weight += proj*overlap;
			}
		}
	}

	for (auto it = inactive_pixels.begin(); it != inactive_pixels.end(); ++it) {
		/* check if it is covered by the edge */
		Point2d disp = (Point2d)*it + offset - A;
		double transverse_dist = disp.ddot(normal);
		if (transverse_dist > spatial_thresh) continue;
		/* compute cost if covered by the edge*/
		double axial_dist = disp.ddot(line_dir);
		//double pixel_weight = compute_pixel_weight(transverse_dist);
		double d_min_pos = axial_dist - 0.5, d_max_pos = axial_dist + 0.5;
		double overlap = 1;
		if (d_max_pos > length) overlap -= d_max_pos - length;
		if (d_min_pos < 0) overlap -= 0 - d_min_pos;
		overlap = max(overlap, 0.);
		//double weighted_overlap = overlap*pixel_weight;
		double weighted_overlap = overlap;
		sum_cost += weighted_overlap;
		sum_overlap_weight += weighted_overlap;
	}
	double norm_grad_t_1 = cv::norm(grad_t_1), norm_grad_t_2 = cv::norm(grad_t_2);
	normal = (norm_grad_t_1 >= norm_grad_t_2) ? normal : -normal;
	grad_center_offset /= max(grad_weight, 1e-5);
	mean_angle_diff /= max(grad_weight, 1e-5);
	unit_edge_cost = compute_unit_grad_cost(sum_cost / max(sum_overlap_weight, 1e-5), mean_angle_diff, grad_center_offset);
}

// calculate the weight on each edge, i.e. the discontinuity along the edge
void Edge::weight(Matrix<bool> & used, const Matrix<double> & I_m, const Matrix<double> & I_t, Point2d & offset, Split_Params & split_params, const vector<double> & cdf_grad_m)
{
	/*weight(v1->pt, v2->pt, width, region, normal, grad_center_offset, used, I_m, I_t,
		offset, split_params, cdf_grad_m, unit_edge_cost, mean_angle_diff, grad_weight, type != OUTER_EDGE);*/
	weight(v1->pt, v2->pt, width, region, normal, grad_center_offset, used, I_m, I_t,
		offset, split_params, cdf_grad_m, unit_edge_cost, mean_angle_diff, grad_weight, false);
}


void Inner_Edge::line_equation(double &a, double &b, double &c)
{
	if (type == INNER_EDGE) {
		Segment* s_ref = (*rays.begin())->parent;
		a = s_ref->a, b = s_ref->b, c = s_ref->c;
	}
}


void Outer_Edge::line_equation(double &a, double &b, double &c)
{
	if (type == OUTER_EDGE) {
		if (boundary == TOP_IMAGE || boundary == BOTTOM_IMAGE) {
			a = 0, b = 1, c = -v1->pt.y;
		} else if (boundary == LEFT_IMAGE || boundary == RIGHT_IMAGE) {
			a = 1, b = 0, c = -v1->pt.x;
		}
	}
}


void Outer_Edge::time_range(double & t_1, double & t_2, const Point2d & O, int intersected_index)
{ // Point O should be the reference point at time = 0.

	bool v1_intersected = false, v2_intersected = false;
	for (list<IndexedEvent *>::iterator it_e = v1->events.begin(); it_e != v1->events.end(); it_e++) {
		IndexedEvent* e = (*it_e);
		if (e->intersected < 0 && e->intersected == intersected_index) {
			t_1 = e->t_intersected;
			v1_intersected = true;
			break;
		}
	}

	if (!v1_intersected) {
		switch (boundary) {
		case TOP_IMAGE: t_1 = v1->pt.x; break;
		case BOTTOM_IMAGE: t_1 = v1->pt.x; break;
		case LEFT_IMAGE: t_1 = v1->pt.y; break;
		case RIGHT_IMAGE: t_1 = v1->pt.y; break;
		case POLYGON_OUTER: t_1 = dir_v1v2.ddot(v1->pt - O); break;
		case POLYGON_INNER: t_1 = dir_v1v2.ddot(v1->pt - O); break;
		}
	}

	for (list<IndexedEvent *>::iterator it_e = v2->events.begin(); it_e != v2->events.end(); it_e++) {
		IndexedEvent* e = (*it_e);
		if (e->intersected < 0 && e->intersected == intersected_index) {
			v2_intersected = true;
			t_2 = e->t_intersected;
			break;
		}
	}

	if (!v2_intersected) {
		switch (boundary) {
		case TOP_IMAGE: t_2 = v2->pt.x; break;
		case BOTTOM_IMAGE: t_2 = v2->pt.x; break;
		case LEFT_IMAGE: t_2 = v2->pt.y; break;
		case RIGHT_IMAGE: t_2 = v2->pt.y; break;
		case POLYGON_OUTER: t_2 = dir_v1v2.ddot(v2->pt - O); break;
		case POLYGON_INNER: t_2 = dir_v1v2.ddot(v2->pt - O); break;
		}
	}
}

void Artificial_Edge::line_equation(double & _a, double & _b, double & _c)
{
	_a = a; b = _b; c = _c;
}


HalfEdge::HalfEdge(Edge* _e, bool _v1_v2)
{
	e = _e;
	v1_v2 = _v1_v2;
	f = nullptr;
}


HalfEdge::~HalfEdge()
{

}


void HalfEdge::set(Face *_f)
{
	f = _f;
}


HalfEdge* HalfEdge::opposite()
{
	return (e->v1_v2 == this ? e->v2_v1 : e->v1_v2);
}



void HalfEdge::intersects_if_extended(HalfEdge* h, list<HalfEdge *> & intersectable_halfedges, Point2d & intersection_point, HalfEdge* & intersected)
{
	std::set<Segment *> support_extended;
	Inner_Edge* e = static_cast<Inner_Edge*>(h->e);
	e->get_supporting_segments(support_extended);
	Segment* s_extended = (*support_extended.begin());

	// Finds the point from which we extend the ray, finds its direction as well
	Point2d h_v = s_extended->finalBarycenter;
	Vec2d h_dir = Vec2d();

	Vec2d u = (h->v1_v2 ? h->e->v2->pt - h->e->v1->pt : h->e->v1->pt - h->e->v2->pt);
	Vec2d s_extended_dir = Vec2d(s_extended->finalEnd2 - s_extended->finalEnd1);
	if (u.ddot(s_extended->finalDirection) > 0) {
		h_dir = s_extended->finalDirection;
	} else {
		h_dir = -s_extended->finalDirection;
	}

	double eps = 1e-4;
	double t_min = FLT_MAX;
	for (list<HalfEdge *>::iterator it_h = intersectable_halfedges.begin() ; it_h != intersectable_halfedges.end(); ++it_h) {
	
		// Gets an edge
		HalfEdge* h_arg = (*it_h);
		Edge* arg = h_arg->e;

		// Gets its equation
		double a, b, c;
		arg->line_equation(a, b, c);

		// Extends rays and sees when its intersects the supporting line of the edge
		// Given that this function will only be used in the merging process...
		double den = (a * h_dir[0] + b * h_dir[1]);
		if (fabs(den) < 1e-6) {
			continue;
		} else {
			double t = -(a * h_v.x + b * h_v.y + c) / den;
			if (t < t_min) {
				intersection_point = Point2d(h_v.x + t * h_dir[0], h_v.y + t * h_dir[1]);
				
				Vertex *arg_v1 = arg->v1, *arg_v2 = arg->v2;
				double dx = arg_v2->pt.x - arg_v1->pt.x;
				double dy = arg_v2->pt.y - arg_v1->pt.y;
				bool check_x = true, check_y = true;
				
				if (fabs(dx) > 1e-6) {
					double u_x = (intersection_point.x - arg_v1->pt.x) / dx;
					check_x = (-eps <= u_x && u_x <= 1 + eps);
				}

				if (fabs(dy) > 1e-6) {
					double u_y = (intersection_point.y - arg_v1->pt.y) / dy;
					check_y = (-eps <= u_y && u_y <= 1 + eps);
				}

				if (check_x && check_y) {
					t_min = t;
					intersected = h_arg;
					return;
				}
			}
		}
	}
}


int HalfEdge::is_clockwise(list<HalfEdge *> & _edges)
{
	int clockwise = 1;
	// get the order of the contour (shoelace formula)
	
	double twice_area = 0;
	for (list<HalfEdge *>::iterator it_h = _edges.begin(); it_h != _edges.end(); ++it_h) {
		Vertex* v1 = (*it_h)->v1_v2 ? (*it_h)->e->v1 : (*it_h)->e->v2;
		Vertex* v2 = (*it_h)->v1_v2 ? (*it_h)->e->v2 : (*it_h)->e->v1;
		twice_area += (v2->pt.x - v1->pt.x)*(v2->pt.y + v1->pt.y);
	}
	// connect v_end-->v_start if the edges do not form a closed loop
	Vertex* v_start = _edges.front()->v1_v2 ? _edges.front()->e->v1: _edges.front()->e->v2;
	Vertex* v_end = _edges.back()->v1_v2 ? _edges.back()->e->v2 : _edges.back()->e->v1;
	if (v_start != v_end) {
		twice_area += (v_start->pt.x - v_end->pt.x)*(v_start->pt.y + v_end->pt.y);
	}

	if (twice_area > 0) { // clockwise
		clockwise = 1;
	}
	else if (twice_area < 0) { // counter-clockwise
		clockwise = 0;
	}
	else { // line
		clockwise = -1;
	}
	return clockwise;
}


Face::Face(int & _id_face, vector<list<HalfEdge *>> & _edges)
{
	id_face = _id_face++;
	label = INVALID_LABEL;
	check_split = true;

	// Copies the list of edges
	edges.resize(_edges.size());
	for (int id = 0; id < _edges.size(); ++id) {
		edges[id] = list<HalfEdge *>(_edges[id].begin(), _edges[id].end());
	}

	// Reciprocally, for every edge we insert a reference to this face
	for (int id = 0; id < edges.size(); ++id) {
		for (list<HalfEdge *>::iterator it = edges[id].begin(); it != edges[id].end(); it++) {
			(*it)->set(this);
		}
	}

	// We are now going to identify the list of vertices with the help of the half-edges
	list_vertices();

	classify();

	// Enable this to see non-simple faces
    //is_simple(true);

	pixels = vector<Point2i>();

	active_gradients = vector<set<Point2i, pointcomp>>();

	semantic_probabilities = vector<double>();

	detected_segments = set<Segment *, costcomp>();
	split_gain = -10000;

	split_vertex_outer_edge_location = map<Vertex*, pair<int, int>>();
	split_vertex_chains = vector<vector<Vertex *>>();
	split_edges_property = map<pair<Vertex *, Vertex *>, Split_Edge_Property>();
}


Face::Face(int & _id_face, list<HalfEdge *> & _edges)
{
	id_face = _id_face++;
	label = INVALID_LABEL;
	check_split = true;

	// Copies the list of edges
	edges.resize(1);
	edges[0] = list<HalfEdge *>(_edges.begin(), _edges.end());

	// Reciprocally, for every edge we insert a reference to this face
	for (int id = 0; id < edges.size(); ++id) {
		for (list<HalfEdge *>::iterator it = edges[id].begin(); it != edges[id].end(); it++) {
			(*it)->set(this);
		}
	}

	// We are now going to identify the list of vertices with the help of the half-edges
	list_vertices();

	classify();

	// Enable this to see non-simple faces
	//is_simple(true);

	pixels = vector<Point2i>();

	active_gradients = vector<set<Point2i, pointcomp>>();

	semantic_probabilities = vector<double>();

	detected_segments = set<Segment *, costcomp>();
	split_gain = -10000;

	split_vertex_outer_edge_location = map<Vertex*, pair<int, int>>();
	split_vertex_chains = vector<vector<Vertex *>>();
	split_edges_property = map<pair<Vertex *, Vertex *>, Split_Edge_Property>();
}


Face::~Face()
{
	// Removes references to this facet in every container
	for (int id = 0; id < edges.size(); ++id) {
		for (list<HalfEdge *>::iterator it_h = edges[id].begin(); it_h != edges[id].end(); it_h++) {
			(*it_h)->set(nullptr);
		}
	}
	edges.clear();

	for (int id = 0; id < vertices.size(); ++id) {
		for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = vertices[id].begin(); it_v != vertices[id].end(); it_v++) {
			(*it_v).first->remove(this);
		}
	}
	vertices.clear();
	
	pixels.clear();

	active_gradients.clear();
	active_lab_gradients.clear();

	for (set<Segment *, costcomp>::iterator it_s = detected_segments.begin(); it_s != detected_segments.end(); ++it_s) {
		Segment* s = *it_s;
		if (s->to_delete) {
			delete s;
			s = nullptr;
		}
	}
	detected_segments.clear();
	split_vertex_outer_edge_location.clear();
	split_vertex_chains.clear();
	split_edges_property.clear();

	for (int i = 0; i < split_vertex_chains.size(); ++i) {
		for (int j = 0; j < split_vertex_chains[i].size(); ++j) {
			Vertex *v = split_vertex_chains[i][j];
			if (v->v_prev == nullptr && v->v_next == nullptr) { // do not delete facet vertices (those still exist in the current partition)
				delete v;
				v = nullptr;
			}
		}
	}
}


bool Face::is_valid()
{
	bool valid_loop = true;
	for (int id = 0; id < edges.size(); ++id) {
		HalfEdge* first_h = edges[id].front();
		Vertex* start_v = first_h->v1_v2 ? first_h->e->v1 : first_h->e->v2;
		HalfEdge* last_h = edges[id].back();
		Vertex* end_v = last_h->v1_v2 ? last_h->e->v2 : last_h->e->v1;
		if (end_v != start_v) {
			valid_loop = false;
			break;
		}
	}
	return valid_loop;
}


// if edges cross each other, the face is considered not simple
bool Face::is_simple(bool verbose)
{
	assert(vertices.size() == edges.size());

	for (int id = 0; id < vertices.size(); ++id) {
		assert(vertices[id].size() == edges[id].size());
		std::map<int, int> vertices_occurences;
		for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = vertices[id].begin(); it_v != vertices[id].end(); it_v++) {
			vertices_occurences.insert(std::make_pair((*it_v).first->id_vertex, 0));
		}

		for (list<HalfEdge *>::iterator it_h = edges[id].begin(); it_h != edges[id].end(); it_h++) {
			int v1 = (*it_h)->e->v1->id_vertex;
			int v2 = (*it_h)->e->v2->id_vertex;
			vertices_occurences[v1]++;
			vertices_occurences[v2]++;
		}

		for (std::map<int, int>::iterator it_m = vertices_occurences.begin(); it_m != vertices_occurences.end(); it_m++) {
			if (it_m->second != 2) {
				vector<double> x;
				vector<double> y;
				vector<int> l;

				for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = vertices[id].begin(); it_v != vertices[id].end(); it_v++) {
					Vertex* v = it_v->first;
					x.push_back(v->pt.x);
					y.push_back(v->pt.y);
					l.push_back(v->id_vertex);
				}
				FILE* file = fopen("error.R", "w");
				fprintf(file, "x <- c("); for (int i = 0; i < x.size() - 1; i++) fprintf(file, "%lf, ", x[i]); fprintf(file, "%lf)\n", x[x.size() - 1]);
				fprintf(file, "y <- c("); for (int i = 0; i < y.size() - 1; i++) fprintf(file, "%lf, ", y[i]); fprintf(file, "%lf)\n", y[y.size() - 1]);
				fprintf(file, "l <- c("); for (int i = 0; i < l.size() - 1; i++) fprintf(file, "%i, ", l[i]); fprintf(file, "%i)\n", l[l.size() - 1]);
				fclose(file);
				if (verbose) {
					std::cout << "****** ERROR : Face " << id_face << " is not simple" << std::endl;
				}
				return false;
			}
		}
	}

    return true;
}


void Face::list_vertices()
{
	vertices.resize(edges.size());
	for (int id = 0; id < edges.size(); ++id) {
		vertices[id].clear();
		list<HalfEdge *>::iterator it_h = edges[id].begin();
		HalfEdge* h_curr = (*it_h);
		Vertex *v = nullptr, *v_2 = nullptr;
		if (h_curr->v1_v2) {
			v = h_curr->e->v1;
			//v_2 = h_curr->e->v2;
		}
		else {
			v = h_curr->e->v2;
			//v_2 = h_curr->e->v1;
		}
		v->add(this);
		vertices[id].push_back(std::make_pair(v, Vertex_Type::UNDETERMINED));
		while (true) {
			it_h++;
			if (it_h != edges[id].end()) {
				h_curr = (*it_h);
				if (h_curr->v1_v2) {
					v = h_curr->e->v1;
				}
				else {
					v = h_curr->e->v2;
				}
				v->add(this);
				vertices[id].push_back(std::make_pair(v, Vertex_Type::UNDETERMINED));
			}
			else {
				break;
			}
		}
	}
	
}

void Face::compute_merged_histogram(const list<vector<double>> & hists, vector<double> & merged_hist)
{
	merged_hist= vector<double>(hists.back().size(), 0);
	for (auto& histogram: hists) {
		for (int i = 0; i < histogram.size(); ++i) merged_hist[i] += histogram[i];
	}
}

void Face::compute_feature(Matrix<double> & I_prob, Matrix<uchar> & I, const Point2d & offset, int & max_label)
{
	if (this->pixels.empty()) this->find_pixels_inside_facet(offset);

	// create a semantic_probabilities vector to record class probability histogram
	histogram(I_prob, I, this->pixels, this->semantic_probabilities, max_label);
}

void Face::histogram(Matrix<double> & I_prob, Matrix<uchar> & I, vector<Point2i> & pixels, 
	vector<double> & hist, int & max_label)
{
	int n_pixels = pixels.size();
	int n_bins = max_label + 1;
	hist.clear();
	hist.resize(n_bins, 0.);
	
	// Each pixel contributes to the histogram in the following way:
	// Percentage of intersection of the pixel intensity block is added to bins of the 1d histogram
	for (auto it_p = pixels.begin(); it_p != pixels.end(); ++it_p) {
		uint i = jclamp(0, int(I.rows) - 1 - it_p->y, I.rows - 1); // rows
		uint j = uint(it_p->x); // columns
		// Class info
		double sum_foreground = 0;
		for (uint k = 0; k < I_prob.channels; k++) {
			// label id = 1 + channel index of I_prob
			hist[k + 1] += I_prob(i,j,k);
			sum_foreground += I_prob(i, j, k);
		}
		// compute bin height for the background label 0
		hist[0] += max(1 - sum_foreground, 0.);
	}
}


void Face::estimate_label(double alpha, int area, double diag_length) {
	bool ignore_boundary_conditions = (alpha < 1e-6);
	estimate_label(alpha, area, diag_length, ignore_boundary_conditions);
}

void Face::estimate_label(double & alpha, int & area, double & diag_length, bool ignore_boundary_conditions) {
	
	/*compute energy for each label*/
	vector<float> energy(semantic_probabilities.size(), 0);
	vector<float> prior_edge(semantic_probabilities.size(), 0);
	int n_pixels = pixels.size();
	for (int i = 0; i < semantic_probabilities.size(); i++) {
		energy[i] += (1 - alpha)*(n_pixels - semantic_probabilities[i]) / float(area);
		//energy[i] += (1 - alpha)*hist[i] / double(area);
	}
	if (!ignore_boundary_conditions) {
		for (int j = 0; j < edges.size(); j++) {
			for (auto it = edges[j].begin(); it != edges[j].end(); ++it) {
				Face* f_adj = (*it)->opposite()->f;
				float boundary_energy_edge_consistent = fidelity_edge(true, (*it)->e->unit_edge_cost, (*it)->e->length);
				float boundary_energy_edge_inconsistent = fidelity_edge(false, (*it)->e->unit_edge_cost, (*it)->e->length);
				for (auto it_h = prior_edge.begin(); it_h != prior_edge.end(); ++it_h) (*it_h) += boundary_energy_edge_inconsistent;
				if (f_adj != nullptr && f_adj->label != INVALID_LABEL) {

					auto it_hist = prior_edge.begin();
					std::advance(it_hist, f_adj->label);
					(*it_hist) += boundary_energy_edge_consistent - boundary_energy_edge_inconsistent;
				}
			}
		}
		for (int i = 0; i < energy.size(); i++) energy[i] += alpha*(prior_edge[i] / float(diag_length));
	}

	/*find minimum energy*/
	double min_energy = energy[0];
	label = 0;
	for (int i = 0; i < semantic_probabilities.size(); i++) {
		if (energy[i] < min_energy) {
			min_energy = energy[i];
			label = i;
		}
	}
}


void Face::estimate_label(vector<double> & hist, int & n_pixels, int & label, 
	double & alpha, int& area, double& diag_length, vector<vector<Adj_Edge_Info>*> & adj_edges_info, bool ignore_boundary_conditions) {
	
	/*compute fidelity for each label*/
	vector<float> energy(hist.size(), 0);
	for (int i = 0; i < hist.size(); i++) {
		energy[i] += (1 - alpha)*(n_pixels - hist[i]) / float(area);
	}

	if (!ignore_boundary_conditions) {
		vector<float> prior_edge(hist.size(), 0);
		for (int j = 0; j < adj_edges_info.size(); ++j) {
			auto it = adj_edges_info[j]->begin();
			auto end_it = adj_edges_info[j]->end();
			while (it != end_it) {
				float boundary_energy_edge_consistent = fidelity_edge(true, it->e->unit_edge_cost, it->e->length);
				float boundary_energy_edge_inconsistent = fidelity_edge(false, it->e->unit_edge_cost, it->e->length);
				for (auto it_h = prior_edge.begin(); it_h != prior_edge.end(); ++it_h) (*it_h) += boundary_energy_edge_inconsistent;
				if (it->label != INVALID_LABEL) {
					auto it_hist = prior_edge.begin();
					std::advance(it_hist, it->label);
					(*it_hist) += boundary_energy_edge_consistent - boundary_energy_edge_inconsistent;
				}
				++it;
			}
		}
		for (int i = 0; i < energy.size(); i++) energy[i] += alpha*(prior_edge[i] / float(diag_length));
	}

	/*find minimum energy*/
	double min_energy = energy[0];
	label = 0;
	for (int i = 0; i < energy.size(); i++) {
		if (energy[i] < min_energy) {
			min_energy = energy[i];
			label = i;
		}
	}
}

void Face::find_active_gradients(Matrix<double> & I_m, Matrix<double> & I_t, Point2d & offset, Split_Params & split_params, bool is_color_image)
{
	vector<set<Point2i, pointcomp>> & bins = is_color_image ? active_gradients : active_lab_gradients;

	// set parameters
	double value_thresh = is_color_image ? split_params.m_value_thresh : split_params.lab_m_value_thresh;
	int n_bins = ceil((0.5 - value_thresh)/ split_params.bin_width);
	int m_rows = static_cast<int>(I_m.rows);
	bins = vector<set<Point2i, pointcomp>>(n_bins, set<Point2i, pointcomp>());

	// apply value_thresh on pixels and save to a histogram with range [value_thresh, 1]
	for (auto it_p = pixels.begin(); it_p != pixels.end(); ++it_p) {
		int i = jclamp(0, m_rows -1 - it_p->y,  m_rows -1);
		int j = it_p->x; // columns
		double value_p = I_m(i, j);
		if (value_p > value_thresh && I_t(i, j) != NOTDEF) {
			int id_bin = floor((value_p - value_thresh) / split_params.bin_width);
			if (id_bin >= n_bins) id_bin = n_bins - 1;
			//std::cout<< "id_bin = "<< id_bin << "   " << bins[id_bin].size();
			bins[id_bin].insert(*it_p);
			//std::cout<< "..." << std::endl;
		}
	}

	// remove pixels that belong to edges
	if (is_color_image) {
		discard_grad_pixels_on_edges(bins, I_m, value_thresh, split_params.bin_width);
	}
}

void Face::discard_grad_pixels_on_edges(vector<set<Point2i, pointcomp>> & bins, Matrix<double> & I_m, double value_thresh, double bin_width) {
	int max_id_bin = bins.size()-1;
	int m_rows = static_cast<int>(I_m.rows);
	for (int i = 0; i < edges.size(); ++i) {
		for (auto edge_it = edges[i].begin(); edge_it != edges[i].end(); ++edge_it) {
			assert((*edge_it)->e->unit_edge_cost >= 0);
			const vector<Point2i> & edge_pixels = (*edge_it)->e->region;
			for (auto it_p = edge_pixels.begin(); it_p != edge_pixels.end(); ++it_p) {
				int i = m_rows - 1 - it_p->y; // rows
				int j = it_p->x; // columns
				int id_bin = jclamp(0, floor((I_m(i, j, 0) - value_thresh) / bin_width), max_id_bin);
				bins[id_bin].erase(*it_p);
			}
		}
	}
}

void Face::delete_invalid_segments() {
	for (auto it = detected_segments.begin(); it != detected_segments.end(); ++it) {
		if ((*it)->to_delete) {
			Segment* s_to_delete = *it;
			it = detected_segments.erase(it);
			delete s_to_delete;
			--it;
		}
	}
}

void Face::find_pixels_inside_facet(const Point2d & offset)
{
	// The coordinates system: origin (0.0, 0.0) is set at the center of pixel (rows-1, 0), with y-axis pointing upwards.
	// Therefore the pixel (rows-1, 0) spans the area of [-0.5, 0.5) X [-0.5, 0.5). This function records the set of coordinates of pixel centers. 

	pixels.clear();

	// Initialize a grid for the outer contour, covering the integers included in K = [x_min x_max[ x [y_min y_max[
	double x_min, x_max, y_min, y_max;
	list<HalfEdge *>::iterator it_l, it_r;
	find_leftmost_and_rightmost_edges(0, x_min, x_max, y_min, y_max, it_l, it_r);
	// Deal with precision at the boundary
	x_min = x_min < 0 ? -0.5 : x_min;
	y_min = y_min < 0 ? -0.5 : y_min;
	
	int x_0 = int(ceil(x_min));
	int y_0 = (y_min == ceil(y_min) ? int(y_min) + 1 : int(ceil(y_min))); // bottom
	int x_1 = (x_max == floor(x_max) ? int(floor(x_max)) - 1 : int(floor(x_max)));
	int y_1 = int(floor(y_max)); // top

	if (y_0 > y_1 || x_0 > x_1) return;

	int g_x = x_1 - x_0 + 1;
	int g_y = y_1 - y_0 + 1;
	Mat grid(g_y, g_x, CV_8U, cv::Scalar::all(1));
	update_pixel_grid(0, grid, x_0, y_1, x_min, x_max, y_min, y_max, offset, it_l, it_r);

	for (int id = 1; id < this->edges.size(); ++id) {
		find_leftmost_and_rightmost_edges(id, x_min, x_max, y_min, y_max, it_l, it_r);
		update_pixel_grid(id, grid, x_0, y_1, x_min, x_max, y_min, y_max, offset, it_l, it_r);
	}

	// Now use the grid to find pixels inside the facet
	int n_pixels = 0;
	for (int y = 0; y < g_y; ++y) {
		for (int x = 0; x < g_x; ++x) {
			if (grid.at<uchar>(y, x) == 1) ++n_pixels;
		}
	}
	pixels.reserve(n_pixels);
	for (int y = y_0; y <= y_1; ++y) {
		for (int x = x_0; x <= x_1; ++x) {
			if (grid.at<uchar>(y_1 - y, x - x_0) == 1) pixels.push_back(Point2i(x, y));
		}
	}
	grid.release();
}


void Face::find_leftmost_and_rightmost_edges(int contour_id, double & x_min, double & x_max, double & y_min, double & y_max, list<HalfEdge *>::iterator & it_l, list<HalfEdge *>::iterator & it_r)
{
	// We are searching for the leftmost and rightmost half-edges of the face
	list<HalfEdge *>::iterator it_h = edges[contour_id].begin();
	
	// Invariant : v1 and v2 represent the vertices of the previously crossed half-edge
	Vertex* v1 = nullptr;
	Vertex* v2 = nullptr;
	if ((*it_h)->v1_v2) {
		v1 = (*it_h)->e->v1;
		v2 = (*it_h)->e->v2;
	} else {
		v2 = (*it_h)->e->v1;
		v1 = (*it_h)->e->v2;
	}
	
	x_min = (v1->pt.x < v2->pt.x ? v1->pt.x : v2->pt.x);
	x_max = (v1->pt.x > v2->pt.x ? v1->pt.x : v2->pt.x);
	y_min = (v1->pt.y < v2->pt.y ? v1->pt.y : v2->pt.y);
	y_max = (v1->pt.y > v2->pt.y ? v1->pt.y : v2->pt.y);
	it_l = it_h, it_r = it_h;

	bool shift_l, shift_r;
	if (x_min == v2->pt.x) {
		shift_l = true;
		shift_r = false;
	} else if (x_max == v2->pt.x) {
		shift_l = false;
		shift_r = true;
	}
	it_h++;
	Vertex* v = nullptr;
	while (it_h != edges[contour_id].end()) {
		HalfEdge* h = *it_h;
		if (h->e->v1 == v1 || h->e->v1 == v2) {
			// This vertex is common to the current half-edge and the previous one
			// This means that h->e->v2 is the targeted vertex
			// We are going to save it in the variable that doesn't contain h->e->v1
			if (v1 != h->e->v1) {
				v1 = h->e->v2;
				v = v1;
			} else {
				v2 = h->e->v2;
				v = v2;
			}
		} else {
			// h->e->v1 is the next new vertex, which means that h->e->v2 is the one
			// that is shared by the current and the previous half-edges. So we are
			// going to save h->e->v1 in the variable that doesn't contain h->e->v2
			if (h->e->v2 != v1) {
				v1 = h->e->v1;
				v = v1;
			} else {
				v2 = h->e->v1;
				v = v2;
			}
		}

		// If we have found the leftmost or the rightmost vertex so far,
		// we remember its coordinates and an iterator on the current half-edge
		if (v->pt.x < x_min) {
			x_min = v->pt.x;
			it_l = it_h;
			shift_l = true;
		}
		if (v->pt.x > x_max) {
			x_max = v->pt.x;
			it_r = it_h;
			shift_r = true;
		}
		y_min = MIN(y_min, v->pt.y);
		y_max = MAX(y_max, v->pt.y);
		++it_h;
	}

	if (shift_l) {
		if (++it_l == edges[contour_id].end()) it_l = edges[contour_id].begin();
	}
	if (shift_r) {
		if (++it_r == edges[contour_id].end()) it_r = edges[contour_id].begin();
	}
}

// Record pixels in grid
void Face::update_pixel_grid(int contour_id, Mat & pixel_grid, int pixel_grid_x_0, int pixel_grid_y_1, double x_min, double x_max, double y_min, double y_max, const Point2d & offset, list<HalfEdge *>::iterator it_l, list<HalfEdge *>::iterator it_r, bool verbose)
{
	// The coordinates system: origin (0.0, 0.0) is set at the center of pixel (rows-1, 0), with y-axis pointing upwards.
	// Therefore the pixel (rows-1, 0) spans the area of [-0.5, 0.5) X [-0.5, 0.5). 

	// Deal with precision at the boundary
	x_min = x_min < 0 ? -0.5 : x_min;
	y_min = y_min < 0 ? -0.5 : y_min;

	// Defines a grid covering the integers included in K = [x_min x_max[ x [y_min y_max[
	int x_0 = int(ceil(x_min));
	int y_0 = (y_min == ceil(y_min) ? int(y_min) + 1 : int(ceil(y_min))); // bottom
	int x_1 = (x_max == floor(x_max) ? int(floor(x_max)) - 1 : int(floor(x_max)));
	int y_1 = int(floor(y_max)); // top
	if (y_0 > y_1 || x_0 > x_1) return;

	int g_x = x_1 - x_0 + 1;
	int g_y = y_1 - y_0 + 1;
	Mat subgrid(g_y, g_x, CV_8U, cv::Scalar::all(1));

	if (contour_id == 0) { // clockwise

						   // First half-loop : from left to right
		list<HalfEdge *>::iterator it_h = it_l;
		while (it_h != it_r) {
			Point2d p1 = (*it_h)->e->v1->pt;
			Point2d p2 = (*it_h)->e->v2->pt;
			int xl_min = 0, xl_max = 0;
			if (p1.x < p2.x) {
				xl_min = int(ceil(p1.x));
				xl_max = (p2.x == floor(p2.x) ? int(p2.x) - 1 : int(floor(p2.x)));
			}
			else {
				xl_min = int(ceil(p2.x));
				xl_max = (p1.x == floor(p1.x) ? int(p1.x) - 1 : int(floor(p1.x)));
			}
			// Finds the equation of the line (v1 v2)
			if (fabs(p2.x - p1.x) > 1e-6) {
				double a = (p2.y - p1.y) / (p2.x - p1.x);
				double b = p1.y - a * p1.x;
				// Discards every point above the line (keep points on the line)
				for (int x = xl_min; x <= xl_max; x++) {
					int yl_min = (a * x + b == ceil(a * x + b)) ? int(a * x + b) + 1 : int(ceil(a * x + b));
					int yl_max = y_1;
					for (int y = yl_min; y <= yl_max; y++) {
						subgrid.at<uchar>(y_1 - y, x - x_0) = 1 - subgrid.at<uchar>(y_1 - y, x - x_0);
					}
				}
			}
			else {
				double y12_max = p1.y > p2.y ? p1.y : p2.y;
				int yl_min = (y12_max == ceil(y12_max)) ? int(y12_max) + 1 : int(ceil(y12_max));
				int yl_max = y_1;
				for (int x = xl_min; x <= xl_max; x++) {
					for (int y = yl_min; y <= yl_max; y++) {
						subgrid.at<uchar>(y_1 - y, x - x_0) = 1 - subgrid.at<uchar>(y_1 - y, x - x_0);
					}
				}
			}
			// Iterates
			if (++it_h == edges[contour_id].end()) it_h = edges[contour_id].begin();
		}

		// Second half-loop : from right to left
		// Same algorithm as before, except that this time, we discard pixels below (v1 v2)
		it_h = it_r;
		while (it_h != it_l) {
			Point2d p1 = (*it_h)->e->v1->pt;
			Point2d p2 = (*it_h)->e->v2->pt;
			int xl_min = 0, xl_max = 0;
			if (p1.x < p2.x) {
				xl_min = int(ceil(p1.x));
				xl_max = (p2.x == floor(p2.x) ? int(p2.x) - 1 : int(floor(p2.x)));
			}
			else {
				xl_min = int(ceil(p2.x));
				xl_max = (p1.x == floor(p1.x) ? int(p1.x) - 1 : int(floor(p1.x)));
			}
			// Finds the equation of the line (v1 v2)
			if (fabs(p2.x - p1.x) > 1e-6) {
				double a = (p2.y - p1.y) / (p2.x - p1.x);
				double b = p1.y - a * p1.x;
				// Discards every point below the line (remove points on the line)
				for (int x = xl_min; x <= xl_max; x++) {
					int yl_min = y_0;
					int yl_max = int(floor(a * x + b));
					for (int y = yl_min; y <= yl_max; y++) {
						subgrid.at<uchar>(y_1 - y, x - x_0) = 1 - subgrid.at<uchar>(y_1 - y, x - x_0);
					}
				}
			}
			else {
				int yl_min = y_0;
				int yl_max = int(floor(max(p1.y, p2.y)));
				for (int x = xl_min; x <= xl_max; x++) {
					for (int y = yl_min; y <= yl_max; y++) {
						subgrid.at<uchar>(y_1 - y, x - x_0) = 1 - subgrid.at<uchar>(y_1 - y, x - x_0);
					}
				}
			}
			// Iterates
			if (++it_h == edges[contour_id].end()) it_h = edges[contour_id].begin();
		}
		// Now copy subgrid to pixel_grid
		for (int y = y_0; y <= y_1; ++y) {
			for (int x = x_0; x <= x_1; ++x) {
				pixel_grid.at<uchar>(pixel_grid_y_1 - y, x - pixel_grid_x_0) = subgrid.at<uchar>(y_1 - y, x - x_0);
			}
		}
	}
	else { // anti-clockwise

		   // First half-loop : from left to right
		list<HalfEdge *>::iterator it_h = it_l;
		while (it_h != it_r) {
			Point2d p1 = (*it_h)->e->v1->pt;
			Point2d p2 = (*it_h)->e->v2->pt;
			int xl_min = 0, xl_max = 0;
			if (p1.x < p2.x) {
				xl_min = int(ceil(p1.x));
				xl_max = (p2.x == floor(p2.x) ? int(p2.x) - 1 : int(floor(p2.x)));
			}
			else {
				xl_min = int(ceil(p2.x));
				xl_max = (p1.x == floor(p1.x) ? int(p1.x) - 1 : int(floor(p1.x)));
			}
			// Finds the equation of the line (v1 v2)
			if (fabs(p2.x - p1.x) > 1e-6) {
				double a = (p2.y - p1.y) / (p2.x - p1.x);
				double b = p1.y - a * p1.x;
				// Discards every point below the line (remove points on the line)
				for (int x = xl_min; x <= xl_max; x++) {
					int yl_min = y_0;
					int yl_max = int(floor(a * x + b));
					for (int y = yl_min; y <= yl_max; y++) {
						subgrid.at<uchar>(y_1 - y, x - x_0) = 1 - subgrid.at<uchar>(y_1 - y, x - x_0);
					}
				}
			}
			else {
				int yl_min = y_0;
				int yl_max = int(floor(max(p1.y, p2.y)));
				for (int x = xl_min; x <= xl_max; x++) {
					for (int y = yl_min; y <= yl_max; y++) {
						subgrid.at<uchar>(y_1 - y, x - x_0) = 1 - subgrid.at<uchar>(y_1 - y, x - x_0);
					}
				}
			}
			// Iterates
			if (++it_h == edges[contour_id].end()) it_h = edges[contour_id].begin();
		}

		// Second half-loop : from right to left
		// Same algorithm as before, except that this time, we discard pixels below (v1 v2)
		it_h = it_r;
		while (it_h != it_l) {
			Point2d p1 = (*it_h)->e->v1->pt;
			Point2d p2 = (*it_h)->e->v2->pt;
			int xl_min = 0, xl_max = 0;
			if (p1.x < p2.x) {
				xl_min = int(ceil(p1.x));
				xl_max = (p2.x == floor(p2.x) ? int(p2.x) - 1 : int(floor(p2.x)));
			}
			else {
				xl_min = int(ceil(p2.x));
				xl_max = (p1.x == floor(p1.x) ? int(p1.x) - 1 : int(floor(p1.x)));
			}
			// Finds the equation of the line (v1 v2)
			if (fabs(p2.x - p1.x) > 1e-6) {
				double a = (p2.y - p1.y) / (p2.x - p1.x);
				double b = p1.y - a * p1.x;
				// Discards every point above the line (keep points on the line)
				for (int x = xl_min; x <= xl_max; x++) {
					int yl_min = (a * x + b == ceil(a * x + b)) ? int(a * x + b) + 1 : int(ceil(a * x + b));
					int yl_max = y_1;
					for (int y = yl_min; y <= yl_max; y++) {
						subgrid.at<uchar>(y_1 - y, x - x_0) = 1 - subgrid.at<uchar>(y_1 - y, x - x_0);
					}
				}
			}
			else {
				double y12_max = p1.y > p2.y ? p1.y : p2.y;
				int yl_min = (y12_max == ceil(y12_max)) ? int(y12_max) + 1 : int(ceil(y12_max));
				int yl_max = y_1;
				for (int x = xl_min; x <= xl_max; x++) {
					for (int y = yl_min; y <= yl_max; y++) {
						subgrid.at<uchar>(y_1 - y, x - x_0) = 1 - subgrid.at<uchar>(y_1 - y, x - x_0);
					}
				}
			}
			// Iterates
			if (++it_h == edges[contour_id].end()) it_h = edges[contour_id].begin();
		}

		// Now remove pixels in subgrid from pixel_grid
		for (int y = y_0; y <= y_1; ++y) {
			for (int x = x_0; x <= x_1; ++x) {
				if (subgrid.at<uchar>(y_1 - y, x - x_0) == 1) pixel_grid.at<uchar>(pixel_grid_y_1 - y, x - pixel_grid_x_0) = 0;
			}
		}
	}
	subgrid.release();
}



void Face::print_details()
{
	std::string name = "face_" + std::to_string(id_face) + ".R";
	FILE* file = fopen(name.c_str(), "w");
	for (int id = 0; id < vertices.size(); ++id) {
		vector<double> x;
		vector<double> y;
		vector<int> l;
		for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = vertices[id].begin(); it_v != vertices[id].end(); it_v++) {
			Vertex* v = it_v->first;
			x.push_back(v->pt.x);
			y.push_back(v->pt.y);
			l.push_back(v->id_vertex);
		}
		fprintf(file, "Contour %i\n", id);
		fprintf(file, "x <- c("); for (int i = 0; i < x.size() - 1; i++) fprintf(file, "%lf, ", x[i]); fprintf(file, "%lf)\n", x[x.size() - 1]);
		fprintf(file, "y <- c("); for (int i = 0; i < y.size() - 1; i++) fprintf(file, "%lf, ", y[i]); fprintf(file, "%lf)\n", y[y.size() - 1]);
		fprintf(file, "l <- c("); for (int i = 0; i < l.size() - 1; i++) fprintf(file, "%i, ", l[i]); fprintf(file, "%i)\n", l[l.size() - 1]);
	}
	fclose(file);
	
}


void Face::get_supporting_segments(set<Segment *> & supporting_segments)
{
	supporting_segments.clear();
	for (auto it_c = edges.begin(); it_c != edges.end(); ++it_c) {
		for (list<HalfEdge *>::iterator it_h = it_c->begin(); it_h != it_c->end(); it_h++) {
			Edge* e = (*it_h)->e;
			if (e->type == INNER_EDGE) {
				Inner_Edge* i_e = static_cast<Inner_Edge*>(e);
				i_e->get_supporting_segments(supporting_segments);
			}
		}
	}
}


void Face::get_neighboring_faces(set<Face *> & neighboring_faces)
{
	neighboring_faces.clear();
	for (auto it_c = edges.begin(); it_c != edges.end(); ++it_c) {
		for (list<HalfEdge *>::iterator it_h = it_c->begin(); it_h != it_c->end(); it_h++) {
			Face* neighbor = (*it_h)->opposite()->f;
			if (neighbor != nullptr) {
				neighboring_faces.insert(neighbor);
			}
		}
	}
}


void Face::classify()
{
	for (int id = 0; id < edges.size(); ++id) {
		vector<HalfEdge *> v_edges(edges[id].begin(), edges[id].end());
		vector<pair<Vertex *, Vertex_Type> > v_vertices(vertices[id].begin(), vertices[id].end());
		uint n = uint(v_vertices.size());

		for (uint i = 0; i < vertices[id].size(); i++) {
			Vertex* v = v_vertices[i].first;
			int v_c = v->connectivity();

			// Accesses the previous and next edges
			Edge* e_1 = v_edges[(i + n - 1) % n]->e;
			Edge* e_2 = v_edges[i]->e;

			// Exhibits their respective support segments
			// The idea is that an inner edge is created by extension of one segment or more
			// which is not the case of an outer edge which relies on the image borders and has no support

			if (e_1->type == Edge_Type::INNER_EDGE && e_2->type == Edge_Type::INNER_EDGE) {

				set<Segment *> support_1, support_2;
				Inner_Edge* ie_1 = static_cast<Inner_Edge*>(e_1);
				Inner_Edge* ie_2 = static_cast<Inner_Edge*>(e_2);
				ie_1->get_supporting_segments(support_1);
				ie_2->get_supporting_segments(support_2);

				// v is a corner if and only if s_1 and s_2 contain colinear segments
				bool detected_colinear_segments = false;
				for (set<Segment *>::iterator it_s1 = support_1.begin(); it_s1 != support_1.end(); it_s1++) {
					Segment* s_1 = (*it_s1);
					for (set<Segment *>::iterator it_s2 = support_2.begin(); it_s2 != support_2.end(); it_s2++) {
						Segment* s_2 = (*it_s2);
						if (s_1 == s_2 || (s_1->node_colinear != nullptr && s_1->node_colinear == s_2->node_colinear)) {
							detected_colinear_segments = true;
							break;
						}
					}
					if (detected_colinear_segments) break;
				}

				if (detected_colinear_segments) {
					v_vertices[i].second = (v_c == 2 ? NON_CORNER_BIVALENT : NON_CORNER_TRIVALENT);
				}
				else {
					v_vertices[i].second = CORNER;
				}

			}
			else if (e_1->type == Edge_Type::OUTER_EDGE && e_2->type == Edge_Type::OUTER_EDGE) {

				// Case when edges e_1 and e_2 are borders
				Outer_Edge* oe_1 = static_cast<Outer_Edge*>(e_1);
				Outer_Edge* oe_2 = static_cast<Outer_Edge*>(e_2);

				if (oe_1->boundary == oe_2->boundary) {
					v_vertices[i].second = (v_c == 2 ? NON_CORNER_BIVALENT : NON_CORNER_TRIVALENT);
				}
				else {
					v_vertices[i].second = CORNER;
				}

			}
			else {
				// v is the intersection of an inner and an other segment
				// It must be a cornerfface
				v_vertices[i].second = CORNER;
			}
		}

		vertices[id] = list<pair<Vertex*, Vertex_Type> >(v_vertices.begin(), v_vertices.end());
	}
	
}



void Face::add_non_corner(HalfEdge* h, Vertex* v, HalfEdge* h_1, HalfEdge* h_2, double min_edge_len)
{
	bool added(false);
	for (int id = 0; id < edges.size(); ++id) {
		list<pair<Vertex*, Vertex_Type> >::iterator it_v = vertices[id].begin();
		list<HalfEdge *>::iterator it_h = edges[id].begin();

		// Localizes the edge to replace
		int edge_id = 0;
		while (it_h != edges[id].end() && (*it_h) != h) {
			++it_v;
			++it_h;
			++edge_id;
		}
		
		if (it_h != edges[id].end()) {
			++it_v;

			// Replaces h = (a b) by h_1 = (a v) and h_2 = (v b)
			it_h = edges[id].erase(it_h);

			edges[id].insert(it_h, h_1);
			edges[id].insert(it_h, h_2);
			h_1->f = this;
			h_2->f = this;

			// Inserts v
			vertices[id].insert(it_v, std::make_pair(v, NON_CORNER_TRIVALENT));
			v->add(this);


			// Adjust edge index for undocked vertices on facet edges
			map<Vertex*, pair<int, int>>::iterator it_v = split_vertex_outer_edge_location.begin();
			while (it_v != split_vertex_outer_edge_location.end()) {
				int & intersect_contour_id = it_v->second.first;
				int & intersect_edge_id = it_v->second.second;
				
				bool remove = false;
				if (intersect_contour_id == id) {
					if (intersect_edge_id > edge_id) {
						++intersect_edge_id;
					}
					else if (intersect_edge_id == edge_id) {
						Vertex* B = h_2->v1_v2 ? h_2->e->v2 : h_2->e->v1;
						// Check if intersection point coincide with the newly added vertex v
						Point2d p_it = it_v->first->pt;
						if (it_v->first == v) {
							remove = true;
							it_v = split_vertex_outer_edge_location.erase(it_v);
						} 
						else if (it_v->first->v_prev == nullptr && it_v->first->v_next == nullptr && cv::norm(p_it - v->pt) < min_edge_len) {
							remove = true;
							Vertex* v_to_delete = it_v->first;
							for (int i = 0; i < split_vertex_chains.size(); ++i) {
								if (split_vertex_chains[i][0] == v_to_delete) {
									Vertex* v_next = split_vertex_chains[i][1];
									// update split vertices list and split edge property map
									split_edges_property[std::make_pair(v, v_next)] = split_edges_property.at(std::make_pair(v_to_delete, v_next));
									split_edges_property.erase(std::make_pair(v_to_delete, v_next));
									split_vertex_chains[i][0] = v;
								}
								else if (split_vertex_chains[i].back() == v_to_delete) {
									int n_chainsize = split_vertex_chains[i].size();
									Vertex* v_prev = split_vertex_chains[i][n_chainsize-2];
									// update split vertices list and split edge property map
									split_edges_property[std::make_pair(v_prev, v)] = split_edges_property.at(std::make_pair(v_prev, v_to_delete));
									split_edges_property.erase(std::make_pair(v_prev, v_to_delete));
									split_vertex_chains[i].back() = v;
								}
							}
							// delete old vertex and remove from map
							delete v_to_delete;
							v_to_delete = nullptr;
							it_v = split_vertex_outer_edge_location.erase(it_v);
						}
						else {
							// Check if intersection point lies between the later edge h_2 = (A, B)
							if (fabs(v->pt.x - B->pt.x) > 1e-3) {
								bool between_vB = (((p_it.x > B->pt.x) && (p_it.x < v->pt.x)) || ((p_it.x > v->pt.x) && (p_it.x < B->pt.x)));
								if (between_vB) ++intersect_edge_id;
							}
							else {
								bool between_vB = (((p_it.y > B->pt.y) && (p_it.y < v->pt.y)) || ((p_it.y > v->pt.y) && (p_it.y < B->pt.y)));
								if (between_vB) ++intersect_edge_id;
							}
						}
					}
				}
				if (!remove) ++it_v;
			}

		}
	}
	
}



void Face::remove_non_corners(Vertex* v1, Edge* e1, Vertex* v2, HalfEdge* h)
{
	// We assume that v1 and v2 are included in the list of vertices of the facet.
	// We assume that all vertices between v1 and v2 are non corner, bivalent vertices that we are going to remove.
	// Finally, the sequence of edges and vertices between v1 and v2 is replaced with a single halfedge h.
	
	for (int id = 0; id < edges.size(); ++id) {
		list<pair<Vertex*, Vertex_Type> >::iterator it_v = vertices[id].begin();
		list<HalfEdge *>::iterator it_h = edges[id].begin();
		int n_edges_initial = edges[id].size();
		int v1_pos_initial = 0;
		int n_removed = 0;
		set<Vertex*> vertices_between = set<Vertex*>();
		bool insert_last_position = false;

		// First of all we have to find v1 at the right location (Note that v1 can appear multiple times)
		HalfEdge* h1 = e1->v1_v2->f == this ? e1->v1_v2 : e1->v2_v1;
		while (*it_h != h1 && it_h != edges[id].end()) {
			++it_v; ++it_h; ++v1_pos_initial;
		}

		if (it_h == edges[id].end()) continue;

		list<pair<Vertex*, Vertex_Type> >::iterator it_v1 = it_v;
		// it_v is now pointing to v1, and it_h to the halfedge (v1 w)

		Vertex* w = (*it_h)->v1_v2 ? (*it_h)->e->v2 : (*it_h)->e->v1;
		it_h = edges[id].erase(it_h);
		++n_removed;
		++it_v;
		// it_v now pointing to w, and it_h to the halfedge (w w_next)
		if (it_v == vertices[id].end() && it_h == edges[id].end()) {
			insert_last_position = true;
			it_h = edges[id].begin();
			it_v = vertices[id].begin();
		}
		else if ((*it_h)->e->v1 != w && (*it_h)->e->v2 != w) {
			std::cout << "Error!" << std::endl;
		}
		// it_v is now pointing to w, and it_h to the halfedge (w w_next)

		// We start removing all the vertices pointing by it_v, starting from w, until finding v2
		while (it_v->first != v2) {
			w = (*it_h)->v1_v2 ? (*it_h)->e->v2 : (*it_h)->e->v1;
			it_h = edges[id].erase(it_h);
			++n_removed;
			vertices_between.insert(it_v->first);
			it_v->first->remove(this);
			it_v = vertices[id].erase(it_v);
			if (it_v == vertices[id].end() && it_h == edges[id].end()) {
				insert_last_position = true;
				it_h = edges[id].begin();
				it_v = vertices[id].begin();
			}
			else if ((*it_h)->e->v1 != w && (*it_h)->e->v2 != w) {
				std::cout << "Error!" << std::endl;
			}
		}

		// In the end, we insert h
		// We need to find the edge iterator that was originally occupied by e1
		list<HalfEdge *>::iterator it_h_insert = edges[id].begin();
		it_v = vertices[id].begin();
		while (it_v != it_v1) {
			++it_v; ++it_h_insert;
		}
		edges[id].insert(it_h_insert, h);

		h->f = this;
		v1->add(this);
		v2->add(this);

		// Adjust edge index for undocked vertices on facet edges
		int v1_shiftback = v1_pos_initial + n_removed > n_edges_initial ? v1_pos_initial + n_removed - n_edges_initial : 0;
		for (map<Vertex*, pair<int, int>>::iterator it_v = split_vertex_outer_edge_location.begin(); it_v != split_vertex_outer_edge_location.end(); ++it_v) {
			Vertex* v_undocked = it_v->first;
			int & intersect_contour_id = it_v->second.first;
			int & intersect_edge_id = it_v->second.second;

			if (intersect_contour_id == id) {
				if (v1_pos_initial + n_removed <= n_edges_initial) {
					if (intersect_edge_id >= v1_pos_initial && intersect_edge_id < v1_pos_initial + n_removed) intersect_edge_id = v1_pos_initial - v1_shiftback;
					else if (intersect_edge_id >= v1_pos_initial + n_removed) intersect_edge_id -= n_removed - 1;
				}
				else {
					if (intersect_edge_id >= v1_pos_initial || intersect_edge_id < (v1_pos_initial + n_removed) % n_edges_initial) intersect_edge_id = v1_pos_initial - v1_shiftback;
					else intersect_edge_id -= v1_shiftback;
				}
			}
		}
		
		// update split vertices list and split edge property map if the vertex is to be deleted, add to the edge index map at the location of the newly added edge
		for (int i = 0; i < split_vertex_chains.size(); ++i) {
			int n_chainsize = split_vertex_chains[i].size();
			if (vertices_between.count(split_vertex_chains[i][0]) > 0) {
				// If the end vertex of the chain is to be deleted, we create a duplicate to replace the old pointer.
				Vertex* v_old = split_vertex_chains[i][0];
				Vertex* v_copy = new Vertex(v_old->id_vertex, v_old->pt.x, v_old->pt.y);
				split_vertex_chains[i][0] = v_copy;
				Vertex* v_next = split_vertex_chains[i][1];
				split_edges_property[std::make_pair(v_copy, v_next)] = split_edges_property.at(std::make_pair(v_old, v_next));
				split_edges_property.erase(std::make_pair(v_old, v_next));
				split_vertex_outer_edge_location.erase(v_old);
				split_vertex_outer_edge_location[v_copy] = std::make_pair(id, v1_pos_initial - v1_shiftback);
			}
			if (vertices_between.count(split_vertex_chains[i].back()) > 0) {
				Vertex* v_old = split_vertex_chains[i].back();
				Vertex* v_copy = new Vertex(v_old->id_vertex, v_old->pt.x, v_old->pt.y);
				split_vertex_chains[i].back() = v_copy;
				Vertex* v_prev = split_vertex_chains[i][n_chainsize - 2];
				split_edges_property[std::make_pair(v_prev, v_copy)] = split_edges_property.at(std::make_pair(v_prev, v_old));
				split_edges_property.erase(std::make_pair(v_prev, v_old));
				split_vertex_outer_edge_location.erase(v_old);
				split_vertex_outer_edge_location[v_copy] = std::make_pair(id, v1_pos_initial - v1_shiftback);
			}
		}

		// safety check
		{
			list<pair<Vertex*, Vertex_Type> >::iterator it_vt = vertices[id].begin();
			list<HalfEdge *>::iterator it_ht = edges[id].begin();
			while (it_vt != vertices[id].end() && it_ht != edges[id].end()) {
				Vertex* v = it_vt->first;
				HalfEdge* h = (*it_ht);
				if (h->v1_v2) {
					assert(h->e->v1 == v);
				}
				else {
					assert(h->e->v2 == v);
				}
				++it_vt;
				++it_ht;
			}
		}
	}

}


void Face::set_as_bivalent(Vertex* v)
{
	bool found_vertex = false;
	typedef pair<Vertex *, Vertex_Type> Vertex_Element;
	for (int id = 0; id < vertices.size(); ++id) {
		for (list<Vertex_Element>::iterator it_v = vertices[id].begin(); it_v != vertices[id].end(); it_v++) {
			if (it_v->first == v) {
				it_v->second = NON_CORNER_BIVALENT;
				found_vertex = true;
				break;
			}
		}
	}
	
	assert(found_vertex);
}


void Face::set_as_trivalent(Vertex* v)
{
	bool found_vertex = false;
	typedef pair<Vertex *, Vertex_Type> Vertex_Element;
	for (int id = 0; id < vertices.size(); ++id) {
		for (list<Vertex_Element>::iterator it_v = vertices[id].begin(); it_v != vertices[id].end(); it_v++) {
			if (it_v->first == v) {
				it_v->second = NON_CORNER_TRIVALENT;
				found_vertex = true;
				break;
			}
		}
	}
	
	assert(found_vertex);
}


void Face::merged_halfedges(Face* f_1, Face* f_2, vector<HalfEdges_Location> & l_12, vector<HalfEdges_Location> & l_21, vector<list<HalfEdge *>> & merged)
{
	// We want to obtain a list of edges of a merged facet (f_1 U f_2) in such a way that the order of the edges form clockwise nested polygon for the outer contour and anti-clockwise polygons for inner contours, 
	// given a sequence of edges to remove in (f_1 U f_2).
	// Here, the sequence of halfedges are given with respect to f_1.
	merged.clear();

	vector<list<HalfEdge *>> merged_contours = vector<list<HalfEdge *>>();
	HalfEdges_Location hl_12, hl_21;

	for (int i_12 = 0; i_12 < l_12.size(); ++i_12) {
		merged_contours.clear();
		int id_1 = l_12[i_12].first, id_2 = -1;
		hl_12 = l_12[i_12];
		vector<HalfEdge *> edges_1(f_1->edges[id_1].begin(), f_1->edges[id_1].end());
		HalfEdge* h = edges_1[hl_12.second.front().first];
		for (int i_21 = 0; i_21 < l_21.size(); ++i_21) {
			int j = l_21[i_21].first;
			for (list<HalfEdge *>::iterator it_h2 = f_2->edges[j].begin(); it_h2 != f_2->edges[j].end(); ++it_h2) {
				if ((*it_h2)->e == h->e) {
					id_2 = j;
					hl_21 = l_21[i_21];
					break;
				}
			}
			if (id_2 != -1) break;
		}
		vector<HalfEdge *> edges_2(f_2->edges[id_2].begin(), f_2->edges[id_2].end());

		// the number of halfedges in the merged cell
		int n_12 = edges_1.size() + edges_2.size();
		for (list<pair<int, int> >::iterator it_h1 = hl_12.second.begin(); it_h1 != hl_12.second.end(); ++it_h1) {
			n_12 -= 2 * (it_h1->second);
		}

		// skip if the merged contour is to be empty
		if (n_12 == 0) continue;
		// form contours if the merging shall return non-empty contours
		int i, j;
		bool start_with_1 = TRUE;
		int length = 0;
		while (length < n_12) {
			if (start_with_1) {
				i = init_start_halfedge(edges_1, hl_12.second, merged_contours);
				if (i < edges_1.size()) {
					merged_contour(i, f_1, f_2, edges_1, edges_2, hl_12, hl_21, merged_contours);
				}
				else {
					start_with_1 = FALSE;
				}
			}
			else {
				j = init_start_halfedge(edges_2, hl_21.second, merged_contours);
				if (j < edges_2.size()) {
					merged_contour(j, f_2, f_1, edges_2, edges_1, hl_21, hl_12, merged_contours);
				}
			}
			length = 0;
			for (int id = 0; id < merged_contours.size(); ++id) {
				length += merged_contours[id].size();
			}
			// std::cout << "halfedges in merged contour = " << length << "/" << n_12 << std::endl;
		}

		if (length != n_12) {
			std::cout << "ERROR!!!" << std::endl;
			merged_contours.clear();
			int i, j;
			bool start_with_1 = TRUE;
			int length = 0;
			while (length < n_12) {
				if (start_with_1) {
					i = init_start_halfedge(edges_1, hl_12.second, merged_contours);
					if (i < edges_1.size()) {
						merged_contour(i, f_1, f_2, edges_1, edges_2, hl_12, hl_21, merged_contours);
					}
					else {
						start_with_1 = FALSE;
					}
				}
				else {
					j = init_start_halfedge(edges_2, hl_21.second, merged_contours);
					if (j < edges_2.size()) {
						merged_contour(j, f_2, f_1, edges_2, edges_1, hl_21, hl_12, merged_contours);
					}
					else {
						break;
					}
				}
				length = 0;
				for (int id = 0; id < merged_contours.size(); ++id) {
					length += merged_contours[id].size();
				}
				std::cout << "halfedges in merged contour = " << length << "/" << n_12 << std::endl;
			}
		}

		if (merged_contours.size() > 0) {
			if (id_1 == 0 && id_2 == 0) {
				merged.insert(merged.begin(), merged_contours.begin(), merged_contours.end());
			}
			else {
				merged.insert(merged.end(), merged_contours.begin(), merged_contours.end());
			}
		}
	}
	
	// add non-overlapping contours
	if (f_1->edges.size() > l_12.size()) {
		for (int id = 0; id < f_1->edges.size(); ++id) {
			bool to_add(true);
			for (int i_12 = 0; i_12 < l_12.size(); ++i_12) {
				if (id == l_12[i_12].first) {
					to_add = false;
					break;
				}
			}
			if (to_add) {
				if (id == 0) {
					merged.insert(merged.begin(), f_1->edges[id]);
				}
				else {
					merged.push_back(f_1->edges[id]);
				}
			}
		}
	}
	if (f_2->edges.size() > l_21.size()) {
		for (int id = 0; id < f_2->edges.size(); ++id) {
			bool to_add(true);
			for (int i_21 = 0; i_21 < l_21.size(); ++i_21) {
				if (id == l_21[i_21].first) {
					to_add = false;
					break;
				}
			}
			if (to_add) {
				if (id == 0) {
					merged.insert(merged.begin(), f_2->edges[id]);
				}
				else {
					merged.push_back(f_2->edges[id]);
				}
			}
		}
	}

}


int Face::init_start_halfedge(vector<HalfEdge *> & edges, list<pair<int, int>> & h, vector<list<HalfEdge *>> & merged_contours)
{
	uint n = edges.size();
	// initialize starting position
	int start_pos = 0;
	bool to_remove = FALSE;
	// choose a starting halfedge that is not already included in merged
	while (start_pos < n) {
		bool found(false);
		for (int id = 0; id < merged_contours.size(); ++id) {
			list<HalfEdge *>::iterator it_m = merged_contours[id].begin();
			while (it_m != merged_contours[id].end()) {
				if (*it_m == edges[start_pos]) break;
				++it_m;
			}
			found = (it_m != merged_contours[id].end());
			if (found) break;
		}

		if (!found) {
			// check if the halfedge is on the removal list
			for (list<pair<int, int> >::iterator it_h = h.begin(); it_h != h.end(); ++it_h) {
				if (it_h->first + it_h->second > n) {
					to_remove = !((start_pos < it_h->first) && (start_pos >= (it_h->first + it_h->second) % n));
				}
				else {
					to_remove = ((start_pos >= it_h->first) && (start_pos < it_h->first + it_h->second));
				}
				if (to_remove) {
					start_pos = max( int((it_h->first + it_h->second) % n) , start_pos+1);
					break;
				}
			}
			if (!to_remove) {
				return start_pos;
			}
		}
		else {
			++start_pos;
		}
	}
	return start_pos;
}



void Face::merged_contour(int start_pos, Face* f_1, Face* f_2, vector<HalfEdge *> & edges_1, vector<HalfEdge *> & edges_2, HalfEdges_Location & hl_12, HalfEdges_Location & hl_21, vector<list<HalfEdge *>> & merged_contours)
{
	list<HalfEdge *> halfedges_12 = list<HalfEdge *>();

	int id_1 = hl_12.first, id_2 = hl_21.first;
	list<pair<int, int>> & h_12 = hl_12.second; 
	list<pair<int, int>> & h_21 = hl_21.second;
	vector<pair<Vertex*, Vertex_Type>> vertices_1(f_1->vertices[id_1].begin(), f_1->vertices[id_1].end());
	vector<pair<Vertex*, Vertex_Type>> vertices_2(f_2->vertices[id_2].begin(), f_2->vertices[id_2].end());
	int n_1 = edges_1.size();
	int n_2 = edges_2.size();


	int i, j;
	i = start_pos;
	halfedges_12.push_back(edges_1[i]);
	i = (i + 1) % n_1;
	while (i != start_pos) {
		bool invalid_jump = FALSE;
		bool jump_to_2 = FALSE;
		for (list<pair<int, int> >::iterator it_h1 = h_12.begin(); it_h1 != h_12.end(); ++it_h1) {
			if (i == it_h1->first) {
				jump_to_2 = TRUE;
				vector<pair<double, HalfEdge *>> v1_dirs, v2_dirs;
				v1_dirs = vertices_1[i].first->directions;

				for (vector<pair<double, HalfEdge *>>::iterator it_d = v1_dirs.begin(); it_d != v1_dirs.end(); ++it_d) {
					if (it_d->second->f == f_2) {
						HalfEdge *h2 = it_d->second;
						// check if the jump is valid, i.e. not jumping to the a halfedge that is already in the to-be-added list, or already in the merged list, or an overlapping halfedge to be removed
						list<HalfEdge *>::iterator it_h = halfedges_12.begin();
						while (it_h != halfedges_12.end()) {
							if (*it_h == h2) {
								invalid_jump = true;
								break;
							}
							++it_h;
						}

						for (int id = 0; id < merged_contours.size(); ++id) {
							if (invalid_jump) break;
							list<HalfEdge *>::iterator it_m = merged_contours[id].begin();
							while (it_m != merged_contours[id].end()) {
								if (*it_m == h2) {
									invalid_jump = true;
									break;
								}
								++it_m;
							}
						}

						if (!invalid_jump) {
							j = 0;
							// locate the halfedge starting at vertices_1[i]
							while (edges_2[j] != h2) ++j;
							for (list<pair<int, int> >::iterator it_h2 = h_21.begin(); it_h2 != h_21.end(); ++it_h2) {
								if (it_h2->first + it_h2->second > n_2) {
									invalid_jump = !((j < it_h2->first) && (j >= (it_h2->first + it_h2->second) % n_2));
								}
								else {
									invalid_jump = ((j >= it_h2->first) && (j < it_h2->first + it_h2->second));
								}
								if (invalid_jump) {
									jump_to_2 = FALSE;
									break;
								}
							}
						}

						if (invalid_jump) {
							jump_to_2 = FALSE;
							break;
						}
						else {
							// if found a valid jump onto f_2
							halfedges_12.push_back(edges_2[j]);
							j = (j + 1) % n_2;
							// loop through f_2 until meeting a halfedge to be removed
							while (edges_2[j] != h2) {
								for (list<pair<int, int> >::iterator it_h2 = h_21.begin(); it_h2 != h_21.end(); ++it_h2) {
									if (j == it_h2->first) {
										v2_dirs = vertices_2[j].first->directions;
										jump_to_2 = FALSE;
									}
								}
								if (jump_to_2) {
									halfedges_12.push_back(edges_2[j]);
									j = (j + 1) % n_2;
								}
								else {
									break;
								}
							}
						}
						
					}
				}
				
				if (invalid_jump) {
					// no valid jump to f_2 exists
					Vertex* vi_jump = nullptr;
					for (list<pair<int, int> >::iterator it_h12 = h_12.begin(); it_h12 != h_12.end(); ++it_h12) {
						int i_jump = (it_h12->first + it_h12->second) % n_1;
						vi_jump = edges_1[i_jump]->v1_v2 ? edges_1[i_jump]->e->v1 : edges_1[i_jump]->e->v2;
						if (vi_jump == vertices_1[i].first) {
							i = i_jump; break;
						}
					}
					assert(vi_jump == vertices_1[i].first);
				}
				else {
					// relocate position of i to jump back to f_1 after looping through part of f_2
					if (!jump_to_2) {
						double alpha_d = 100;
						for (vector<pair<double, HalfEdge *>>::iterator it_d = v2_dirs.begin(); it_d != v2_dirs.end(); ++it_d) {
							if (it_d->second->f == f_1  && it_d->first <= alpha_d) {
								HalfEdge *h1 = it_d->second;
								int i_jump = 0;
								while (edges_1[i_jump] != h1) i_jump = (i_jump + 1) % n_1;
								// check if the jump is valid, i.e. not jumping to a halfedge to be removed (inclusion check is later)
								for (list<pair<int, int> >::iterator it_h12 = h_12.begin(); it_h12 != h_12.end(); ++it_h12) {
									if (it_h12->first + it_h12->second > n_1) {
										invalid_jump = !((i_jump < it_h12->first) && (i_jump >= (it_h12->first + it_h12->second) % n_1));
									}
									else {
										invalid_jump = ((i_jump >= it_h12->first) && (i_jump <  it_h12->first + it_h12->second));
									}
								}
								if (!invalid_jump) {
									i = i_jump; alpha_d = it_d->first;
								}
									
							}
						}
					}
				}

			}
		}
		bool found(false);
		list<HalfEdge *>::iterator it_h12 = halfedges_12.begin();
		while (it_h12 != halfedges_12.end()) {
			if (*it_h12 == edges_1[i]) {
				found = true;
				break;
			}
			++it_h12;
		}
		for (int id = 0; id < merged_contours.size(); ++id) {
			if (found) break;
			list<HalfEdge *>::iterator it_m = merged_contours[id].begin();
			while (it_m != merged_contours[id].end()) {
				if (*it_m == edges_1[i]) {
					found = true;
					break;
				}
				++it_m;
			}
		}
		if (found) break;
		if (!jump_to_2) {
			halfedges_12.push_back(edges_1[i]);
		}
		i = (i + 1) % n_1;
	}

	// get the order of the contour (shoelace formula)
	int clockwise = HalfEdge::is_clockwise(halfedges_12);
	if (clockwise == 1) {
		merged_contours.insert(merged_contours.begin(), halfedges_12);
	}
	else {
		merged_contours.push_back(halfedges_12);
	}
}