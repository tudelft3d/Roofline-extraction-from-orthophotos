#include "segment_ray.h"
#include "defs.h"
#include "means.h"


Segment::Segment(uint _index, double x1, double y1, double x2, double y2, double _width, double _log_nfa, bool _is_artificial, bool _is_color_edge)
: support_cost(0.)
{
    index = _index;

    // Sets the two endpoints of the segment and computes the current direction
    // of the segment such that alpha = mes(direction) \in [0, pi]
    end1 = Point2d(x1, y1);
    end2 = Point2d(x2, y2);
    barycenter = Point2d((x1 + x2) / 2, (y1 + y2) / 2);
    length = cv::norm(end2 - end1);
	width = _width;

    referencing_coordinates = barycenter;

    direction = (end2 - end1)/length;
    if (direction[1] < 0 || (direction[1] == 0 && direction[0] < 0)) direction = -direction;

    alpha = atan2(direction[1], direction[0]) * 180 / PI;
	if (alpha < 0) alpha += 180;

    // Current coordinates are equel to the initial coordinates
    finalBarycenter = interBarycenter = barycenter;
    finalEnd1 = interEnd1 = end1;
    finalEnd2 = interEnd2 = end2;
    finalDirection = direction;

    // For the moment this segment is not assigned to any cluster of parallel or colinear segments
    node_parallel = NULL;
    node_colinear = NULL;

    // No perturbation for the moment
    is_dalpha_set = false;
    dalpha = 0;
	theta = alpha;

	a = -sin(theta * PI / 180);
	b = cos(theta * PI / 180);
	c = -a * finalBarycenter.x - b * finalBarycenter.y;

    is_dt_set = false;
    dt = 0;

    is_artificial = _is_artificial;

	is_color_edge = _is_color_edge;
    // By default the segment is not disabled, and it doesn't produce rays for the moment
    is_disabled = false;
    rays = std::make_pair(reinterpret_cast<Ray*>(NULL), reinterpret_cast<Ray*>(NULL));
}

Segment::Segment(uint _index, double x1, double y1, double x2, double y2, double _width, double _support_cost, Point2d _grad_normal, bool _is_artificial, bool _is_color_edge)
{
	index = _index;

	// Sets the two endpoints of the segment and computes the current direction
	// of the segment such that alpha = mes(direction) \in [0, pi]
	end1 = Point2d(x1, y1);
	end2 = Point2d(x2, y2);
	barycenter = Point2d((x1 + x2) / 2, (y1 + y2) / 2);
	length = cv::norm(end2 - end1);
	width = _width;
	referencing_coordinates = barycenter;

	direction = (end2 - end1) / length;
	if (direction[1] < 0 || (direction[1] == 0 && direction[0] < 0)) direction = -direction;

	alpha = atan2(direction[1], direction[0]) * 180 / PI;
	if (alpha < 0) alpha += 180;

	// Current coordinates are equel to the initial coordinates
	finalBarycenter = interBarycenter = barycenter;
	finalEnd1 = interEnd1 = end1;
	finalEnd2 = interEnd2 = end2;
	finalDirection = direction;

	// For the moment this segment is not assigned to any cluster of parallel or colinear segments
	node_parallel = NULL;
	node_colinear = NULL;

	// No perturbation for the moment
	is_dalpha_set = false;
	dalpha = 0;
	theta = alpha;

	a = -sin(theta * PI / 180);
	b = cos(theta * PI / 180);
	c = -a * finalBarycenter.x - b * finalBarycenter.y;

	is_dt_set = false;
	dt = 0;
	is_artificial = _is_artificial;

	support_cost = _support_cost;
	grad_normal = _grad_normal;
	is_color_edge = _is_color_edge;


	// By default the segment is not disabled, and it doesn't produce rays for the moment
	is_disabled = false;
	rays = std::make_pair(reinterpret_cast<Ray*>(NULL), reinterpret_cast<Ray*>(NULL));
}


Segment::~Segment()
{
}


void Segment::set_index(const uint _index)
{
    index = _index;
}

/*
void Segment::set_dalpha(double _dalpha)
{
    // We rotate the segment anti-clockwise about its barycenter
}
*/

void Segment::set_dalpha(double _dalpha, double _theta, double _a, double _b, double _c, Vec2d &_u)
{
    // We rotate the segment by an angle dalpha around its barycenter
    dalpha = _dalpha;
	theta = _theta;

    // Update the supporting line equation
    a = _a, b = _b, c = _c;

    // The position of the barycenter currentBarycenter remains unchanged
    // which is not the case of the ends. The abscissa of the points are
    // computed according to the new angle dalpha, however we want to make
    // sure that the final points belong to the line defined by equation
    // _ax + _by + _c = 0, so we are going to use this equation to find
    // the ordinates of the points
    interBarycenter = barycenter;
    finalDirection = _u;
    if (finalDirection[1] < 0 || (finalDirection[1] == 0 && finalDirection[0] < 0)) {
        finalDirection = -finalDirection;
    }

    if (fabs(finalDirection[0]) > fabs(finalDirection[1])) {
        interEnd1.x = interBarycenter.x - length * finalDirection[0] / 2;
        interEnd2.x = interBarycenter.x + length * finalDirection[0] / 2;
        interEnd1.y = (-c - a * interEnd1.x) / b;
        interEnd2.y = (-c - a * interEnd2.x) / b;
    } else {
        interEnd1.y = interBarycenter.y - length * finalDirection[1] / 2;
        interEnd2.y = interBarycenter.y + length * finalDirection[1] / 2;
        interEnd1.x = (-c - b * interEnd1.y) / a;
        interEnd2.x = (-c - b * interEnd2.y) / a;
    }

    finalBarycenter = interBarycenter;
    finalEnd1 = interEnd1;
    finalEnd2 = interEnd2;
    length = norm(finalEnd2 - finalEnd1);

	is_dalpha_set = true;
}


void Segment::set_dt(double _dt)
{
    // We translate the segment by _dt in the direction of the normal vector
    dt = _dt;

    Vec2d finalNormal = Vec2d(-finalDirection[1], finalDirection[0]);
    finalEnd1 = Point2d(interEnd1.x + dt * finalNormal[0], interEnd1.y + dt * finalNormal[1]);
    finalEnd2 = Point2d(interEnd2.x + dt * finalNormal[0], interEnd2.y + dt * finalNormal[1]);
    finalBarycenter = (finalEnd1 + finalEnd2) / 2;

    c = -a * finalBarycenter.x - b * finalBarycenter.y;

    is_dt_set = true;
}


void Segment::set_dt(double _dt, double _a, double _b, double _c, Vec2d & _u)
{
    // We translate the segment by _dt in the direction of the normal vector
    dt = _dt;

    // We update the equation of the support line
    a = _a, b = _b, c = _c;
    finalDirection = _u;
    if (finalDirection[1] < 0 || (finalDirection[1] == 0 && finalDirection[0] < 0)) {
        finalDirection = -finalDirection;
    }

    Vec2d finalNormal = Vec2d(-finalDirection[1], finalDirection[0]);
    if (fabs(b) > fabs(a)) {
        finalBarycenter.x = interBarycenter.x + dt * finalNormal[0];
        finalEnd1.x = interEnd1.x + dt * finalNormal[0];
        finalEnd2.x = interEnd2.x + dt * finalNormal[0];
        finalBarycenter.y = (-c - a * finalBarycenter.x) / b;
        finalEnd1.y = (-c - a * finalEnd1.x) / b;
        finalEnd2.y = (-c - a * finalEnd2.x) / b;
    } else {
        finalBarycenter.y = interBarycenter.y + dt * finalNormal[1];
        finalEnd1.y = interEnd1.y + dt * finalNormal[1];
        finalEnd2.y = interEnd2.y + dt * finalNormal[1];
        finalBarycenter.x = (-c - b * finalBarycenter.y) / a;
        finalEnd1.x = (-c - b * finalEnd1.y) / a;
        finalEnd2.x = (-c - b * finalEnd2.y) / a;
    }

    is_dt_set = true;
}

void Segment::add_ray(Ray *r)
{
    if (rays.first == NULL) {
        rays.first = r;
    } else if (rays.second == NULL) {
        rays.second = r;
    }
}


void Segment::remove_ray(Ray *r)
{
    if (rays.first == r) {
        rays.first = reinterpret_cast<Ray*>(NULL);
    } else if (rays.second == r) {
        rays.second = reinterpret_cast<Ray*>(NULL);
    }
}


void Segment::update_line_coefficients()
{
    if (finalEnd1.x == finalEnd2.x) {
        a = 1;
        b = 0;
        c = -finalEnd1.x;
    } else {
        double final_alpha = atan2(finalDirection[1], finalDirection[0]);
        a = -sin(final_alpha);
        b = cos(final_alpha);
        c = -a * finalBarycenter.x - b * finalBarycenter.y;
    }
}


void Segment::set_final_extrema(Point2d & _finalEnd1, Point2d & _finalEnd2)
{
	// Arguments are points that belong to the line ax + by + c = 0
	// The direction should remain unchanged
	finalEnd1 = _finalEnd1;
	finalEnd2 = _finalEnd2;

    finalBarycenter = (finalEnd1 + finalEnd2) / 2;
    length = norm(finalEnd2 - finalEnd1);
}


void Segment::move_to_line(double _a, double _b, double _c) {
    a = _a;
    b = _b;
    c = _c;

    if (fabs(b) > fabs(a)) {
        finalBarycenter.y = (-c - a * finalBarycenter.x) / b;
        finalEnd1.y = (-c - a * finalEnd1.x) / b;
        finalEnd2.y = (-c - a * finalEnd2.x) / b;
        length = norm(finalEnd2 - finalEnd1);
        finalDirection = (finalEnd2 - finalEnd1) / length;
    } else {
        finalBarycenter.x = (-c - b * finalBarycenter.y) / a;
        finalEnd1.x = (-c - b * finalEnd1.y) / a;
        finalEnd2.x = (-c - b * finalEnd2.y) / a;
        length = norm(finalEnd2 - finalEnd1);
        finalDirection = (finalEnd2 - finalEnd1) / length;
    }
}


void Segment::bounding_box_coordinates(double & x_min, double & x_max, double & y_min, double & y_max)
{
    if (finalEnd1.x < finalEnd2.x) {
        x_min = finalEnd1.x;
        x_max = finalEnd2.x;
    } else {
        x_min = finalEnd2.x;
        x_max = finalEnd1.x;
    }
    if (finalEnd1.y < finalEnd2.y) {
        y_min = finalEnd1.y;
        y_max = finalEnd2.y;
    } else {
        y_min = finalEnd2.y;
        y_max = finalEnd1.y;
    }
}


void Segment::draw(Matrix<uchar> & background, uchar* rgb)
{
    Point2i p1 = Point2i(int(round(jclamp(0, finalEnd1.x, background.cols - 1))), int(round(jclamp(0, background.rows - finalEnd1.y, background.rows - 1))));
    Point2i p2 = Point2i(int(round(jclamp(0, finalEnd2.x, background.cols - 1))), int(round(jclamp(0, background.rows - finalEnd2.y, background.rows - 1))));
    background.line(uint(p1.y), uint(p1.x), uint(p2.y), uint(p2.x), rgb);
}


void Segment::enable()
{
    is_disabled = false;
}

void Segment::disable()
{
    is_disabled = true;
}

void Segment::enable(vector<Segment *> &T)
{
    for (uint i = 0 ; i < T.size() ; i++) {
        T[i]->enable();
    }
}


void Segment::clear_segments(std::vector<Segment *> & T)
{
    for (int i = 0; i < (int)T.size(); i++) delete T[i];
    T.clear();
}


void Segment::draw_segment(Matrix<uchar> & I, Matrix<uchar> & J, Segment* s)
{
    J = Matrix<uchar>(I.rows, I.cols, 3);
    for (uint i = 0; i < J.rows; i++) {
        for (uint j = 0; j < J.cols; j++) {
            for (uint c = 0; c < 3; c++) J(i, j, c) = I(i, j, c);
        }
    }

    uchar rgb[3] = { 255, 0, 0 };
    Point2i p1 = Point2i(int(round(jclamp(0, s->finalEnd1.x, J.cols - 1))), int(round(jclamp(0, J.rows - s->finalEnd1.y, J.rows - 1))));
    Point2i p2 = Point2i(int(round(jclamp(0, s->finalEnd2.x, J.cols - 1))), int(round(jclamp(0, J.rows - s->finalEnd2.y, J.rows - 1))));
    J.line(uint(p1.y), uint(p1.x), uint(p2.y), uint(p2.x), rgb);
}


bool Segment::print_segment(Segment* s, Matrix<uchar> & background, std::string & name)
{
    Matrix<uchar> J;
    draw_segment(background, J, s);
    return J.write_uchar(name);
}


void Segment::draw_segments(Matrix<uchar> & I, Matrix<uchar> & J, std::vector<Segment *> &T)
{
    J = Matrix<uchar>(I.rows, I.cols, 3);
    for (uint i = 0; i < J.rows; i++) {
        for (uint j = 0; j < J.cols; j++) {
            for (uint c = 0; c < 3; c++) J(i, j, c) = I(i, j, c);
        }
    }

    for (int i = 0; i < int(T.size()); i++) {
        uchar rgb[3] = { 255, 0, 0 };

        Point2i p1 = Point2i(int(round(jclamp(0, T[i]->finalEnd1.x, J.cols - 1))), int(round(jclamp(0, J.rows - T[i]->finalEnd1.y, J.rows - 1))));
        Point2i p2 = Point2i(int(round(jclamp(0, T[i]->finalEnd2.x, J.cols - 1))), int(round(jclamp(0, J.rows - T[i]->finalEnd2.y, J.rows - 1))));
        J.line(uint(p1.y), uint(p1.x), uint(p2.y), uint(p2.x), rgb);

    }
}


bool Segment::print_segments(std::vector<Segment *> &T, Matrix<uchar> & background, std::string & name)
{
    Matrix<uchar> J;
    draw_segments(background, J, T);
    return J.write_uchar(name);
}




Ray::Ray(uint _index, Segment* _parent, Point2d & _O, Point2d & _A, const Vec2d & _OA, double _incidence_angle, uint _ttl)
{
    index = _index;
    parent = _parent;
    parent->add_ray(this);

    O = _O;
    A = _A;
	OA = _OA;
	incidence_angle = _incidence_angle; // is in the opposite direction of O-->A

	double dx = A.x - O.x, dy = A.y - O.y;
    initial_length = sqrt(dx * dx + dy * dy);

    t = 0.0;
    t_set = false;
    ttl = _ttl;

	// Idea : there is a primary condition for stopping a ray that is based either on the number of intersections or the distance ran
	// When this primary condition becauses 'false', we may set secondary condition to 'true' if the gradient is high enough
	// The ray is stopped iff both conditions are set to false
	primary_condition = true;
	secondary_condition = false;
	stopped = false;
	t_swap = FLT_MAX;

	t_collinear_boundaries = list<pair<double, double>>();
}


Ray::~Ray()
{
    parent->remove_ray(this);
}


void Ray::set_time(double _t)
{
    t = _t;
    t_set = true;
}



bool Ray::is_time_set()
{
    return t_set;
}



bool Ray::intersects_orthogonal_ray(Ray* r)
{
    double d = 2.0;

    double alpha_1, alpha_2;
    if (parent->node_parallel != NULL && r->parent->node_parallel != NULL) {
        alpha_1 = parent->node_parallel->angle;
        alpha_2 = r->parent->node_parallel->angle;
    } else {
        alpha_1 = atan2(parent->finalDirection[1], parent->finalDirection[0]);
        alpha_2 = atan2(r->parent->finalDirection[1], r->parent->finalDirection[0]);
    }

    double dalpha = alpha_2 - alpha_1;
    if (dalpha > 0) {
        return (fabs(dalpha - 90) < d);
    } else {
        return (fabs(-dalpha - 90) < d);
    }
}


void Ray::stop()
{
	primary_condition = false;
	secondary_condition = false;
	stopped = true;
}


void Ray::should_ray_be_stopped()
{
	stopped = !primary_condition && !secondary_condition;
}


bool Ray::has_been_stopped()
{
    return stopped;
}


void Ray::build_rays(std::vector<Segment*> & segments, std::vector<Ray*> & rays, unsigned int ttl, double length_thresh)
{
    clear_rays(rays);
    uint index = 0;

	// A segment s_i makes an angle theta \in [0, PI] with the horizontal line
	// This angle theta is used to compute the line equation of the segment, that is ax + by + c = 0
	// with a = -sin(theta) and b = cos(theta). n = (a, b) is orthogonal to the segment and u = (-b, a)
	// gives its direction.

    for (unsigned int i = 0; i < segments.size() ; i++) {
        Segment* s_i = segments[i];
		uint ttl_i = s_i->length > length_thresh ? ttl : 1;
		s_i->rays = std::make_pair(reinterpret_cast<Ray*>(NULL), reinterpret_cast<Ray*>(NULL));
        if (!s_i->is_disabled) {
            Point2d & O = s_i->finalBarycenter;
            Point2d & A = s_i->finalEnd1;
            Point2d & B = s_i->finalEnd2;

			Ray *r_OA = NULL, *r_OB = NULL;

			// We want to determine if its u or -u that we should assign as direction vector to the rays OA and OB.
			// To this end we should check if u.x and AB.x go in the same direction.

			Vec2d u = Vec2d(-s_i->b, s_i->a);
			Vec2d AB = Vec2d(B.x - A.x, B.y - A.y);
			bool same_orientation = (fabs(u[0]) > 1e-6 ? u[0] * AB[0] > 0 : u[1] * AB[1] > 0);
			if (same_orientation) {
				r_OA = new Ray(index + 0, s_i, O, A, -u, (s_i->theta - 180) * PI / 180, ttl_i);
				r_OB = new Ray(index + 1, s_i, O, B, u, s_i->theta * PI / 180, ttl_i);
			} else {
				r_OA = new Ray(index + 0, s_i, O, A, u, s_i->theta * PI / 180, ttl_i);
				r_OB = new Ray(index + 1, s_i, O, B, -u, (s_i->theta - 180) * PI / 180, ttl_i);
			}
            rays.push_back(r_OA);
            rays.push_back(r_OB);
            index += 2;
        }
    }
}


void Ray::clear_rays(std::vector<Ray*> & rays)
{
    for (unsigned int i = 0 ; i < rays.size() ; i++) {
        delete rays[i];
    }
    rays.clear();
}



Node_Colinear_Segments::Node_Colinear_Segments(double _ordinate)
    : ordinate (_ordinate)
{
    colinear_segments.clear();
}


Node_Colinear_Segments::Node_Colinear_Segments(const Node_Colinear_Segments &node)
{
    ordinate = node.ordinate;
    colinear_segments = list<Segment *>(node.colinear_segments.begin(), node.colinear_segments.end());
}


Node_Colinear_Segments::~Node_Colinear_Segments()
{
    delete_references_to_colinear_node();
    colinear_segments.clear();
}


void Node_Colinear_Segments::delete_references_to_colinear_node()
{
    for (list<Segment *>::iterator it_s = colinear_segments.begin() ; it_s != colinear_segments.end() ; it_s++) {
        (*it_s)->node_colinear = NULL;
    }
}

void Node_Colinear_Segments::add(Segment *s)
{
    colinear_segments.push_back(s);
    s->node_colinear = this;
}


void Node_Colinear_Segments::remove(Segment *s)
{
    list<Segment *>::iterator it_s = colinear_segments.begin();
    while (it_s != colinear_segments.end()) {
        if ((*it_s) == s) {
            it_s = colinear_segments.erase(it_s);
            s->node_colinear = NULL;
            break;
        }
        ++it_s;
    }
}




Node_Parallel_Segments::Node_Parallel_Segments(double _angle)
    : angle (_angle)
{
    frame_origin = Point2d();
    colinear_segments = map<double, Node_Colinear_Segments*>();
    other_segments = list<Segment *>();

    parallel_segments = list<Segment *>();
}


Node_Parallel_Segments::Node_Parallel_Segments(const Node_Parallel_Segments &node)
{
    angle = node.angle;
    frame_origin = node.frame_origin;
    colinear_segments = map<double, Node_Colinear_Segments*>(node.colinear_segments.begin(), node.colinear_segments.end());
    other_segments = list<Segment *>(node.other_segments.begin(), node.other_segments.end());

    parallel_segments = list<Segment *>(node.parallel_segments.begin(), node.parallel_segments.end());
}


Node_Parallel_Segments::~Node_Parallel_Segments()
{
    // Deletes references to the current node in all segments that are assigned to it
    delete_references_to_parallel_node();

    // Destroys objects
    for (std::map<double, Node_Colinear_Segments*>::iterator it_m = colinear_segments.begin(); it_m != colinear_segments.end(); it_m++) {
        delete it_m->second;
    }
    colinear_segments.clear();

    other_segments.clear();
    parallel_segments.clear();
}


void Node_Parallel_Segments::delete_references_to_parallel_node()
{
    for (map<double, Node_Colinear_Segments*>::iterator it_m = colinear_segments.begin() ; it_m != colinear_segments.end() ; it_m++) {
        Node_Colinear_Segments* node = it_m->second;
        for (list<Segment *>::iterator it_s = node->colinear_segments.begin() ; it_s != node->colinear_segments.end() ; it_s++) {
            if (*it_s) (*it_s)->node_parallel = NULL;
        }
    }
    for (list<Segment *>::iterator it_s = other_segments.begin() ; it_s != other_segments.end() ; it_s++) {
        if (*it_s) (*it_s)->node_parallel = NULL;
    }
    for (list<Segment *>::iterator it_s = parallel_segments.begin() ; it_s != parallel_segments.end() ; it_s++) {
        if (*it_s) (*it_s)->node_parallel = NULL;
    }
}


void Node_Parallel_Segments::delete_colinear_nodes()
{
    if (parallel_segments.size() == 0) {
        parallel_segments = std::list<Segment *>(other_segments.begin(), other_segments.end());
        other_segments.clear();

        // We delete all colinear nodes and put them again in parallel_segments
        // This way we can run the second step of the regularization process again
        for (map<double, Node_Colinear_Segments *>::iterator it_m = colinear_segments.begin() ; it_m != colinear_segments.end() ; it_m++) {
            Node_Colinear_Segments* node = it_m->second;
            list<Segment *>::iterator it_s = node->colinear_segments.begin();
            while (it_s != node->colinear_segments.end()) {
                (*it_s)->node_colinear = NULL;
                parallel_segments.push_back(*it_s);
                it_s = node->colinear_segments.erase(it_s);
            }
            delete node;
        }
        colinear_segments.clear();
    }
}


void Node_Parallel_Segments::create_colinear_node(double _ordinate)
{
    if (colinear_segments.find(_ordinate) == colinear_segments.end()) {
        colinear_segments[_ordinate] = new Node_Colinear_Segments(_ordinate);
    }
}


void Node_Parallel_Segments::add(Segment *s)
{
    parallel_segments.push_back(s);
    s->node_parallel = this;
}


void Node_Parallel_Segments::remove(Segment *s)
{
    list<Segment *>::iterator it_s = parallel_segments.begin();
    while (it_s != parallel_segments.end()) {
        if ((*it_s) == s) {
            it_s = parallel_segments.erase(it_s);
            s->node_parallel = NULL;
            break;
        }
        ++it_s;
    }
}


void Node_Parallel_Segments::assign_to_colinear_node(double _ordinate, Segment* s)
{
    if (colinear_segments.find(_ordinate) != colinear_segments.end()) {
        colinear_segments[_ordinate]->add(s);
    }
}


void Node_Parallel_Segments::assign_to_colinear_node(double _ordinate, list<Segment *> & ls)
{
    if (colinear_segments.find(_ordinate) != colinear_segments.end()) {
        for (list<Segment *>::iterator it_s = ls.begin(); it_s != ls.end(); it_s++) {
            colinear_segments[_ordinate]->add(*it_s);
        }
    }
}


void Node_Parallel_Segments::assign_to_other(Segment* s)
{
    other_segments.push_back(s);
}


void Node_Parallel_Segments::assign_to_other(list<Segment *> & ls)
{
    for (list<Segment *>::iterator it_s = ls.begin(); it_s != ls.end(); it_s++) {
        other_segments.push_back((*it_s));
    }
}

