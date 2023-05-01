#ifndef SEGMENT_RAY_H
#define SEGMENT_RAY_H

#include <opencv2/core.hpp>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <utility>
#include <string>

#include "quadtree_point.h"
#include "matrix.h"


using std::vector;
using std::list;
using std::map;
using std::set;
using std::pair;

using cv::Point2i;
using cv::Point2d;
using cv::Vec2d;
using cv::Size2i;


typedef enum {
    TOP_IMAGE = -1,
    BOTTOM_IMAGE = -2,
    LEFT_IMAGE = -3,
    RIGHT_IMAGE = -4,
	INVALID_BORDER = -5,
	POLYGON_OUTER = -6,
	POLYGON_INNER = -7
} Image_Boundary;


class Ray;

class Node_Parallel_Segments;

class Node_Colinear_Segments;



class Segment {
public:
    Segment(uint _index, double _x1, double _y1, double _x2, double _y2, double _width, double _log_nfa, bool _is_artificial, bool _is_color_edge);

	Segment(uint _index, double _x1, double _y1, double _x2, double _y2, double _width, double _support_cost, Point2d _normal, bool _is_artificial, bool _is_color_edge);

    ~Segment();

    void set_index(const uint _index);

    void set_dalpha(double _dalpha);

    void set_dalpha(double _dalpha, double _theta, double _a, double _b, double _c, Vec2d & _u);

    void set_dt(double _dt);

    void set_dt(double _dt, double _a, double _b, double _c, Vec2d & _u);

    void add_ray(Ray *r);

    void remove_ray(Ray *r);

    void update_line_coefficients();

    void set_final_extrema(Point2d & _finalEnd1, Point2d & _finalEnd2);

    void move_to_line(double _a, double _b, double _c);

    void bounding_box_coordinates(double & x_min, double & x_max, double & y_min, double & y_max);

    void draw(Matrix<uchar> & background, uchar* color);

    void enable();

    void disable();

    static void enable(vector<Segment *> & T);


    static void clear_segments(vector<Segment *> & T);

    static void draw_segment(Matrix<uchar> & I, Matrix<uchar> & J, Segment* s);
    static void draw_segments(Matrix<uchar> & I, Matrix<uchar> & J, vector<Segment *> &T);

    static bool print_segment(Segment* s, Matrix<uchar> & background, std::string & name);
    static bool print_segments(vector<Segment *> &T, Matrix<uchar> & background, std::string & name);

public:
    int index;

    Point2d barycenter;
    Point2d end1;
    Point2d end2;
    Vec2d direction;
    double length;
	double width;

    double alpha;
	double support_cost;
	vector<Point2i> pixels;
	Point2d grad_normal;
	bool is_color_edge;
	bool to_delete = false;

	Point2d referencing_coordinates;

    Node_Parallel_Segments* node_parallel;
    Node_Colinear_Segments* node_colinear;

    double dalpha;
    bool is_dalpha_set;
	double theta;

    double dt;
    bool is_dt_set;

    double a, b, c;

    Point2d interBarycenter;
    Point2d interEnd1;
    Point2d interEnd2;

    Point2d finalBarycenter;
    Point2d finalEnd1;
    Point2d finalEnd2;
    Vec2d finalDirection;

    bool is_artificial;
    bool is_disabled;
    pair<Ray*, Ray*> rays;
};



class Ray {
public:
    Ray(uint _index, Segment* _parent, Point2d & _O, Point2d & _A, const Vec2d & _OA, double _incidence_angle, uint _ttl);

    ~Ray();

	inline Ray* opposite() {
		return (parent->rays.first == this ? parent->rays.second : parent->rays.first);
	}

    void set_time(double _t);

    bool is_time_set();

    bool intersects_orthogonal_ray(Ray* r);

	void stop();

	void should_ray_be_stopped();

    bool has_been_stopped();

    static void build_rays(std::vector<Segment *> & segments, std::vector<Ray *> & rays, unsigned int ttl, double length_thresh);

    static void clear_rays(std::vector<Ray *> & rays);

public:
    uint index;
    Segment* parent;
    Point2d O;
    Point2d A;
    Vec2d OA;
	double incidence_angle; // is in the opposite direction of O-->A
    double initial_length;

    double t;
    bool t_set;

    int ttl;
	bool primary_condition;
	bool secondary_condition;
	double t_swap;

	set<int> collinear_boundaries;

	list<pair<double, double>> t_collinear_boundaries;

private:
	bool stopped;
};



class Node_Colinear_Segments
{
public:
    Node_Colinear_Segments(double _ordinate);

    Node_Colinear_Segments(const Node_Colinear_Segments & node);

    ~Node_Colinear_Segments();

    void delete_references_to_colinear_node();

    void add(Segment *s);

    void remove(Segment *s);

public:
    double ordinate;
    list<Segment *> colinear_segments;
};



class Node_Parallel_Segments
{
public:
    Node_Parallel_Segments(double _angle);

    Node_Parallel_Segments(const Node_Parallel_Segments & node);

    ~Node_Parallel_Segments();

    void delete_references_to_parallel_node();

    void delete_colinear_nodes();

    void create_colinear_node(double _ordinate);

    void add(Segment *s);

    void remove(Segment *s);

    void assign_to_colinear_node(double _ordinate, Segment* s);

    void assign_to_colinear_node(double _ordinate, list<Segment *> & ls);

    void assign_to_other(Segment* s);

    void assign_to_other(list<Segment *> & ls);

public:
    double angle;
    Point2d frame_origin;
    map<double, Node_Colinear_Segments*> colinear_segments;
    list<Segment *> other_segments;

    list<Segment *> parallel_segments;
};

#endif // SEGMENT_RAY_H
