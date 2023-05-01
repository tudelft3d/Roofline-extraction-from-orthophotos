#include "regularization_angles.h"
#include "geometry.h"
#include <cmath>


Regularization_Angles::Regularization_Angles()
{

}


Regularization_Angles::~Regularization_Angles()
{

}


void Regularization_Angles::build_rtree(vector<Segment *> & segments, Boost_RTree & rtree_segments, bool include_artificial)
{
    double x_min, x_max, y_min, y_max;
    for (unsigned int i = 0 ; i < segments.size() ; i++) {

        // For each segment, we create a box that represents the bounding box of this segment
        Segment* s_i = segments[i];
        if (!s_i->is_artificial || (s_i->is_artificial && include_artificial)) {

            // Gets the bounding box's coordinates
            if (s_i->end1.x < s_i->end2.x) {
                x_min = s_i->end1.x; x_max = s_i->end2.x;
            } else {
                x_min = s_i->end2.x; x_max = s_i->end1.x;
            }
            if (s_i->end1.y < s_i->end2.y) {
                y_min = s_i->end1.y; y_max = s_i->end2.y;
            } else {
                y_min = s_i->end2.y; y_max = s_i->end1.y;
            }

            // We insert this box in the r-tree
            Boost_Box s_box(Boost_Point(x_min, y_min), Boost_Point(x_max, y_max));
            rtree_segments.insert(std::make_pair(s_box, i));
        }
    }
}


void Regularization_Angles::search_neighborhood(Segment* s, double D, Boost_RTree & rtree_segments, vector<Segment *> & segments, list<Segment *> & neighborhood, bool include_artificial)
{
    neighborhood.clear();

    // First, we search for segments whose bounding boxes' edges are located at a distance lower than D from the bounding box of s
    double x_min, x_max, y_min, y_max;
    vector<Boost_Value> possible_neighbors;

    if (s->end1.x < s->end2.x) {
        x_min = s->end1.x; x_max = s->end2.x;
    } else {
        x_min = s->end2.x; x_max = s->end1.x;
    }
    if (s->end1.y < s->end2.y) {
        y_min = s->end1.y; y_max = s->end2.y;
    } else {
        y_min = s->end2.y; y_max = s->end1.y;
    }
    Boost_Box query(Boost_Point(x_min - D, y_min - D), Boost_Point(x_max + D, y_max + D));
    rtree_segments.query(bgi::intersects(query), std::back_inserter(possible_neighbors));

    // Second, we refine the list of possible neighbors, by selecting the ones located at a distance lower than D from s
    for (unsigned int i = 0 ; i < possible_neighbors.size() ; i++) {
        Segment* s_i = segments[possible_neighbors[i].second];
        if (!s_i->is_artificial || (s_i->is_artificial && include_artificial)) {
            if (s_i == s || Geometry::distance_initial_coordinates(s, s_i) < D) {
                neighborhood.push_back(s_i);
            }
        }
    }
}


void Regularization_Angles::rotate(Segment_Regularization_Tree *tree)
{
    for (map<double, Node_Parallel_Segments*>::iterator it_m = tree->parallel_segments.begin(); it_m != tree->parallel_segments.end(); it_m++) {
        double theta = it_m->first;
        Node_Parallel_Segments* subtree = it_m->second;

        // Each group of parallel segments has a normal vector that we compute with alpha
        Vec2d v_dir = Vec2d(cos(theta * PI / 180), sin(theta * PI / 180));
        Vec2d v_nor = Vec2d(-v_dir[1], v_dir[0]);
        double a = v_nor[0];
        double b = v_nor[1];

        // Rotates the segments with precision
        for (list<Segment *>::iterator it_s = subtree->parallel_segments.begin(); it_s != subtree->parallel_segments.end(); it_s++) {
            Segment* s = (*it_s);

            // Computes the equation of the supporting line of the rotated segment
            double c = -a * s->barycenter.x - b * s->barycenter.y;
            s->set_dalpha(theta - s->alpha, theta, a, b, c, v_dir);
        }
    }
}


#if NOT_MEASURING_PERFORMANCES
void Regularization_Angles::make_layer(Kinetic_Model *model)
{
	model->clear_line_items(model->L_rega);

    // Initializes a uniform distribution
    std::default_random_engine & generator = model->generator;
    std::uniform_int_distribution<int> uniform_dist(50, 255);

    // Copies the background image
    Matrix<uchar> & I = model->I;

    // Prints segments by cluster, with a random color per cluster
    for (map<double, Node_Parallel_Segments*>::iterator it_m = model->tree->parallel_segments.begin(); it_m != model->tree->parallel_segments.end(); it_m++) {
        uchar r = uchar(uniform_dist(generator));
        uchar g = uchar(uniform_dist(generator));
        uchar b = uchar(uniform_dist(generator));

        if (it_m->second->parallel_segments.size() > 0) {
            // If the second phase of the regularization has not been applied yet
            for (list<Segment *>::iterator it_s = it_m->second->parallel_segments.begin(); it_s != it_m->second->parallel_segments.end(); it_s++) {
                Segment* s = (*it_s);
                Point2d p1 = Point2d(jclamp(0, s->interEnd1.x, I.cols - 1), jclamp(0, I.rows - s->interEnd1.y, I.rows - 1));
                Point2d p2 = Point2d(jclamp(0, s->interEnd2.x, I.cols - 1), jclamp(0, I.rows - s->interEnd2.y, I.rows - 1));
                model->add_line(model->L_rega, p1.x, p1.y, p2.x, p2.y, r, g, b);
            }

        } else {
            // If the second phase of the regularization has already been applied
            for (list<Segment *>::iterator it_s = it_m->second->other_segments.begin() ; it_s != it_m->second->other_segments.end() ; it_s++) {
                Segment* s = (*it_s);
                Point2d p1 = Point2d(jclamp(0, s->interEnd1.x, I.cols - 1), jclamp(0, I.rows - s->interEnd1.y, I.rows - 1));
                Point2d p2 = Point2d(jclamp(0, s->interEnd2.x, I.cols - 1), jclamp(0, I.rows - s->interEnd2.y, I.rows - 1));
                model->add_line(model->L_rega, p1.x, p1.y, p2.x, p2.y, r, g, b);
            }
            for (map<double, Node_Colinear_Segments *>::iterator it_m2 = it_m->second->colinear_segments.begin() ; it_m2 != it_m->second->colinear_segments.end() ; it_m2++) {
                for (list<Segment *>::iterator it_s = it_m2->second->colinear_segments.begin() ; it_s != it_m2->second->colinear_segments.end() ; it_m2++) {
                    Segment* s = (*it_s);
                    Point2d p1 = Point2d(jclamp(0, s->interEnd1.x, I.cols - 1), jclamp(0, I.rows - s->interEnd1.y, I.rows - 1));
                    Point2d p2 = Point2d(jclamp(0, s->interEnd2.x, I.cols - 1), jclamp(0, I.rows - s->interEnd2.y, I.rows - 1));
                    model->add_line(model->L_rega, p1.x, p1.y, p2.x, p2.y, r, g, b);
                }
            }
        }
    }

    // Prints other segments in black
    for (list<Segment *>::iterator it_s = model->tree->other_segments.begin(); it_s != model->tree->other_segments.end(); it_s++) {
        Segment* s = (*it_s);
        Point2d p1 = Point2d(jclamp(0, s->interEnd1.x, I.cols - 1), jclamp(0, I.rows - s->interEnd1.y, I.rows - 1));
        Point2d p2 = Point2d(jclamp(0, s->interEnd2.x, I.cols - 1), jclamp(0, I.rows - s->interEnd2.y, I.rows - 1));
        model->add_line(model->L_rega, p1.x, p1.y, p2.x, p2.y, 0, 0, 0);
    }
}
#endif

double Regularization_Angles::dalpha_max_const(Segment *s, Parameters *params)
{
    return params->rega_angle_const;
}


double Regularization_Angles::dalpha_max_offset(Segment *s, Parameters *params)
{
    return 180 * atan(2 * params->rega_angle_offset / s->length) / PI;
}


double Regularization_Angles::dalpha_max_lsd(Segment *s, Parameters *params)
{
	return 1; // 180 * atan2(0.5 * s->width, 0.5 * s->length) / PI;
}