#include "regularization_ordinates.h"
#include <random>


Regularization_Ordinates::Regularization_Ordinates()
{
}


Regularization_Ordinates::~Regularization_Ordinates()
{
}


void Regularization_Ordinates::transform_coordinates(double theta, Node_Parallel_Segments* cluster, double & x_min, double & x_max, double & y_min, double & y_max)
{
    x_min = FLT_MAX; x_max = -FLT_MAX;
    y_min = FLT_MAX; y_max = -FLT_MAX;
    cluster->frame_origin = cluster->parallel_segments.front()->barycenter;

    // We are going to compute the coordinates of the segments in the frame defined by
    // the center of first segment and the cluster, and the vectors I(cos(theta), sin(theta)) and J(-sin(theta), cos(theta))
    for (list<Segment *>::iterator it_s = cluster->parallel_segments.begin(); it_s != cluster->parallel_segments.end(); it_s++) {

        Segment* s = (*it_s);
        double x = (s->interBarycenter.x - cluster->frame_origin.x) * cos(PI * theta / 180.0) + (s->interBarycenter.y - cluster->frame_origin.y) * sin(PI * theta / 180.0);
        double y = (s->interBarycenter.y - cluster->frame_origin.y) * cos(PI * theta / 180.0) - (s->interBarycenter.x - cluster->frame_origin.x) * sin(PI * theta / 180.0);

        // Updates the extrema
        if (x < x_min) x_min = x;
        if (x > x_max) x_max = x;
        if (y < y_min) y_min = y;
        if (y > y_max) y_max = y;
        s->referencing_coordinates = Point2d(x, y);
    }
}


void Regularization_Ordinates::build_rtree(Node_Parallel_Segments *cluster, Boost_RTree &rtree_segments)
{
    for (list<Segment *>::iterator it_s = cluster->parallel_segments.begin() ; it_s != cluster->parallel_segments.end() ; it_s++) {
        Segment* s = (*it_s);

        Point2d & coords = s->referencing_coordinates;
        double & length = s->length;

        // Creation of a box that represents the segment's bounding box with coordinates expressed in the reoriented frame
        Boost_Box s_box(Boost_Point(coords.x - length / 2, coords.y - 0.001), Boost_Point(coords.x + length / 2, coords.y + 0.001));
        rtree_segments.insert(std::make_pair(s_box, s->index));
    }
}


void Regularization_Ordinates::search_neighborhood(Segment* s, double dx, double dy, Boost_RTree & rtree_segments, vector<Segment *> & segments, list<Segment *> & neighborhood)
{
    neighborhood.clear();

    // Unlike what's done in Regularization_Angles we only search for segments intersecting a box defined by s, dx, dy
    vector<Boost_Value> neighbors;

    Point2d ref_coords = s->referencing_coordinates;
    Boost_Box query(Boost_Point(ref_coords.x - (s->length + dx) / 2, ref_coords.y - dy / 2),
                    Boost_Point(ref_coords.x + (s->length + dx) / 2, ref_coords.y + dy / 2));
    rtree_segments.query(bgi::intersects(query), std::back_inserter(neighbors));

    for (unsigned int i = 0 ; i < neighbors.size() ; i++) {
        Segment* s_i = segments[neighbors[i].second];
        neighborhood.push_back(s_i);
    }
}


void Regularization_Ordinates::translate(Node_Parallel_Segments* node)
{
	map<double, Node_Colinear_Segments*>::iterator it_m = node->colinear_segments.begin();
	while (it_m != node->colinear_segments.end()) {
        double dt = it_m->first;
        Node_Colinear_Segments* subnode = it_m->second;

		if (subnode->colinear_segments.empty()) {
			it_m = node->colinear_segments.erase(it_m);
			continue;
		}

        // Gets the longest segment
        double l_max = -FLT_MAX;
        Segment* s_longest = NULL;
        for (list<Segment *>::iterator it_s = subnode->colinear_segments.begin() ; it_s != subnode->colinear_segments.end() ; it_s++) {
            if ((*it_s)->length > l_max) {
                l_max = (*it_s)->length;
                s_longest = (*it_s);
            }
        }

        // if (s_longest == NULL) continue;

        // Translates the longest segment and gets the line equation
        s_longest->set_dt(dt - s_longest->referencing_coordinates.y);
        double a = s_longest->a, b = s_longest->b, c = s_longest->c;
        Vec2d dir = s_longest->finalDirection;

        // Translates the other segments, so that they rest upon the line ax + by + c = 0
        for (list<Segment *>::iterator it_s = subnode->colinear_segments.begin(); it_s != subnode->colinear_segments.end(); it_s++) {
            if ((*it_s) != s_longest) {
                (*it_s)->set_dt(dt - (*it_s)->referencing_coordinates.y, a, b, c, dir);
            }
        }

		it_m++;
    }
}


#if NOT_MEASURING_PERFORMANCES
void Regularization_Ordinates::make_layer(Kinetic_Model *model)
{
	model->clear_line_items(model->L_regp);

    // Initializes a uniform distribution and copies the background image
    default_random_engine & generator = model->generator;
    std::uniform_int_distribution<int> uniform_dist(50, 255);

    Segment_Regularization_Tree* tree = model->tree;
    Matrix<uchar> & I = model->I;

    // Loops on the list of clusters of parallel segments
    // For each cluster, there exists lists of colinear segments, and a list of segments that are not colinear with anyone else
    // We print each cluster of colinear segments with a random color, others are in black
    for (map<double, Node_Parallel_Segments*>::iterator it_m1 = tree->parallel_segments.begin(); it_m1 != tree->parallel_segments.end(); it_m1++) {
        Node_Parallel_Segments* subtree = it_m1->second;
        for (map<double, Node_Colinear_Segments*>::iterator it_m2 = subtree->colinear_segments.begin(); it_m2 != subtree->colinear_segments.end(); it_m2++) {
            uchar r = uchar(uniform_dist(generator));
            uchar g = uchar(uniform_dist(generator));
            uchar b = uchar(uniform_dist(generator));
            for (list<Segment *>::iterator it_s = it_m2->second->colinear_segments.begin(); it_s != it_m2->second->colinear_segments.end(); it_s++) {
                Segment* s = (*it_s);
                Point2d p1 = Point2d(jclamp(0, s->finalEnd1.x, I.cols - 1), jclamp(0, I.rows - s->finalEnd1.y, I.rows - 1));
                Point2d p2 = Point2d(jclamp(0, s->finalEnd2.x, I.cols - 1), jclamp(0, I.rows - s->finalEnd2.y, I.rows - 1));
                model->add_line(model->L_regp, p1.x, p1.y, p2.x, p2.y, r, g, b);
            }
        }
        for (list<Segment *>::iterator it_s = subtree->other_segments.begin(); it_s != subtree->other_segments.end(); it_s++) {
            Segment* s = (*it_s);
            Point2d p1 = Point2d(jclamp(0, s->finalEnd1.x, I.cols - 1), jclamp(0, I.rows - s->finalEnd1.y, I.rows - 1));
            Point2d p2 = Point2d(jclamp(0, s->finalEnd2.x, I.cols - 1), jclamp(0, I.rows - s->finalEnd2.y, I.rows - 1));
            model->add_line(model->L_regp, p1.x, p1.y, p2.x, p2.y, 0, 0, 0);
        }
    }

    // Prints other segments in black
    for (list<Segment *>::iterator it_s = tree->other_segments.begin(); it_s != tree->other_segments.end(); it_s++) {
        Segment* s = (*it_s);
        Point2d p1 = Point2d(jclamp(0, s->finalEnd1.x, I.cols - 1), jclamp(0, I.rows - s->finalEnd1.y, I.rows - 1));
        Point2d p2 = Point2d(jclamp(0, s->finalEnd2.x, I.cols - 1), jclamp(0, I.rows - s->finalEnd2.y, I.rows - 1));
        model->add_line(model->L_regp, p1.x, p1.y, p2.x, p2.y, 0, 0, 0);
    }
}
#endif

double Regularization_Ordinates::dt_max_const(Segment* s, Parameters* params)
{
	return params->regp_trans_const;
}


double Regularization_Ordinates::dt_max_lsd_width(Segment* s, Parameters* params)
{
	return 2; // s->width / 2;
}