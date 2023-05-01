#include "regularization_angles_ms.h"
#include "geometry.h"
#include "defs.h"
#include "trace.h"


Regularization_Angles_Mean_Shift::Regularization_Angles_Mean_Shift()
    : Regularization_Angles()
{

}



Regularization_Angles_Mean_Shift::~Regularization_Angles_Mean_Shift()
{

}


void Regularization_Angles_Mean_Shift::regularize(Kinetic_Model* model)
{
	clock_t t_begin = clock();

    bool include_artificial = false;

    // Cleans the regularization tree, if another execution has been perform before
    model->tree->delete_parallel_nodes();

    // Builds a R-tree from the list of segments
    Boost_RTree rtree_segments, rtree_regularities;
    Regularization_Angles::build_rtree(model->segments, rtree_segments, include_artificial);

    // We try to identify local regularities on the image
    // A local regularity is a 3D point : x- and y- coordinates represent the center of the neighborhood,
    // and the z-coordinate represent an angle formed by a sufficiently high number of segments in this neighborhood
    vector<R_Angle *> regularities;
    identify_regularities(model, rtree_segments, regularities, include_artificial);
    rtree_segments.clear();

    // We have obtained a list of local regularities. We are now trying to identify global regularities at the scale of the image
    // To this end, again, we apply the mean-shift algorithm on the set of local regularities we identified before.
    // This will allow us to regroup points on planes, that will stand for our global regularities.
    build_rtree(regularities, rtree_regularities);
    smooth_regularities(model, rtree_regularities, regularities);

    // Finally, we assign segments to planes of regularity
    assign_segments_to_local_regularities(model, rtree_regularities, regularities);

    // Applies transformations
    rotate(model->tree);

    // Cleans memeory
    for (int r = 0; r < int(regularities.size()); r++) delete regularities[r];
    rtree_regularities.clear();

    model->applied_regularization_angles = true;

#if NOT_MEASURING_PERFORMANCES
    make_layer(model);
#endif

	clock_t t_end = clock();
	trace(model->params->verbose_level, 5, "** Regularization (part 1) done in " + std::to_string(float(t_end - t_begin) / CLOCKS_PER_SEC) + " s.");
}



void Regularization_Angles_Mean_Shift::build_rtree(vector<R_Angle *> & regularities, Boost_RTree & rtree_regularities)
{
    double x_min, x_max, y_min, y_max;
    for (unsigned int i = 0 ; i < regularities.size() ; i++) {

        // For each segment, we create a box that represents the bounding box of this segment
        R_Angle* r_i = regularities[i];
        Segment* s_i = r_i->segment;

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
        Boost_Box r_box(Boost_Point(x_min, y_min), Boost_Point(x_max, y_max));
        rtree_regularities.insert(std::make_pair(r_box, i));
    }
}



void Regularization_Angles_Mean_Shift::search_neighborhood(Segment* s, double D, Boost_RTree & rtree_regularities, vector<R_Angle *> & regularities, list<R_Angle *> & neighborhood)
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
    rtree_regularities.query(bgi::intersects(query), std::back_inserter(possible_neighbors));

    // Second, we refine the list of possible neighbors, by selecting the ones located at a distance lower than D from s
    for (unsigned int i = 0 ; i < possible_neighbors.size() ; i++) {
        R_Angle* r_i = regularities[possible_neighbors[i].second];
        if (r_i->segment == s || Geometry::distance_initial_coordinates(s, r_i->segment) < D) {
            neighborhood.push_back(r_i);
        }
    }
}


void Regularization_Angles_Mean_Shift::search_neighborhood(R_Angle* r, double D, Boost_RTree & rtree_regularities, vector<R_Angle *> & regularities, list<R_Angle *> & neighborhood)
{
    search_neighborhood(r->segment, D, rtree_regularities, regularities, neighborhood);
}



void Regularization_Angles_Mean_Shift::identify_regularities(Kinetic_Model *model, Boost_RTree & rtree_segments, vector<R_Angle *> & regularities, bool include_artificial)
{
    uint index = -1;
    vector<Segment *> & segments = model->segments;
    regularities.reserve(segments.size());

    for (vector<Segment *>::iterator it_s = segments.begin(); it_s != segments.end(); it_s++) {
        Segment* s = (*it_s);
        if (!s->is_artificial || (s->is_artificial && include_artificial)) {
            double value;
            if (propose_regularity_by_mean_shift(model, rtree_segments, s, value, include_artificial)) {
                regularities.push_back(new R_Angle(++index, s, value, 1));
            }
        }
    }
    regularities.shrink_to_fit();
}



bool Regularization_Angles_Mean_Shift::propose_regularity_by_mean_shift(Kinetic_Model* model, Boost_RTree & rtree_segments, Segment* seed, double & value, bool include_artificial)
{
    vector<Segment *> & segments = model->segments;

    // Parameters of the mean shift algorithm
    double eps = model->params->rega_ms_epsilon;
    double sigma = model->params->rega_ms_sigma;
    double D = model->params->rega_ms_distance;
    double gamma = 1.0 / (2 * sigma * sigma);
    int min_terms = model->params->rega_ms_min_terms;

    list<Segment *> neighborhood;
    Regularization_Angles::search_neighborhood(seed, D, rtree_segments, segments, neighborhood, include_artificial);

    // We apply the mean shift algorithm to find the best alpha value for the seed segment
    // and the ones that are, supposedly, almost parallel to it
    double x = seed->alpha;

    // While convergence is not reached
    while (true) {

        // We compute a new value of x
        double num_m = 0;
        double den_m = 0;
        int terms = 0;
        for (list<Segment *>::iterator it_s = neighborhood.begin(); it_s != neighborhood.end(); it_s++) {
            Segment* s = (*it_s);
            if (inside_interval_mod_180(x, s->alpha, sigma)) {
                double k = kernel(gamma, s->alpha - x);
                num_m += k * s->length * s->alpha;
                den_m += k * s->length;
                terms++;
            }
        }

        if (terms < min_terms) {
            // Eliminates propositions that are formulated with an insufficient number of inliers
            return false;
        } else {
            double m = num_m / den_m;
            if (fabs(m - x) < eps) {
                value = m;
                return true;
            } else {
                x = m;
            }
        }
    }
}



void Regularization_Angles_Mean_Shift::smooth_regularities(Kinetic_Model *model, Boost_RTree &rtree_regularities, vector<R_Angle *> &all_regularities)
{
    int t = 0;
    int r_size = int(all_regularities.size());
    int tries = 10 * r_size;
    std::uniform_int_distribution<int> uniform_dist;
    Segment_Regularization_Tree* tree = model->tree;

    vector<R_Angle *> non_smoothed_regularities = vector<R_Angle *>(all_regularities.begin(), all_regularities.end());
    while (non_smoothed_regularities.size() > 0)
    {
        double z;
        int min_terms = 1 + MAX(0, (tries - t) / r_size);
		int seed_index = 0; // uniform_dist(model->generator) % int(non_smoothed_regularities.size());
        if (smooth_by_mean_shift(model, non_smoothed_regularities[seed_index], rtree_regularities, all_regularities, min_terms, z))
        {
            // If we have obtained a smoothed value via the mean-shift algorithm
            // Then we create an entry in the map of regularized segments for this value
            tree->create_parallel_node(z);

            // Converts the vector of regularities into a vector, erasing elements should be faster
            list<R_Angle *> l_non_smoothed_regularities = list<R_Angle *>(non_smoothed_regularities.begin(), non_smoothed_regularities.end());
            list<R_Angle *>::iterator it_r = l_non_smoothed_regularities.begin();
            while (it_r != l_non_smoothed_regularities.end()) {
                R_Angle* r = (*it_r);
                if (fabs(r->angle_observed - z) <= model->params->rega_ms_smooth_dist) {

                    // Removes the box corresponding to r from the R-Tree
                    Segment* s = r->segment;
                    double x_min, x_max, y_min, y_max;
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
                    Boost_Box r_box(Boost_Point(x_min, y_min), Boost_Point(x_max, y_max));
                    rtree_regularities.remove(std::make_pair(r_box, r->index));

                    // Removes the regularity from the list
                    it_r = l_non_smoothed_regularities.erase(it_r);

                    // Smoothes the regularity with the provided value
                    r->smooth(z);

                } else {
                    it_r++;
                }
            }
            non_smoothed_regularities = vector<R_Angle *>(l_non_smoothed_regularities.begin(), l_non_smoothed_regularities.end());
            l_non_smoothed_regularities.clear();
        }
        t++;
    }

    // The R-Tree has been progressively emptied : we need to rebuild it again
    build_rtree(all_regularities, rtree_regularities);
}



bool Regularization_Angles_Mean_Shift::smooth_by_mean_shift(Kinetic_Model* model, R_Angle* seed, Boost_RTree & rtree_regularities, vector<R_Angle *> &all_regularities, int min_terms, double & z_0)
{
    double eps = model->params->rega_ms_epsilon;
    double sigma = model->params->rega_ms_sigma;
    double D = model->params->rega_ms_distance;
    double gamma = 1.0 / (2 * sigma * sigma);

    list<R_Angle *> neighborhood;
    search_neighborhood(seed, D, rtree_regularities, all_regularities, neighborhood);

    // We apply the mean shift algorithm to find the best alpha value for the seed segment
    // and the ones that are, supposedly, almost parallel to it
    double z = seed->angle_observed;
    bool first_iteration = true;

    // While convergence is not reached
    while (true) {

        // We compute a new value of x
        double num_m = 0;
        double den_m = 0;
        int terms = 0;
        for (list<R_Angle *>::iterator it_r = neighborhood.begin(); it_r != neighborhood.end(); it_r++) {
            R_Angle* r = (*it_r);
            if (inside_interval_mod_180(z, r->angle_observed, sigma)) {
                double k = kernel(gamma, r->angle_observed - z);
                num_m += k * r->segment->length * r->angle_observed;
                den_m += k * r->segment->length;
                terms++;
            }
        }

        if (first_iteration) {
            if (terms < min_terms) return false;
            first_iteration = false;
        }
        double m = num_m / den_m;
        if (fabs(m - z) < eps) {
            z_0 = m;
            return true;
        } else {
            z = m;
        }
    }
}



void Regularization_Angles_Mean_Shift::assign_segments_to_local_regularities(Kinetic_Model* model, Boost_RTree & rtree_regularities, vector<R_Angle *> & all_regularities)
{
    double (Regularization_Angles::*dalpha_function)(Segment*, Parameters*) = NULL;
	switch (model->params->rega_angle_function) {
	case 0: dalpha_function = &Regularization_Angles::dalpha_max_const; break;
	case 1: dalpha_function = &Regularization_Angles::dalpha_max_offset; break;
	case 2: dalpha_function = &Regularization_Angles::dalpha_max_lsd; break;
	}

    vector<Segment *> & segments = model->segments;
    Segment_Regularization_Tree* tree = model->tree;
    double D = model->params->rega_ms_distance;

    // We try to assign each segment to a local regularity
    for (vector<Segment *>::iterator it_s = segments.begin(); it_s != segments.end(); it_s++) {
        Segment* s = (*it_s);

        // Searches for such regularities
        list<R_Angle *> local_regularities;
        search_neighborhood(s, D, rtree_regularities, all_regularities, local_regularities);

        // If the segment is a kind of outlier, skip it
		// if (local_regularities.size() == 0) continue;

        R_Angle* best_regularity = NULL;
        double d_min = D;
        for (list<R_Angle *>::iterator it_r = local_regularities.begin(); it_r != local_regularities.end(); it_r++) {
            R_Angle* r = (*it_r);
            double d = Geometry::distance_initial_coordinates(r->segment, s);

            // If the current segment is sufficiently close to a detected regularity,
            // and if its value is close enough to the value of the local regularity
            if ((d < d_min) && fabs_mod_180(r->angle_smoothed, s->alpha, (this->*dalpha_function)(s, model->params))) {
                d_min = d;
                best_regularity = r;
            }
        }

        if (best_regularity != NULL) {
            tree->assign_to_parallel_node(best_regularity->angle_smoothed, s);
        } else {
            tree->assign_to_other(s);
        }
    }

	// Forces segments to be assigned to a node
	double alpha_eps = 0.5;
	list<Segment *>::iterator it_s = tree->other_segments.begin(); 
	while (it_s != tree->other_segments.end()) {
		double alpha = (*it_s)->alpha;

		Node_Parallel_Segments* node = nullptr;
		for (map<double, Node_Parallel_Segments*>::iterator it_m = tree->parallel_segments.begin() ; it_m != tree->parallel_segments.end() ; it_m++) {
			double theta = it_m->first;
			for (int k = -1 ; k <= 1 ; k++) {
				if (fabs(theta - alpha + 180 * k) < alpha_eps) {
					node = it_m->second;
					break;
				}
			}
			if (node != nullptr) break;
		}

		if (node != nullptr) {
			node->parallel_segments.push_back(*it_s);
		} else {
			tree->create_parallel_node(alpha);
			tree->assign_to_parallel_node(alpha, *it_s);
		}
		it_s = tree->other_segments.erase(it_s);
	}
	
}
