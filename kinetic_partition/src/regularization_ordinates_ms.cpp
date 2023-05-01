#include "regularization_ordinates_ms.h"
#include "geometry.h"
#include "trace.h"

#define BBOX_EPS 0.001


Regularization_Ordinates_Mean_Shift::Regularization_Ordinates_Mean_Shift()
	: Regularization_Ordinates()
{ 

}


Regularization_Ordinates_Mean_Shift::~Regularization_Ordinates_Mean_Shift()
{

}



void Regularization_Ordinates_Mean_Shift::regularize(Kinetic_Model* model)
{
	clock_t t_begin = clock();

    Segment_Regularization_Tree* tree = model->tree;

    map<double, Node_Parallel_Segments*>::iterator it_c = tree->parallel_segments.begin();
    while (it_c != tree->parallel_segments.end()) {

        // Accesses the cluster
        double theta = it_c->first;
        Node_Parallel_Segments* cluster = it_c->second;
        cluster->delete_colinear_nodes();

        if (cluster->parallel_segments.size() == 0) {
            it_c = tree->parallel_segments.erase(it_c);
            continue;
        }

        // Transforms the referencing coordinates of the segments included in the current cluster
        // From now on they will be expressed in a local frame
        vector<R_Ordinate *> regularities;
        double x_min, x_max, y_min, y_max;
        transform_coordinates(theta, cluster, x_min, x_max, y_min, y_max);

        Boost_RTree rtree_segments;
        Regularization_Ordinates::build_rtree(cluster, rtree_segments);

        // We would like to encourage colinearity relationships between segments
        // Segments are colinear when they have the same y-coordinate in the local frame
        identify_regularities(model, rtree_segments, cluster, regularities);
        rtree_segments.clear();

        // We have obtained a list of local regularities. We are now trying to identify global regularities at the scale of the image
        // To this end, again, we apply the mean-shift algorithm on the set of local regularities we identified before.
        // This will allow us to regroup points on lines, that will stand for our global regularities.
        Boost_RTree rtree_regularities;
        build_rtree(regularities, rtree_regularities);
        smooth_regularities(model, rtree_regularities, regularities, cluster);

        // Finally, we assign segments to lines of regularity
        assign_segments_to_local_regularities(model, rtree_regularities, regularities, cluster);

        // Applies transformations
        translate(cluster);

        // Cleans memory
        for (int r = 0; r < int(regularities.size()); r++) delete regularities[r];
        rtree_regularities.clear();

        it_c++;
    }

    model->applied_regularization_ordinates = true;

	clock_t t_end = clock();
	trace(model->params->verbose_level, 5, "** Regularization (part 2) done in : " + std::to_string(float(t_end - t_begin) / CLOCKS_PER_SEC) + " s.");
}






void Regularization_Ordinates_Mean_Shift::build_rtree(vector<R_Ordinate *> & regularities, Boost_RTree & rtree_regularities)
{
    for (int i = 0 ; i < regularities.size() ; i++) {
        R_Ordinate* r = regularities[i];

        double & rx = r->segment->referencing_coordinates.x;
        double & ry = r->y_observed;
        double & length = r->segment->length;

        // Creation of a box representing the regularity
        Boost_Box r_box(Boost_Point(rx - length / 2, ry - BBOX_EPS), Boost_Point(rx + length / 2, ry + BBOX_EPS));
        rtree_regularities.insert(std::make_pair(r_box, r->index));
    }
}






void Regularization_Ordinates_Mean_Shift::search_neighborhood(Segment* s, double dx, double dy, Boost_RTree & rtree_regularities, vector<R_Ordinate *> & regularities, list<R_Ordinate *> & neighborhood)
{
    neighborhood.clear();

    vector<Boost_Value> neighbors;
    Point2d ref_coords = s->referencing_coordinates;
    Boost_Box query(Boost_Point(ref_coords.x - (s->length + dx) / 2, ref_coords.y - dy / 2),
                    Boost_Point(ref_coords.x + (s->length + dx) / 2, ref_coords.y + dy / 2));
    rtree_regularities.query(bgi::intersects(query), std::back_inserter(neighbors));

    for (unsigned int i = 0 ; i < neighbors.size() ; i++) {
        R_Ordinate* r_i = regularities[neighbors[i].second];
        neighborhood.push_back(r_i);
    }
}



void Regularization_Ordinates_Mean_Shift::search_neighborhood(R_Ordinate* r, double dx, double dy, Boost_RTree & rtree_regularities, vector<R_Ordinate *> & regularities, list<R_Ordinate *> & neighborhood)
{
    neighborhood.clear();

    vector<Boost_Value> neighbors;
    double & rx = r->segment->referencing_coordinates.x;
    double & ry = r->y_observed;
    double & length = r->segment->length;

    Boost_Box query (Boost_Point(rx - (length + dx) / 2, ry - dy / 2), Boost_Point(rx + (length + dx) / 2, ry + dy / 2));
    rtree_regularities.query(bgi::intersects(query), std::back_inserter(neighbors));

    for (unsigned int i = 0 ; i < neighbors.size() ; i++) {
        R_Ordinate* r_i = regularities[neighbors[i].second];
        neighborhood.push_back(r_i);
    }
}



void Regularization_Ordinates_Mean_Shift::identify_regularities(Kinetic_Model* model, Boost_RTree & rtree_segments, Node_Parallel_Segments* cluster, vector<R_Ordinate *> & regularities)
{
    uint index = -1;
    regularities.reserve(cluster->parallel_segments.size());
    for (list<Segment *>::iterator it_s = cluster->parallel_segments.begin(); it_s != cluster->parallel_segments.end(); it_s++) {
        double location, ordinate_observed;
        bool process_completed = propose_regularity_by_mean_shift(model, rtree_segments, (*it_s), location, ordinate_observed);
        if (process_completed) {
            regularities.push_back(new R_Ordinate(++index, (*it_s), ordinate_observed, 1));
        }
    }
    regularities.shrink_to_fit();
}



bool Regularization_Ordinates_Mean_Shift::propose_regularity_by_mean_shift(Kinetic_Model *model, Boost_RTree &rtree_segments, Segment *seed, double &location, double &ordinate_observed)
{
    // Parameters of the mean shift algorithm
    double eps = model->params->regp_ms_epsilon;
    double sigma = model->params->regp_ms_sigma;
    double dx = model->params->regp_ms_distx, dy = model->params->regp_ms_disty;

    double gamma = 1.0 / (2 * sigma * sigma);
    int min_terms = model->params->regp_ms_min_terms;

    list<Segment *> neighborhood;
    Regularization_Ordinates::search_neighborhood(seed, dx, dy, rtree_segments, model->segments, neighborhood);

    // We apply the mean shift algorithm to find the best alpha value for the seed segment
    // and the ones that are, supposedly, almost parallel to it
    double y = seed->referencing_coordinates.y;

    // While convergence is not reached
    while (true) {

        // We compute a new value of x
        double num_m = 0;
        double den_m = 0;
        int terms = 0;
        for (list<Segment *>::iterator it_s = neighborhood.begin(); it_s != neighborhood.end(); it_s++) {
            Segment* s = (*it_s);
            if (fabs(y - s->referencing_coordinates.y) <= sigma) {
                double k = kernel(gamma, s->referencing_coordinates.y - y);
                num_m += k * s->referencing_coordinates.y;
                den_m += k;
                ++terms;
            }
        }

        if (terms < min_terms) {
            // Eliminates propositions that are formulated with an insufficient number of inliers
            return false;
        } else {
            double m = num_m / den_m;
            if (fabs(m - y) < eps) {
                location = seed->referencing_coordinates.x;
                ordinate_observed = m;
                return true;
            } else {
                y = m;
            }
        }
    }
}



void Regularization_Ordinates_Mean_Shift::smooth_regularities(Kinetic_Model *model, Boost_RTree &rtree_regularities, vector<R_Ordinate *> & all_regularities, Node_Parallel_Segments *cluster)
{
    int t = 0;
    int r_size = int(all_regularities.size());
    int tries = 10 * r_size;
    std::uniform_int_distribution<int> uniform_dist;

    vector<R_Ordinate *> non_smoothed_regularities = vector<R_Ordinate *>(all_regularities.begin(), all_regularities.end());
    while (non_smoothed_regularities.size() > 0)
    {
        double y;
        int min_terms = 1 + MAX(0, (tries - t) / r_size);
		int seed_index = 0; // uniform_dist(model->generator) % int(non_smoothed_regularities.size());
        if (smooth_by_mean_shift(model, non_smoothed_regularities[seed_index], rtree_regularities, all_regularities, min_terms, y)) {

            // If we could obtain a smoothed value via the mean-shift algorithm
            // Then we create an entry in the map of regularized segments for this value
            cluster->create_colinear_node(y);

            // Converts the vector of regularities into a vector, erasing elements should be faster
            list<R_Ordinate *> l_non_smoothed_regularities = list<R_Ordinate *>(non_smoothed_regularities.begin(), non_smoothed_regularities.end());
            list<R_Ordinate *>::iterator it_r = l_non_smoothed_regularities.begin();
            while (it_r != l_non_smoothed_regularities.end()) {
                R_Ordinate* r = (*it_r);
                if (fabs(r->y_observed - y) <= model->params->regp_ms_smooth_dist) {

                    // Removes the box corresponding to r from the R-Tree
                    double & rx = r->segment->referencing_coordinates.x;
                    double & ry = r->y_observed;
                    double & length = r->segment->length;
                    Boost_Box r_box(Boost_Point(rx - length / 2, ry - BBOX_EPS), Boost_Point(rx + length / 2, ry + BBOX_EPS));
                    rtree_regularities.remove(std::make_pair(r_box, r->index));

                    // Removes the regularity from the list
                    it_r = l_non_smoothed_regularities.erase(it_r);

                    // Smoothes the regularity with the provided value
                    r->smooth(y);

                } else {
                    it_r++;
                }
            }

            non_smoothed_regularities = vector<R_Ordinate *>(l_non_smoothed_regularities.begin(), l_non_smoothed_regularities.end());
            l_non_smoothed_regularities.clear();
        }
        t++;
    }

    // The R-tree has been progressively emptied : we need to rebuild it again
    build_rtree(all_regularities, rtree_regularities);
}



bool Regularization_Ordinates_Mean_Shift::smooth_by_mean_shift(Kinetic_Model *model, R_Ordinate *seed, Boost_RTree &rtree_regularities, vector<R_Ordinate *> & all_regularities, int min_terms, double &ordinate_smoothed)
{
    double eps = model->params->regp_ms_epsilon;
    double sigma = model->params->regp_ms_sigma;
    double dx = model->params->regp_ms_distx, dy = model->params->regp_ms_disty;
    double gamma = 1.0 / (2 * sigma * sigma);

    list<R_Ordinate *> neighborhood;
    search_neighborhood(seed, dx, dy, rtree_regularities, all_regularities, neighborhood);

    // We apply the mean shift algorithm to find the best alpha value for the seed segment
    // and the ones that are, supposedly, almost parallel to it
    double y = seed->y_observed;
    bool first_iteration = true;

    // While convergence is not reached
    while (true) {

        // We compute a new value of x
        double num_m = 0;
        double den_m = 0;
        int terms = 0;
        for (list<R_Ordinate *>::iterator it_r = neighborhood.begin(); it_r != neighborhood.end(); it_r++) {
            R_Ordinate* r = (*it_r);
            if (y - sigma <= r->y_observed && r->y_observed <= y + sigma) {
                double k = kernel(gamma, r->y_observed - y);
                num_m += k * r->y_observed;
                den_m += k;
                terms++;
            }
        }

        if (first_iteration) {
            if (terms < min_terms) return false;
            first_iteration = false;
        }
        double m = num_m / den_m;
        if (fabs(m - y) < eps) {
            ordinate_smoothed = m;
            return true;
        } else {
            y = m;
        }
    }
}



void Regularization_Ordinates_Mean_Shift::assign_segments_to_local_regularities(Kinetic_Model* model, Boost_RTree &rtree_regularities, vector<R_Ordinate *> & all_regularities, Node_Parallel_Segments* cluster)
{
	double (Regularization_Ordinates_Mean_Shift::*dt_function)(Segment*, Parameters*) = NULL;
    dt_function = &Regularization_Ordinates_Mean_Shift::dt_max_const;

    list<Segment *>::iterator it_s = cluster->parallel_segments.begin();
    while (it_s != cluster->parallel_segments.end()) {

        Segment* s = (*it_s);
        list<R_Ordinate *> local_regularities;

        double dx = model->params->regp_ms_distx;
        double dy = model->params->regp_ms_disty;
		search_neighborhood(s, dx, dy, rtree_regularities, all_regularities, local_regularities);

        if (local_regularities.size() == 0) {
            cluster->assign_to_other(s);
        } else {
            R_Ordinate* best_regularity = NULL;
            double d_min = (this->*dt_function)(s, model->params);

            for (list<R_Ordinate *>::iterator it_r = local_regularities.begin(); it_r != local_regularities.end(); it_r++) {
                R_Ordinate* r = (*it_r);

                // If the segment is sufficiently close to a measured regularity,
                // and if its value is close enough to the value of this regularity
                double d = fabs(r->y_smoothed - s->referencing_coordinates.y);
                if (d < d_min) {
                    d_min = d;
                    best_regularity = r;
                }
            }

            if (best_regularity != NULL) {
                cluster->assign_to_colinear_node(best_regularity->y_smoothed, s);
            } else {
                cluster->assign_to_other(s);
            }
        }

        it_s = cluster->parallel_segments.erase(it_s);
    }

	// Forces segments to be assigned to a node
	double y_eps = 1;
	it_s = cluster->other_segments.begin();
	while (it_s != cluster->other_segments.end()) {
		Segment* s = (*it_s);
		Node_Colinear_Segments* node = nullptr;
		for (map<double, Node_Colinear_Segments*>::iterator it_m = cluster->colinear_segments.begin() ; it_m != cluster->colinear_segments.end() ; it_m++) {
			if (fabs(s->referencing_coordinates.y - it_m->first) < y_eps) {
				node = it_m->second;
				break;
			}
		}
		if (node != nullptr) {
			node->add(s);
		} else {
			cluster->create_colinear_node(s->referencing_coordinates.y);
			cluster->assign_to_colinear_node(s->referencing_coordinates.y, s);
		}
		it_s = cluster->other_segments.erase(it_s);
	}
}
