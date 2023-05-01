#pragma once
#include "defs.h"
#include "kinetic_model.h"
#include "regularization_ordinates.h"
#include "r_ordinate.h"

#include <vector>
#include <random>

using std::vector;
using std::default_random_engine;
using std::uniform_int_distribution;



class Regularization_Ordinates_Mean_Shift : public Regularization_Ordinates
{
public:
	Regularization_Ordinates_Mean_Shift();

	~Regularization_Ordinates_Mean_Shift();

    void regularize(Kinetic_Model* model);

private:
    void build_rtree(vector<R_Ordinate *> & regularities, Boost_RTree & rtree_regularities);

    void search_neighborhood(Segment* s, double dx, double dy, Boost_RTree & rtree_regularities, vector<R_Ordinate *> & regularities, list<R_Ordinate *> & neighborhood);

    void search_neighborhood(R_Ordinate* r, double dx, double dy, Boost_RTree & rtree_regularities, vector<R_Ordinate *> & regularities, list<R_Ordinate *> & neighborhood);


    void identify_regularities(Kinetic_Model* model, Boost_RTree & rtree_segments, Node_Parallel_Segments* cluster, vector<R_Ordinate *> & regularities);

    bool propose_regularity_by_mean_shift(Kinetic_Model* model, Boost_RTree & rtree_segments, Segment* seed, double & location, double & ordinate_observed);

    void smooth_regularities(Kinetic_Model* model, Boost_RTree & rtree_regularities, vector<R_Ordinate *> & all_regularities, Node_Parallel_Segments* cluster);

    bool smooth_by_mean_shift(Kinetic_Model* model, R_Ordinate* seed, Boost_RTree & rtree_regularities, vector<R_Ordinate *> &all_regularities, int min_terms, double & ordinate_smoothed);

    void assign_segments_to_local_regularities(Kinetic_Model* model, Boost_RTree & rtree_regularities, vector<R_Ordinate *> & all_regularities, Node_Parallel_Segments* cluster);

    double kernel(double k, double x) {
        return exp(-k * x * x);
    }
};
