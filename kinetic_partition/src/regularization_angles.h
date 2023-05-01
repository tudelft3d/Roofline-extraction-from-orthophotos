#ifndef REGULARIZATION_ANGLES_H
#define REGULARIZATION_ANGLES_H

#include "defs.h"
#include "kinetic_model.h"
#include "segment_ray.h"
#include "segment_tree.h"


class Regularization_Angles
{
public:
    Regularization_Angles();

    virtual ~Regularization_Angles();

    virtual void regularize(Kinetic_Model* model) = 0;


    void build_rtree(vector<Segment *> & segments, Boost_RTree & rtree_segments, bool include_artificial);

    void search_neighborhood(Segment* s, double D, Boost_RTree & rtree_segments, vector<Segment *> & segments, list<Segment *> & neighborhood, bool include_artificial);

    void rotate(Segment_Regularization_Tree* tree);

#if NOT_MEASURING_PERFORMANCES
    void make_layer(Kinetic_Model* model);
#endif

    double dalpha_max_const(Segment* s, Parameters* params);

    double dalpha_max_offset(Segment* s, Parameters* params);

	double dalpha_max_lsd(Segment* s, Parameters* params);
};

#endif // REGULARIZATION_ANGLES_H
