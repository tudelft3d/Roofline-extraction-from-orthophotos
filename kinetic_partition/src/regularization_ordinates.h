#pragma once
#include "defs.h"
#include "kinetic_model.h"


class Regularization_Ordinates
{
public:
	Regularization_Ordinates();

	virtual ~Regularization_Ordinates();

	virtual void regularize(Kinetic_Model* model) = 0;
	
#if NOT_MEASURING_PERFORMANCES
	void make_layer(Kinetic_Model* model);
#endif

protected:
	void transform_coordinates(double theta, Node_Parallel_Segments* cluster, double & x_min, double & x_max, double & y_min, double & y_max);

	void build_rtree(Node_Parallel_Segments* cluster, Boost_RTree & rtree_segments);

	void search_neighborhood(Segment* s, double dx, double dy, Boost_RTree & rtree_segments, vector<Segment *> & segments, list<Segment *> & neighborhood);

	void translate(Node_Parallel_Segments* node);

public:
	double dt_max_const(Segment* s, Parameters* params);

	double dt_max_lsd_width(Segment* s, Parameters* params);
};

