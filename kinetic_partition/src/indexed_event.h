#pragma once
#include "segment_ray.h"

class IndexedEvent
{
public:
	IndexedEvent(int _intersectant, int _intersected, double _t_intersectant, double _t_intersected, bool _is_colliding_colinear);

	IndexedEvent(int _intersectant, Image_Boundary _intersected, double _t_intersectant, double _t_intersected);

	IndexedEvent(int _intersectant, Image_Boundary _intersected, int outer_edge_index, double _t_intersectant, double _t_intersected);

	IndexedEvent(const IndexedEvent & E);
	
	~IndexedEvent();

    static void clear_schedule(IndexedEvent *&schedule);
	
public:
	bool is_colliding_ray;
	bool is_colliding_colinear;
	int intersectant;
	int intersected;
	int intersected_index;
	double t_intersectant;
	double t_intersected;
	
	IndexedEvent* prev_i;
	IndexedEvent* next_i;
	IndexedEvent* prev_j;
	IndexedEvent* next_j;

	IndexedEvent* previous;
	IndexedEvent* next;
};


inline bool before_indexed_event(IndexedEvent* i, IndexedEvent* j) {
	if (i->t_intersectant < j->t_intersectant) {
		return true;
	}
	else if (i->t_intersectant == j->t_intersectant) {
		return (i->t_intersected < j->t_intersected);
	}
	else {
		return false;
	}

	//return  (i->t_intersectant < j->t_intersectant);
}


inline bool sorted_by_index_of_intersected_segment(IndexedEvent* e_i, IndexedEvent* e_j)
{
    return ((e_i->intersected << 1) < (e_j->intersected << 1));
}
