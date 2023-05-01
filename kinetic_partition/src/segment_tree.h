#pragma once
#include "segment_ray.h"
#include <list>
#include <map>

using std::list;
using std::map;


class Segment_Regularization_Tree
{
public:
    Segment_Regularization_Tree();

    Segment_Regularization_Tree(const Segment_Regularization_Tree & tree);

    ~Segment_Regularization_Tree();

    void delete_parallel_nodes();

	void create_parallel_node(double _angle);

	void assign_to_parallel_node(double _angle, Segment* s);

	void assign_to_parallel_node(double _angle, list<Segment *> & ls);

	void assign_to_other(Segment* s);

	void assign_to_other(list<Segment *> & ls);

    bool empty();

public:
	map<double, Node_Parallel_Segments*> parallel_segments;
	list<Segment *> other_segments;
};
