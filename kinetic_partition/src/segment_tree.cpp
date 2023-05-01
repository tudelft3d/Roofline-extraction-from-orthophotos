#include "segment_tree.h"


Segment_Regularization_Tree::Segment_Regularization_Tree()
{
	parallel_segments = map<double, Node_Parallel_Segments*>();
	other_segments = list<Segment *>();
}


Segment_Regularization_Tree::Segment_Regularization_Tree(const Segment_Regularization_Tree & tree)
{
    parallel_segments = map<double, Node_Parallel_Segments*>(tree.parallel_segments.begin(), tree.parallel_segments.end());
    other_segments = list<Segment *>(other_segments.begin(), other_segments.end());
}


Segment_Regularization_Tree::~Segment_Regularization_Tree()
{
    delete_parallel_nodes();
}


void Segment_Regularization_Tree::delete_parallel_nodes()
{
    for (std::map<double, Node_Parallel_Segments*>::iterator it_m = parallel_segments.begin() ; it_m != parallel_segments.end() ; it_m++) {
        delete it_m->second;
    }
    parallel_segments.clear();
    other_segments.clear();
}


void Segment_Regularization_Tree::create_parallel_node(double _angle)
{
	if (parallel_segments.find(_angle) == parallel_segments.end()) {
		parallel_segments[_angle] = new Node_Parallel_Segments(_angle);
	}
}


void Segment_Regularization_Tree::assign_to_parallel_node(double _angle, Segment* s)
{
	if (parallel_segments.find(_angle) != parallel_segments.end()) {
		parallel_segments[_angle]->add(s);
	}
}


void Segment_Regularization_Tree::assign_to_parallel_node(double _angle, list<Segment *> & ls)
{
	if (parallel_segments.find(_angle) != parallel_segments.end()) {
		for (list<Segment *>::iterator it_s = ls.begin(); it_s != ls.end(); it_s++) {
			parallel_segments[_angle]->add(*it_s);
		}
	}
}


void Segment_Regularization_Tree::assign_to_other(Segment* s)
{
	other_segments.push_back(s);
}


void Segment_Regularization_Tree::assign_to_other(list<Segment *> & ls)
{
	for (list<Segment *>::iterator it_s = ls.begin(); it_s != ls.end(); it_s++) {
		other_segments.push_back((*it_s));
	}
}


bool Segment_Regularization_Tree::empty()
{
    return (parallel_segments.size() == 0 && other_segments.size() == 0);
}
