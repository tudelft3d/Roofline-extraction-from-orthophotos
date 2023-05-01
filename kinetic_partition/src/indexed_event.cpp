#include "indexed_event.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <ctime>


IndexedEvent::IndexedEvent(int _intersectant, int _intersected, double _t_intersectant, double _t_intersected, bool _is_colliding_colinear) {
	is_colliding_colinear = _is_colliding_colinear;
	is_colliding_ray = true;
	intersectant = _intersectant;
	intersected = _intersected;
	intersected_index = _intersected;
	t_intersectant = _t_intersectant;
	t_intersected = _t_intersected;
	prev_i = NULL;
	next_i = NULL;
	prev_j = NULL;
	next_j = NULL;
	previous = NULL;
	next = NULL;
}


IndexedEvent::IndexedEvent(int _intersectant, Image_Boundary _intersected, double _t_intersectant, double _t_intersected) {
	is_colliding_colinear = false;
	is_colliding_ray = false;
	intersectant = _intersectant;
	intersected = _intersected;
	intersected_index = -1 - _intersected;
	t_intersectant = _t_intersectant;
	t_intersected = _t_intersected;
	prev_i = NULL;
	next_i = NULL;
	prev_j = NULL;
	next_j = NULL;
	previous = NULL;
	next = NULL;
}

IndexedEvent::IndexedEvent(int _intersectant, Image_Boundary _intersected, int outer_edge_index, double _t_intersectant, double _t_intersected) {
	is_colliding_colinear = false;
	is_colliding_ray = false;
	intersectant = _intersectant;
	intersected = _intersected;
	intersected_index = outer_edge_index;
	t_intersectant = _t_intersectant;
	t_intersected = _t_intersected;
	prev_i = NULL;
	next_i = NULL;
	prev_j = NULL;
	next_j = NULL;
	previous = NULL;
	next = NULL;
}

IndexedEvent::IndexedEvent(const IndexedEvent & E)
{
	is_colliding_colinear = E.is_colliding_colinear;
	is_colliding_ray = E.is_colliding_ray;
	intersectant = E.intersectant;
	intersected = E.intersected;
	intersected_index = E.intersected_index;
	t_intersectant = E.t_intersectant;
	t_intersected = E.t_intersected;

	prev_i = NULL;
	next_i = NULL;
	prev_j = NULL;
	next_j = NULL;
	previous = NULL;
	next = NULL;
}


IndexedEvent::~IndexedEvent()
{
}


void IndexedEvent::clear_schedule(IndexedEvent* &schedule)
{
    IndexedEvent *event = schedule, *next_event = NULL;
    while (event != NULL) {
        next_event = event->next;
        delete event;
        event = next_event;
    }
    schedule = NULL;
}
