#pragma once
#include <iostream>
#include <list>
#include <vector>
#include "segment_ray.h"
#include "quadtree_point.h"

//using std::list;
using std::vector;


class Quadtree
{
public:
	Quadtree(double _x_min, double _x_max, double _y_min, double _y_max, Quadtree* _parent = NULL)
	{
		is_node = false;

		x_min = _x_min;
		x_max = _x_max;
		y_min = _y_min;
		y_max = _y_max;

		parent = _parent;

		north_west = NULL;
		north_east = NULL;
		south_west = NULL;
		south_east = NULL;
	}


    template<typename T>
    Quadtree(double _x_min, double _x_max, double _y_min, double _y_max, std::list<T *> & elements, Quadtree* _parent = NULL)
    {
        is_node = false;

        x_min = _x_min;
        x_max = _x_max;
        y_min = _y_min;
        y_max = _y_max;

		parent = _parent;

        north_west = NULL;
        north_east = NULL;
        south_west = NULL;
        south_east = NULL;

        for (typename std::list<T *>::iterator it = elements.begin() ; it != elements.end() ; it++) {
            this->add(dynamic_cast<Quadtree_Point*>(*it));
        }
    }


	template<typename T>
	Quadtree(double _x_min, double _x_max, double _y_min, double _y_max, vector<T *> & elements, Quadtree* _parent = NULL)
	{
		is_node = false;

		x_min = _x_min;
		x_max = _x_max;
		y_min = _y_min;
		y_max = _y_max;

		parent = _parent;
		
		north_west = NULL;
		north_east = NULL;
		south_west = NULL;
		south_east = NULL;

		for (int i = 0; i < int(elements.size()); i++) {
			this->add(dynamic_cast<Quadtree_Point*>(elements[i]));
		}
	}


    ~Quadtree() {
        if (is_node) {
            delete north_west;
            delete north_east;
            delete south_west;
            delete south_east;
        }
    }

public:
    template <typename T, typename U>
    void search(U & pt, double D, std::list<T *> & found_vertices)
    {
        double a_min = pt.x - D, a_max = pt.x + D;
        double b_min = pt.y - D, b_max = pt.y + D;
		if (((in(x_min, a_min, x_max) || in(x_min, a_max, x_max)) || (in(a_min, x_min, a_max) || in(a_min, x_max, a_max)))
			&& ((in(y_min, b_min, y_max) || in(y_min, b_max, y_max)) || (in(b_min, y_min, b_max) || in(b_min, y_max, b_max)))) {

			if (is_node) {
				// If the quadtree intercepts the area of search, and if this tree is a node,
				// we call the function recursively on each of the four subtrees
				this->north_west->search(pt, D, found_vertices);
				this->north_east->search(pt, D, found_vertices);
				this->south_west->search(pt, D, found_vertices);
				this->south_east->search(pt, D, found_vertices);

			} else {
				// Otherwise, if the quadtree is a leaf, we loop on the list of points
				// If a vertex is contained in the area of search, then we append it to the list of vertices
                for (std::list<Quadtree_Point *>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
					Quadtree_Point* v = (*it_v);
					if (norm(Vec2d(v->keypoint().x - pt.x, v->keypoint().y - pt.y)) < D) {
						found_vertices.push_back(dynamic_cast<T *>(v));
					}
				}
			}
		}
    }


	template <typename T>
    void search(double _x_min, double _x_max, double _y_min, double _y_max, std::list<T *> & found_vertices)
	{
		// First we check that the area covered by the quadtree intercepts the rectangle defined by parameters
		if (((in(x_min, _x_min, x_max) || in(x_min, _x_max, x_max)) || (in(_x_min, x_min, _x_max) || in(_x_min, x_max, _x_max)))
			&& ((in(y_min, _y_min, y_max) || in(y_min, _y_max, y_max)) || (in(_y_min, y_min, _y_max) || in(_y_min, y_max, _y_max)))) {

			if (is_node) {
				// If the quadtree intercepts the area of search, and if this tree is a node,
				// we call the function recursively on each of the four subtrees
				this->north_west->search(_x_min, _x_max, _y_min, _y_max, found_vertices);
				this->north_east->search(_x_min, _x_max, _y_min, _y_max, found_vertices);
				this->south_west->search(_x_min, _x_max, _y_min, _y_max, found_vertices);
				this->south_east->search(_x_min, _x_max, _y_min, _y_max, found_vertices);

			} else {
				// Otherwise, if the quadtree is a leaf, we loop on the list of points
				// If a vertex is contained in the area of search, then we append it to the list of vertices
                for (std::list<Quadtree_Point *>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
					Quadtree_Point* v = (*it_v);
					Point2d pt = v->keypoint();
					if (_x_min <= pt.x && pt.x <= _x_max && _y_min <= pt.y && pt.y <= _y_max) {
						found_vertices.push_back(dynamic_cast<T *>(v));
					}
				}
			}
		}
	}


    void add(Quadtree_Point* v0)
    {
        if (is_node) {
            // If the current object is a node, we insert the vertex into the right subtree
            this->assign(v0->keypoint())->add(v0);

        } else {
			double dx = x_max - x_min;
			double dy = y_max - y_min;
			bool capacity_exceeded = (int(vertices.size()) > 20);
            
			// If the current node hasn't exceed its maximal capacity, or if it can no longer be divided,
            // we insert the vertex
            if (!capacity_exceeded || (capacity_exceeded && (dy < 0.5 || dx < 0.5))) {
                vertices.push_back(v0);

            } else {
                is_node = true;
                // Otherwise, we divide the leaf into four subtrees
                north_west = new Quadtree(x_min, x_min + dx / 2, y_min + dy / 2, y_max, this);
                north_east = new Quadtree(x_min + dx / 2, x_max, y_min + dy / 2, y_max, this);
                south_west = new Quadtree(x_min, x_min + dx / 2, y_min, y_min + dy / 2, this);
                south_east = new Quadtree(x_min + dx / 2, x_max, y_min, y_min + dy / 2, this);

                // We divide the list of points and assign them to one of the subtrees according to their coordinates
                for (std::list<Quadtree_Point *>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
                    Quadtree_Point* v = (*it_v);
                    this->assign(v->keypoint())->add(v);
                }
                vertices.clear();

                // Finally, we add v0 to the corresponding node
                this->assign(v0->keypoint())->add(v0);
            }
        }
    }


	void remove(Quadtree_Point* v, bool destroy)
	{
		// This function should be called on the root of the subtree
		if (parent != NULL) return;

		// First we search for the subtree where v is located
		// Once it is done, we remove v from the list of vertices
		
		const Point2d & pt = v->keypoint();
		Quadtree* q = this;
		while (q->is_node) {
			q = q->assign(pt);
		}

        for (std::list<Quadtree_Point *>::iterator it_v = q->vertices.begin(); it_v != q->vertices.end(); it_v++) {
			if ((*it_v) == v) {
				if (destroy) delete v;
				it_v = q->vertices.erase(it_v);
				break;
			}
		}

		// Second, we check if q is now an empty leaf
		if (q->is_empty_leaf()) {
			Quadtree* q_p = q->parent;
			// - q_curr == NULL : we have reached the root of the quadtree
			// - q_prev == q_curr : q_curr has not been deleted because one of the four children is not an empty leaf
			while (q_p != NULL && q != q_p) {
				q = remove(q_p);
				if (q != NULL) {
					q_p = q->parent;
				} else {
					q_p = NULL;
				}
			}
		}
	}

private:
    template<typename U>
    Quadtree* assign(const U & pt)
    {
        double dx = x_max - x_min;
        double dy = y_max - y_min;
        if (pt.x < x_min + double(dx) / 2) {
            if (pt.y < y_min + double(dy) / 2) {
                return south_west;
            } else {
                return north_west;
            }
        } else {
            if (pt.y < y_min + double(dy) / 2) {
                return south_east;
            } else {
                return north_east;
            }
        }
    }


	static Quadtree* remove(Quadtree* q)
	{
		assert(q != NULL);

		// Given a quadtree, we check if it is a node with four empty leaves
		// If so, this quadtree is turned into a leaf, and we return a pointer 
		// to this quadtree's parent, so that the same check can be performed on it
		if (q->is_node && q->children_are_empty_leaves()) {

			// If the four subtrees of q are empty leaves,
			// then we obtain a pointer to the parent tree
			Quadtree* q_p = q->parent;

			// Deletes the subtrees and returns a pointer to the parent tree
			delete q->north_west;
			delete q->north_east;
			delete q->south_west;
			delete q->south_east;
			q->is_node = false;
			q->north_west = q->north_east = q->south_west = q->south_east = NULL;
			return q_p;
		} else {
			// We do nothing
			return q;
		}
	}


	inline bool is_empty_leaf() {
		return (!is_node && (vertices.size() == 0));
	}


	inline bool children_are_empty_leaves() {
		return (north_west->is_empty_leaf() && north_east->is_empty_leaf() && south_west->is_empty_leaf() && south_east->is_empty_leaf());
	}


    inline bool in(double a, double x, double b) {
        return ((a <= x) && (x <= b));
    }

public:
    bool is_node;

    double x_min;
    double x_max;
    double y_min;
    double y_max;

	Quadtree* parent;

    Quadtree* north_west;
    Quadtree* north_east;
    Quadtree* south_west;
    Quadtree* south_east;

    std::list<Quadtree_Point *> vertices;
};
