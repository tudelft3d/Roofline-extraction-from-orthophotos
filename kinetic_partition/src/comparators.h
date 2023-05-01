#pragma once

using std::vector;
using std::pair;

class pointcomp {
public:
	template <class point_type>
	bool operator()(const point_type * pl, const point_type * pr) const
	{
		if (pl->x < pr->x) { // the comparator should not contain >= or <=
			return 1;
		}
		else if (pl->x == pr->x) {
			return (pl->y < pr->y);
		}
		else {
			return 0;
		}
	}

	template <class point_type>
	bool operator()(const point_type pl, const point_type pr) const
	{
		if (pl.x < pr.x) {
			return 1;
		}
		else if (pl.x == pr.x) {
			return (pl.y < pr.y);
		}
		else {
			return 0;
		}
	}
};

class costcomp {
public:
	template <class T>
	bool operator()(const T * pl, const T * pr)  const
	{
		if (pl->support_cost < pr->support_cost) {
			return 1;
		}
		else if (pl->support_cost == pr->support_cost) {
			if ((pl->length > pr->length)) return 1;
			else if (pl->length == pr->length) return (pl->index < pr->index);
			else return 0;
		}
		else {
			return 0;
		}
	}
};

class queuecomp {
public:
	template <class T>
	bool operator()(const T pl, const T pr)  const
	{
		if (pl.second > pr.second) {
			return 1;
		}
		else if (pl.second == pr.second) {
			return (pl.first > pr.first);
		}
		else {
			return 0;
		}
	}
};

class bivalentcomp {
public:
	template <class T>
	bool operator()(const T pl, const T pr)  const
	{
		if (pl.second > pr.second) {
			return 1;
		}
		else if (pl.second == pr.second) {
			return (pl.first > pr.first);
		}
		else {
			return 0;
		}
	}
};

class bboxcomp_by_x {
public:
	template <class T>
	bool operator()(const T pl, const T pr)  const
	{
		if ( pl.second > pr.second) {
			return 1;
		}
		else if (pl.second == pr.second) {
			return (pl.first > pr.first);
		}
		else {
			return 0;
		}
	}
};

class bboxcomp_by_y {
public:
	template <class T>
	bool operator()(const T pl, const T pr)  const
	{
		if (pl.second > pr.second) {
		return 1;
		}
		else if (pl.second == pr.second) {
			return (pl.first > pr.first);
		}
		else {
			return 0;
		}
	}
};

class comp_by_second {
public:
	template <class T>
	bool operator()(const T pl, const T pr)  const
	{
		if (pl.second < pr.second) {
			return 1;
		}
		else if (pl.second == pr.second) {
			return (pl.first < pr.first);
		}
		else {
			return 0;
		}
	}
};

class comp_x {
public:
	template <class point_type>
	bool operator()(const point_type * pl, const point_type * pr) const
	{
		if (pl->pt.x < pr->pt.x) { // the comparator should not contain >= or <=
			return 1;
		}
		else {
			return 0;
		}
	}
};

class comp_y {
public:
	template <class point_type>
	bool operator()(const point_type * pl, const point_type * pr) const
	{
		if (pl->pt.y < pr->pt.y) { // the comparator should not contain >= or <=
			return 1;
		}
		else {
			return 0;
		}
	}
};