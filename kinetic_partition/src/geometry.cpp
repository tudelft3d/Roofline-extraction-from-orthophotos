#include "defs.h"
#include "geometry.h"
#include "segment_tree.h"
#include "partition_elements.h"
#include <iostream>

int Geometry::dcmp(double x) {
	if (std::fabs(x) < 1e-6) return 0;
	else return x<0? -1:1;
}

double Geometry::mes(Vec2d & _a, Vec2d & _b) {
    Vec2d a = normalize(_a);
    Vec2d b = normalize(_b);
    double x = a.ddot(b);
    if (x > 1.0 - 1e-7) {
        return 0.0;
    } else if (x < -1.0 + 1e-7) {
        return PI;
    } else {
        return acos(x);
    }
}


double Geometry::distance_initial_coordinates(Segment *s, Segment *t) {
    Vec2d & s_dir = s->direction;
    Vec2d & t_dir = t->direction;
    Point2d & A = s->end1;
    Point2d & B = s->end2;
    Point2d & C = t->end1;
    Point2d & D = t->end2;

    Vec2d s_nor(-s_dir[1], s_dir[0]);
    if (fabs(s_nor.ddot(t_dir)) > 1e-6) {
		// Computes the intersection
		double det = (B.x - A.x) * (C.y - D.y) - (B.y - A.y) * (C.x - D.x);
		double t_AB = ((C.x - A.x) * (C.y - D.y) - (C.y - A.y) * (C.x - D.x)) / det;
		double u_CD = ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) / det;
		if (t_AB >= 0 && t_AB <= 1 && u_CD >= 0 && u_CD <= 1) {
			// Intersection belongs to [AB] and [CD] : distance is null
			return 0;
		}
	}
	Vec2d AB = B - A, AC = C - A, AD = D - A;
	Vec2d CA = A - C, CB = B - C, CD = D - C;
	Vec2d BC = C - B, BD = D - B;
	double l_AB2 = norm(AB) * norm(AB);
	double l_CD2 = norm(CD) * norm(CD);
	double t_C = AC.ddot(AB) / l_AB2;
	double t_D = AD.ddot(AB) / l_AB2;
	double t_A = CA.ddot(CD) / l_CD2;
	double t_B = CB.ddot(CD) / l_CD2;
    double dist = MIN(MIN(norm(AC), norm(AD)), MIN(norm(BC), norm(BD)));
	if (t_C >= 0 && t_C <= 1) {
		Point2d C_prime = A + Point2d(t_C * AB);
        dist = MIN(dist, norm(C_prime - C));
	}
	if (t_D >= 0 && t_D <= 1) {
		Point2d D_prime = A + Point2d(t_D * AB);
        dist = MIN(dist, norm(D_prime - D));
	}
	if (t_A >= 0 && t_A <= 1) {
		Point2d A_prime = C + Point2d(t_A * CD);
        dist = MIN(dist, norm(A_prime - A));
	}
	if (t_B >= 0 && t_B <= 1) {
		Point2d B_prime = C + Point2d(t_B * CD);
        dist = MIN(dist, norm(B_prime - B));
	}
	return dist;
}


double Geometry::distance(cv::Point2d &p, Segment *s)
{
	return (fabs(s->a * p.x + s->b * p.y + s->c) / sqrt(s->a * s->a + s->b * s->b));
}


void Geometry::build_tree_polygon_boundaries(vector<Segment *> & segments, vector<list<HalfEdge *>> & polygon_edges, Segment_Regularization_Tree * tree, 
	map<int, double> & reg_angles, vector<int> & bin_to_angle_index, map<HalfEdge *, list<Segment *>> & boundary_collinear_segments)
{
	tree->delete_parallel_nodes();
	map<int, list<Segment *> > groups;
	map<int, set<pair<HalfEdge *, double>, comp_by_second>> group_edges;
	map<Node_Parallel_Segments *, set<pair<HalfEdge *, double>, comp_by_second>> node_edges;
	uint n = uint(segments.size());

	int alpha_n_bins = bin_to_angle_index.size();
	double alpha_bin_width = 180. / double(alpha_n_bins); // default: 2 degrees
	double y_bin_width = 1.5;

	// Part 1. Discretization of the set of observed orientations
	
	// First, we associate each segment to a group of parallel segments, indexed by quantized angle starting from 0
	map<int, double> angles;
	for (uint i = 0; i < n; ++i) {
		// unit: degree, in [0, 180]
		// Tries to match a segment's orientation with another quantized angle
		int ind_bin = int(round(segments[i]->alpha / alpha_bin_width)) % alpha_n_bins;
		if (ind_bin < 0) ind_bin += alpha_n_bins;
		int related_group = bin_to_angle_index[ind_bin] < 0 ? ind_bin : bin_to_angle_index[ind_bin];
		// add angle of the bin center
		if (reg_angles.count(related_group)==0 || bin_to_angle_index[ind_bin] < 0) angles[ind_bin] = ind_bin*alpha_bin_width;
		else {
			angles[related_group] = reg_angles[related_group];
		}

		groups[related_group].push_back(segments[i]);
	}

	// Second, we associate each boundary edges to, if exist, one parallel group
	for (int i = 0; i < polygon_edges.size(); ++i) {
		for (list<HalfEdge *>::iterator it_e = polygon_edges[i].begin(); it_e != polygon_edges[i].end(); ++it_e) {
			double alpha = (*it_e)->e->get_alpha_1() * 180 / PI;
			if (alpha < 0) alpha += 180;
			if (alpha > 180) alpha -= 180;
			double d_angle = 180;
			int ind_bin = int(round(alpha / alpha_bin_width)) % alpha_n_bins;
			if (ind_bin < 0) ind_bin += alpha_n_bins;
			
			// Tries to match edge orientation with an angle group
			int related_group = bin_to_angle_index[ind_bin];
			for (int k = -1; k <= 1; ++k) {
				d_angle = fabs(alpha - angles[related_group] + k * 180);
				if (d_angle < 0.5*alpha_bin_width) {
					break;
				}
			}
			if (d_angle < 0.5*alpha_bin_width) {
				group_edges[related_group].insert(make_pair(*it_e, d_angle));
			}
		}
	}

	// Builds the first level of the tree
	for (map<int, list<Segment *> >::iterator it_m = groups.begin(); it_m != groups.end(); ++it_m) {
		double theta = angles[it_m->first];

		Vec2d v_dir = Vec2d(cos(theta * PI / 180), sin(theta * PI / 180));
		Vec2d v_nor = Vec2d(-v_dir[1], v_dir[0]);
		double a = v_nor[0], b = v_nor[1];

		tree->create_parallel_node(theta);
		for (list<Segment *>::iterator it_s = it_m->second.begin(); it_s != it_m->second.end(); ++it_s) {
			Segment* s = (*it_s);

			// Computes the equation of the support line of s
			double c = -a * s->barycenter.x - b * s->barycenter.y;
			s->set_dalpha(theta - s->alpha, theta, a, b, c, v_dir);
			tree->assign_to_parallel_node(theta, s);
		}
		if (group_edges.count(it_m->first) != 0) {
			node_edges[tree->parallel_segments[theta]] = group_edges[it_m->first];
		}
	}

	// Part 2. Detects near-colinear segments

	for (map<double, Node_Parallel_Segments *>::iterator it_c = tree->parallel_segments.begin(); it_c != tree->parallel_segments.end(); ++it_c) {

		groups.clear();
		group_edges.clear();
		double theta = it_c->first;
		Node_Parallel_Segments* cluster = it_c->second;
		
		// We transform the coordinates of the barycenters of all segments in the local frame defined by :
		// - the center of the first segment in the cluster
		// - the unit vectors I(cos(theta), sin(theta)) and J(-sin(theta), cos(theta))
		Point2d O = cluster->parallel_segments.front()->interBarycenter;
		Vec2d I = Vec2d(cos(theta * PI / 180), sin(theta * PI / 180));
		if (I[1] < 0 || (I[1] == 0 && I[0] < 0)) I = -I;
		Vec2d J = Vec2d(-I[1], I[0]);

		map<int, double> ordinates;
		map<int, int> map_bin_to_group;
		for (list<Segment *>::iterator it_s = cluster->parallel_segments.begin(); it_s != cluster->parallel_segments.end(); ++it_s) {
			Segment* s = (*it_s);
			double dx = s->interBarycenter.x - O.x, dy = s->interBarycenter.y - O.y;
			double y = dx * J[0] + dy * J[1];

			// Following a process directly adapted from the first part of this function,
			// we want to merge segments whose y is too close to other quantized values of y computed before (near-colinear segments)
			int related_group = -1;

			int ind_bin = int(round(y / y_bin_width)); // it can be negative
			if (map_bin_to_group.count(ind_bin) != 0) {
				related_group = map_bin_to_group.at(ind_bin);
			}

			if (related_group == -1) {
				related_group = groups.size();
				map_bin_to_group[ind_bin] = related_group;
				// add quantized y value
				ordinates[related_group] = ind_bin*y_bin_width;
			}
			groups[related_group].push_back(s);
		}

		// We also transform the coordinates of the related edge centers to the above local frame
		map<int, Vec2d> edge_frame_I;
		map<int, Point2d> edge_frame_center;
		if (node_edges.count(cluster) != 0) {
			/* find the first edge with angle deviation less than 0.05 */
			set<pair<HalfEdge *, double>, comp_by_second> related_edges = node_edges.at(cluster);
			auto it_e = related_edges.begin();
			while (it_e != related_edges.end()) {
				bool removed = false;
				if (it_e->second > 0.05) {
					++it_e; continue;
				}
				Vertex* v1 = it_e->first->e->v1, *v2 = it_e->first->e->v2;
				Point2d edge_center = (v1->pt + v2->pt)/2;
				double dx = edge_center.x - O.x, dy = edge_center.y - O.y;

				double y = dx * J[0] + dy * J[1];
				// Align ordinates of collinear segments with edge ordinate
				double y_d = y_bin_width / 3;
				int related_group = -1;
				for (uint j = 0; j < groups.size(); ++j) {
					double y_j = ordinates[j];
					if (fabs(y_j - y) < y_d) {
						y_d = fabs(y_j - y);
						related_group = j; break;
					}
				}
				if (related_group != -1) {
					if (group_edges.count(related_group) != 0) {
						// There already exists a collinear edge for this collinear group
						Edge * e_selected = group_edges.at(related_group).begin()->first->e;
						// save the extra edge if it is collinear with the selected edge
						for (int k = -2; k <= 2; ++k) {
							if (fabs(it_e->first->e->get_alpha_1() - e_selected->get_alpha_1() + k * PI) < 1e-4) {
								Point2d center_selected = edge_frame_center[related_group];
								Vec2d normal(-edge_frame_I[related_group][1], edge_frame_I[related_group][0]);
								if (fabs(normal.dot(center_selected - v1->pt)) < 1e-3 && fabs(normal.dot(center_selected - v2->pt)) < 1e-3) {
									group_edges[related_group].insert(make_pair(it_e->first, y_d));
									boundary_collinear_segments[it_e->first] = groups[related_group];
									// remove from related_edges
									it_e = related_edges.erase(it_e);
									removed = true;
								}
								break;
							}
						}
					}
					else {
						set<pair<HalfEdge *, double>, comp_by_second> related_edge;
						related_edge.insert(make_pair(it_e->first, y_d));
						group_edges[related_group] = related_edge;
						boundary_collinear_segments[it_e->first] = groups[related_group];

						// align segment orientation with the selected edge
						double alpha = it_e->first->e->get_alpha_1();
						while (alpha < 0) alpha += PI;
						while (alpha > PI) alpha -= PI;
						Vec2d v_dir = Vec2d(cos(alpha), sin(alpha));
						if (v_dir[1] < 0 || (v_dir[1] == 0 && v_dir[0] < 0)) v_dir = -v_dir;
						Vec2d v_nor = Vec2d(-v_dir[1], v_dir[0]);

						// rotate segments
						double new_theta = alpha * 180 / PI;
						double a = v_nor[0], b = v_nor[1];
						for (list<Segment *>::iterator it_s = groups[related_group].begin(); it_s != groups[related_group].end(); ++it_s) {
							Segment* s = (*it_s);
							// Computes the equation of the support line of s
							double c = -a * s->interBarycenter.x - b * s->interBarycenter.y;
							s->set_dalpha(new_theta - s->alpha, new_theta, a, b, c, v_dir);
						}
						// Record the new local frame (keep the center O unchanged)
						edge_frame_center[related_group] = edge_center;
						//edge_frame_I[related_group] = (*groups[related_group].begin())->finalDirection;
						edge_frame_I[related_group] = v_dir;
						// remove from related_edges
						it_e = related_edges.erase(it_e);
						removed = true;
					}
				}
				if (!removed) ++it_e;
			}
			// add the remaining edges to the collinear group
			it_e = related_edges.begin();
			while (it_e != related_edges.end()) {
				Vertex* v1 = it_e->first->e->v1, *v2 = it_e->first->e->v2;
				Point2d edge_center = (v1->pt + v2->pt) / 2;
				double dx = edge_center.x - O.x, dy = edge_center.y - O.y;
				double y = dx * J[0] + dy * J[1];
				// Align ordinates of collinear segments with edge ordinate
				double y_d = y_bin_width / 2;
				int related_group = -1;
				for (uint j = 0; j < groups.size(); ++j) {
					double y_j = ordinates[j];
					if (fabs(y_j - y) < y_d) {
						y_d = fabs(y_j - y);
						related_group = j; break;
					}
				}
				if (related_group != -1 && group_edges.count(related_group) != 0) {
					Edge * e_selected = group_edges.at(related_group).begin()->first->e;
					// save the extra edge if it is collinear with the selected edge
					for (int k = -2; k <= 2; ++k) {
						if (fabs(it_e->first->e->get_alpha_1() - e_selected->get_alpha_1() + k * PI) < 1e-4) {
							Point2d center_selected = edge_frame_center[related_group];
							Vec2d normal(-edge_frame_I[related_group][1], edge_frame_I[related_group][0]);
							if (fabs(normal.dot(center_selected - v1->pt)) < 1e-3 && fabs(normal.dot(center_selected - v2->pt)) < 1e-3) {
								group_edges[related_group].insert(make_pair(it_e->first, y_d));
								boundary_collinear_segments[it_e->first] = groups[related_group];
								// remove from related_edges
								it_e = related_edges.erase(it_e);
							}
							break;
						}
					}
				}
				++it_e;
			}
		}

		// Builds the second level of the tree
		for (map<int, list<Segment *> >::iterator it_m = groups.begin(); it_m != groups.end(); ++it_m) {

			double dt = ordinates[it_m->first];
			Vec2d local_J = Vec2d(-I[1], I[0]);
			cluster->create_colinear_node(dt);

			if (group_edges.count(it_m->first) != 0) {
				Vec2d local_I = edge_frame_I.at(it_m->first);
				local_J = Vec2d(-local_I[1], local_I[0]);
				Point2d edge_center = edge_frame_center.at(it_m->first);
				// Translates all segments
				for (list<Segment *>::iterator it_s = it_m->second.begin(); it_s != it_m->second.end(); ++it_s) {
					Segment* s = (*it_s);
					double dx = s->interBarycenter.x - edge_center.x, dy = s->interBarycenter.y - edge_center.y;
					double y = dx * local_J[0] + dy * local_J[1];
					if (fabs(-y) >= y_bin_width) {
						std::cout << "edge center: " << edge_center << std::endl;
						std::cout << "local_J: " << local_J << std::endl;
						std::cout << "segment center " << s->interBarycenter << std::endl;
						std::cout << "y = " << y << std::endl;
					}
					assert(fabs(-y) <= y_bin_width);
					s->set_dt(- y);
				}
			}
			else {
				// Gets the longest segment
				double l_max = -FLT_MAX;
				Segment* s_longest = NULL;
				for (list<Segment *>::iterator it_s = it_m->second.begin(); it_s != it_m->second.end(); ++it_s) {
					if ((*it_s)->length > l_max) {
						l_max = (*it_s)->length;
						s_longest = (*it_s);
					}
				}
				// Translates the longest segment and gets its line equation
				double dx = s_longest->interBarycenter.x - O.x, dy = s_longest->interBarycenter.y - O.y;
				double y = dx * local_J[0] + dy * local_J[1];
				s_longest->set_dt(dt - y);
				double a = s_longest->a, b = s_longest->b, c = s_longest->c;
				Vec2d dir = s_longest->finalDirection;

				// Translates other segments
				for (list<Segment *>::iterator it_s = it_m->second.begin(); it_s != it_m->second.end(); ++it_s) {
					Segment* s = (*it_s);
					if (s == s_longest) continue;

					dx = s->interBarycenter.x - O.x, dy = s->interBarycenter.y - O.y;
					y = dx * local_J[0] + dy * local_J[1];
					s->set_dt(dt - y, a, b, c, dir);
				}
			}
			
			cluster->assign_to_colinear_node(dt, it_m->second);
		}
	}

}


void Geometry::build_tree_when_regularization_disabled(vector<Segment *> & segments, Segment_Regularization_Tree * tree)
{
	map<int, list<Segment *> > groups;
	uint n = uint(segments.size());

	int alpha_n_bins = 90;
	double alpha_bin_width = 180./double(alpha_n_bins); // default: 2 degrees
	
	double y_bin_width = 0.9;
	// Part 1. Discretization of the set of observed orientations

	// First, we associate each segment to a group of parallel segments, indexed by quantized angle starting from 0
	vector<double> angles;
	vector<int> bin_to_group_index(alpha_n_bins, -1);
	for (uint i = 0; i < n; ++i) {
		int related_group = -1;
		double alpha = segments[i]->alpha; // unit: degree, in [0, 180]
		// Tries to match a segment's orientation with another quantized angle
		int ind_bin = int(round(alpha / alpha_bin_width)) % alpha_n_bins;
		if (ind_bin < 0) ind_bin += alpha_n_bins;
		if (bin_to_group_index[ind_bin] != -1) {
			related_group = bin_to_group_index[ind_bin];
		}

		if (related_group == -1) {
			// Creates a new group with alpha_i
			related_group = int(groups.size());
			bin_to_group_index[ind_bin] = related_group;
			// add angle of the bin center
			angles.push_back(ind_bin*alpha_bin_width);
		}
		groups[related_group].push_back(segments[i]);
	}

	// Builds the first level of the tree
	for (map<int, list<Segment *> >::iterator it_m = groups.begin() ; it_m != groups.end() ; ++it_m) {
		double theta = angles[it_m->first];

		Vec2d v_dir = Vec2d(cos(theta * PI / 180), sin(theta * PI / 180));
        Vec2d v_nor = Vec2d(-v_dir[1], v_dir[0]);
        double a = v_nor[0], b = v_nor[1];

		tree->create_parallel_node(theta);
		for (list<Segment *>::iterator it_s = it_m->second.begin() ; it_s != it_m->second.end() ; ++it_s) {
			Segment* s = (*it_s);

			// Computes the equation of the support line of s
			double c = -a * s->barycenter.x - b * s->barycenter.y;
            s->set_dalpha(theta - s->alpha, theta, a, b, c, v_dir);
			tree->assign_to_parallel_node(theta, s);
		}
	}

	// Part 2. Detects near-colinear segments

	for (map<double, Node_Parallel_Segments *>::iterator it_c = tree->parallel_segments.begin() ; it_c != tree->parallel_segments.end() ; ++it_c) {

		groups.clear();
		double theta = it_c->first;
		Node_Parallel_Segments* cluster = it_c->second;
		Point2d O = cluster->parallel_segments.front()->interBarycenter;

		// We transform the coordinates of the barycenters of all segments in the local frame defined by :
		// - the center of the first segment in the cluster
		// - the unit vectors I(cos(theta), sin(theta)) and J(-sin(theta), cos(theta))
		Vec2d I = Vec2d(cos(theta * PI / 180), sin(theta * PI / 180));
        Vec2d J = Vec2d(-I[1], I[0]);

		vector<double> ordinates;
		map<int, int> map_bin_to_group;
		for (list<Segment *>::iterator it_s = cluster->parallel_segments.begin(); it_s != cluster->parallel_segments.end(); ++it_s) {
			Segment* s = (*it_s);
			double dx = s->interBarycenter.x - O.x, dy = s->interBarycenter.y - O.y;
			double y = dx * J[0] + dy * J[1];

			// Following a process directly adapted from the first part of this function,
			// we want to merge segments whose y is too close to other quantized values of y computed before (near-colinear segments)
			int related_group = -1;

			int ind_bin = int(round(y / y_bin_width)); // it can be negative
			if (map_bin_to_group.count(ind_bin) != 0) {
				related_group = map_bin_to_group.at(ind_bin);
			}

			if (related_group == -1) {
				related_group = groups.size();
				map_bin_to_group[ind_bin] = related_group;
				// add quantized y value
				ordinates.push_back(ind_bin*y_bin_width);
			}
			groups[related_group].push_back(s);
		}

		//for (list<Segment *>::iterator it_s = cluster->parallel_segments.begin() ; it_s != cluster->parallel_segments.end() ; ++it_s) {
		//	Segment* s = (*it_s);
		//	double dx = s->interBarycenter.x - O.x, dy = s->interBarycenter.y - O.y;
		//	double y = dx * J[0] + dy * J[1];

		//	// Following a process directly adapted from the first part of this function,
		//	// we want to merge segments whose y is too close to other values of y computed before (near-colinear segments)
		//	int related_group = -1;
		//	for (uint j = 0 ; j < groups.size() ; ++j) {
		//		double y_j = ordinates[j];
		//		if (fabs(y_j - y) < y_eps) {
		//			related_group = j;
		//			break;
		//		}
		//		if (related_group != -1) break;
		//	}

		//	if (related_group == -1) {
		//		related_group = groups.size();
		//		ordinates.push_back(y);
		//	}
		//	groups[related_group].push_back(s);
		//}

		// Builds the second level of the tree
		for (map<int, list<Segment *> >::iterator it_m = groups.begin() ; it_m != groups.end() ; ++it_m) {
			
			double dt = ordinates[it_m->first];

			// Gets the longest segment
			double l_max = -FLT_MAX;
			Segment* s_longest = NULL;
			for (list<Segment *>::iterator it_s = it_m->second.begin() ; it_s != it_m->second.end() ; ++it_s) {
				if ((*it_s)->length > l_max) {
					l_max = (*it_s)->length;
					s_longest = (*it_s);
				}
			}

			cluster->create_colinear_node(dt);

			// Translates the longest segment and gets its line equation
			double dx = s_longest->interBarycenter.x - O.x, dy = s_longest->interBarycenter.y - O.y;
			double y = dx * J[0] + dy * J[1];

			s_longest->set_dt(dt - y);
			double a = s_longest->a, b = s_longest->b, c = s_longest->c;
			Vec2d dir = s_longest->finalDirection;

			// Translates other segments
			for (list<Segment *>::iterator it_s = it_m->second.begin() ; it_s != it_m->second.end() ; ++it_s) {
				Segment* s = (*it_s);
				if (s == s_longest) continue;

				dx = s->interBarycenter.x - O.x, dy = s->interBarycenter.y - O.y;
				y = dx * J[0] + dy * J[1];
				s->set_dt(dt - y, a, b, c, dir);
			}
			cluster->assign_to_colinear_node(dt, it_m->second);
		}
	}
}


void Geometry::merge_overlapped_segments(Segment_Regularization_Tree *tree)
{
	for (map<double, Node_Parallel_Segments*>::iterator it_m1 = tree->parallel_segments.begin() ; it_m1 != tree->parallel_segments.end() ; it_m1++) {
		Node_Parallel_Segments* node_parallel = it_m1->second;
		for (map<double, Node_Colinear_Segments*>::iterator it_m2 = node_parallel->colinear_segments.begin() ; it_m2 != node_parallel->colinear_segments.end() ; it_m2++) {
			Node_Colinear_Segments* node_colinear = it_m2->second;
			if (node_colinear->colinear_segments.size() == 1) continue;

			bool at_least_one_merge;
			do {
				at_least_one_merge = false;
				list<Segment *>::iterator it_s1, it_s2;
				for (it_s1 = node_colinear->colinear_segments.begin() ; it_s1 != node_colinear->colinear_segments.end() ; ++it_s1) {
					Segment* s1 = (*it_s1);
					if (s1->is_disabled) {
						continue;
					}
					it_s2 = it_s1;
					while (++it_s2 != node_colinear->colinear_segments.end()) {
						Segment* s2 = (*it_s2);
						if (s2->is_disabled) continue;
						Point2d P, Q;
						// Tests if colinear segments s1 and s2 are overlapping
						if (are_overlapping(s1, s2, P, Q)) {
							// s1 takes the coordinates of the merged segment s1 U s2
							s1->set_final_extrema(P, Q);
							if (s2->support_cost < s1->support_cost) s1->support_cost = s2->support_cost;
							s2->disable();
							s2->to_delete = true;
							at_least_one_merge = true;
						}
					}
				}
			} while (at_least_one_merge);
		}
	}
}


void Geometry::disable_segments_outside_boundaries(vector<Segment *> & segments, int rows, int cols)
{
	uint n = uint(segments.size());
	for (uint i = 0 ; i < n ; i++) {
		Point2d & O = segments[i]->finalBarycenter;
		if (O.x <= -0.5 || O.x >= cols-0.5 || O.y <= -0.5 || O.y >= rows-0.5) {
			segments[i]->disable();
		}
	}
}

void Geometry::disable_segments_outside_polygon_boundaries(vector<Segment *> & segments, vector<list<HalfEdge *>> & polygon_edges)
{
	// Removes segments outside polygon boundaries.

	uint n = uint(segments.size());
	for (uint i = 0; i < n; i++) {
		Segment* s_i = segments[i];
		Point2d & s_i_O = s_i->finalBarycenter;

		// Check if the segment finalcenter lies outside the polygon if it's not already disabled
		if (!s_i->is_disabled) {
			bool is_inside_polygon = is_inside_contour(s_i_O, polygon_edges[0]);
			if (is_inside_polygon) { 
				// Continue to check the inner boundaries. The segment is outside the polygon if it lies inside any inner contour.
				for (int id = 1; id < polygon_edges.size(); ++id) {
					if (is_inside_contour(s_i_O, polygon_edges[id])) {
						is_inside_polygon = false;
						break;
					}
				}
			}
			if (!is_inside_polygon) {
				s_i->disable();
			}
		}
	}

}


bool Geometry::OnSegment(Point2d P1,Point2d P2,Point2d Q)
{
    return dcmp((P1-Q).cross(P2-Q))==0&&dcmp((P1-Q).dot(P2-Q))<=0;
}

bool Geometry::is_inside_contour(Point2d & pt, list<HalfEdge *> & contour_edges) {
	bool flag = false;
	for (list<HalfEdge*>::iterator it = contour_edges.begin(); it != contour_edges.end(); ++it) {
		Vertex* A = (*it)->e->v1, *B = (*it)->e->v2;
		if(OnSegment(A->pt,B->pt,pt)) return true;
		if ( (dcmp(A->pt.y-pt.y)>0 != dcmp(B->pt.y-pt.y)>0) && dcmp(pt.x - (pt.y-A->pt.y)*(A->pt.x-B->pt.x)/(A->pt.y-B->pt.y)-A->pt.x)<0 ) 
			flag = !flag;
	}
	return flag;
}

bool Geometry::is_inside_contour(Point2d & pt, list<Point2d*> & contour_pts) {
	bool flag = false;
	for (list<Point2d*>::iterator it = contour_pts.begin(); it != contour_pts.end(); ++it) {
		Point2d* A = (*it);
		Point2d* B = std::next(it) == contour_pts.end() ? *contour_pts.begin(): *std::next(it);
		if(OnSegment(*A,*B,pt)) return true;
		if ( (dcmp(A->y-pt.y)>0 != dcmp(B->y-pt.y)>0) && dcmp(pt.x - (pt.y-A->y)*(A->x-B->x)/(A->y-B->y)-A->x)<0 ) 
			flag = !flag;
	}
	return flag;
}

bool Geometry::is_between(double alpha, double alpha_1, double alpha_2) {
	if (std::abs(alpha_1 - alpha_2) > PI) {
		if (alpha >= MAX(alpha_1, alpha_2) - 1e-5 || alpha <= MIN(alpha_1, alpha_2) + 1e-5) return true;
		else return false;
	}
	else {
		if (alpha >= MIN(alpha_1, alpha_2) - 1e-5 && alpha <= MAX(alpha_1, alpha_2) + 1e-5) return true;
		else return false;
	}
}

void Geometry::disable_segments_overlapping_polygon_boundaries(vector<Segment *> & segments, vector<list<HalfEdge *>> & polygon_edges, double eps)
{
	// 1. Disable segments overlapping with polygon boundaries. We consider a segment to be overlapping with a boundary edge if 
	// the furthest distance between the segment endpoints and the edge is less than the minimal edge length threshold, 
	// and the segment center projects within the edge.
	double collinear_eps = eps;

	for (uint i = 0; i < segments.size(); ++i) {
		Segment* s_i = segments[i];
		Point2d & s_i_O = s_i->finalBarycenter;
		Vec2d & s_i_dir = s_i->finalDirection;
		Vec2d s_i_normal(-s_i_dir[1], s_i_dir[0]);
		list<pair<HalfEdge *, pair<double, double>>> collinear_edges;

		for (int id = 0; id < polygon_edges.size(); ++id) {
			if (s_i->is_disabled) break;
			for (list<HalfEdge *>::iterator it_h = polygon_edges[id].begin(); it_h != polygon_edges[id].end(); ++it_h) {
				double e_length = (*it_h)->e->length;
				Point2d & p1 = (*it_h)->v1_v2 ? (*it_h)->e->v1->pt : (*it_h)->e->v2->pt;
				Point2d & p2 = (*it_h)->v1_v2 ? (*it_h)->e->v2->pt : (*it_h)->e->v1->pt;
				Vec2d & e_normal = (*it_h)->e->normal;

				// Check the first overlapping condition
				bool aligned_normal = (fabs(s_i_normal[0] * e_normal[0] + s_i_normal[1] * e_normal[1]) > 0.9848); // within +-10 degrees
				if (aligned_normal) {
					Point2d & s_i_A = s_i->finalEnd1;
					Point2d & s_i_B = s_i->finalEnd2;
					double dist_B = fabs(e_normal.dot(s_i_B - p1));
					double dist_A = fabs(e_normal.dot(s_i_A - p1));
					// Check the distance from segment center to the edge
					if (dist_A <= eps && dist_B <= eps) {
						Vec2d e_dir = (p2 - p1) / e_length;
						double axial_dist_1 = e_dir.dot(s_i_A - p1), axial_dist_2 = e_dir.dot(s_i_B - p1);
						double axial_dist_O = 0.5*(axial_dist_1 + axial_dist_2);
						// Check if the segment projects in between the edge vertices.
						if (axial_dist_O > 0 && axial_dist_O < e_length) {
							s_i->disable();
							if ((axial_dist_1 >= 0 && axial_dist_2 <= e_length) || (axial_dist_1 <= e_length && axial_dist_2 >= 0)) {
								break;
							}
							if (axial_dist_1 < 0 || axial_dist_1 > e_length) { // create cropped version to replace s_i
								double x1 = s_i->end1.x, y1 = s_i->end1.y;
								double temp = axial_dist_1 < 0 ? s_i_dir.dot(p1 - s_i_O) : s_i_dir.dot(p2 - s_i_O);
								double x2 = s_i_O.x + s_i_dir[0] * temp, y2 = s_i_O.y + s_i_dir[1] * temp;
								Segment* s_i_crop = new Segment(-1, x1, y1, x2, y2, s_i->width, s_i->support_cost, 
									s_i->grad_normal, false, s_i->is_color_edge);
								s_i_crop->to_delete = true;
								segments.push_back(s_i_crop);
							}
							if (axial_dist_2 < 0 || axial_dist_2 > e_length) { // create cropped version to replace s_i
								double temp = axial_dist_2 < 0 ? s_i_dir.dot(p1 - s_i_O) : s_i_dir.dot(p2 - s_i_O);
								double x1 = s_i_O.x + s_i_dir[0] * temp, y1 = s_i_O.y + s_i_dir[1] * temp;
								double x2 = s_i->end2.x, y2 = s_i->end2.y;
								Segment* s_i_crop = new Segment(-1, x1, y1, x2, y2, s_i->width, s_i->support_cost,
									s_i->grad_normal, false, s_i->is_color_edge);
								s_i_crop->to_delete = true;
								segments.push_back(s_i_crop);
							}
						}
					}
				}
			}
		}
	}
	
}

void Geometry::intersection_polygon_boundary(Ray* r_i, vector<list<Outer_Edge *> > & outer_edges, vector<double> & t_i, vector<Image_Boundary> & r_j, vector<int> & outer_edges_index, vector<double> & t_j, double min_edge_len, double corner_eps)
{
	// We set all intersections time for polygon boundaries to be negative by associating an imaginary ray that propagates from v1 to v2 which reaches v2 at time = 0.
	// We save two things: 1) the first intersection with an outer contour; 2) all intersections with inner contours closer than the outer contour intersection.
	// Assumption: segments that overlap with the boundary have been disabled.

	Point2d & B = r_i->A;
	Vec2d & AB = r_i->OA; // unit vector
	
	t_i = { FLT_MAX };
	t_j = { FLT_MAX };
	r_j.clear(); 
	outer_edges_index.clear();
	r_j.resize(1);
	outer_edges_index.resize(1);

	int intersection = 0;
	for (int i = 0; i < outer_edges.size(); ++i) {
		for (list<Outer_Edge *>::iterator it_h = outer_edges[i].begin(); it_h != outer_edges[i].end(); ++it_h) {
			Point2d & C = (*it_h)->v1_2? (*it_h)->v1->pt : (*it_h)->v2->pt;
			Point2d & D = (*it_h)->v1_2? (*it_h)->v2->pt : (*it_h)->v1->pt;

			Vec2d & CD = (*it_h)->dir_v1v2;  // unit vector

			double l_AB = r_i->initial_length;
			double l_CD = (*it_h)->length;

			if (r_i->collinear_boundaries.count((*it_h)->id_edge) > 0) {
				// Collinear case : intersections exist if the ray propagates towards the boundary
				// For outer boundary edge, we pick the closer intersecting vertex. 
				// For inner boundary edge, the ray either intersects this edge twice, or 0 time.

				double x_C = AB.dot(C - B);
				double x_D = AB.dot(D - B);
				//// We determine x such that C = B + x * AB
				//double x_C = fabs(AB[0]) > fabs(AB[1]) ? (C.x - B.x) / AB[0] : (C.y - B.y) / AB[1];
				//// We determine x such that D = B + x * AB
				//double x_D = fabs(AB[0]) > fabs(AB[1]) ? (D.x - B.x) / AB[0] : (D.y - B.y) / AB[1];
				if (x_C >= -l_AB && x_D >= -l_AB) {
					int contour_id = (*it_h)->init_boundary_location.first;
					if (contour_id == 0) { // record nearest intersection with the outer contour
						double x_i = (x_C < x_D) ? x_C : x_D;
						double y_i = (x_C < x_D) ? -l_CD : 0;
						if (x_i < t_i[0]) {
							t_i[0] = x_i;
							t_j[0] = y_i;
							r_j[0] = POLYGON_OUTER;
							outer_edges_index[0] = i;
							++intersection;
						}
					}
					else { // record all intersections with inner contours closer than the computed outer contour intersection
						if (x_C < t_i[0]) {
							t_i.push_back(x_C);
							t_j.push_back(-l_CD);
							r_j.push_back(POLYGON_INNER);
							outer_edges_index.push_back(i);
						}
						if (x_D < t_i[0]) {
							t_i.push_back(x_D);
							t_j.push_back(0);
							r_j.push_back(POLYGON_INNER);
							outer_edges_index.push_back(i);
						}
						pair<double, double> t_collinear = x_C < x_D ? make_pair(x_C, x_D) : make_pair(x_D, x_C);
						r_i->t_collinear_boundaries.push_back(t_collinear);
					}
					if (fabs(B.x + x_C*AB[0] - C.x) >= 1e-2 || fabs(B.y + x_C*AB[1] - C.y) >= 1e-2) {
						std::cout << "fabs(B.x + x_C*AB[0] - C.x) = " << fabs(B.x + x_C*AB[0] - C.x) << std::endl;
						std::cout << "fabs(B.y + x_C*AB[1] - C.y) = " << fabs(B.y + x_C*AB[1] - C.y) << std::endl;
					}
					if (fabs(B.x + x_D*AB[0] - D.x) >= 1e-2 || fabs(B.y + x_D*AB[1] - D.y) >= 1e-2) {
						std::cout << "fabs(B.x + x_D*AB[0] - D.x) = " << fabs(B.x + x_C*AB[0] - C.x) << std::endl;
						std::cout << "fabs(B.y + x_D*AB[1] - D.y) = " << fabs(B.y + x_D*AB[1] - D.y) << std::endl;
					}
					assert(fabs(B.x + x_C*AB[0] - C.x) < 1e-2 && fabs(B.y + x_C*AB[1] - C.y) < 1e-2);
					assert(fabs(B.x + x_D*AB[0] - D.x) < 1e-2 && fabs(B.y + x_D*AB[1] - D.y) < 1e-2);

				}
				else if (x_C < -l_AB && x_D < -l_AB) {
					// No intersection exists
				}
				else { // non-overlapping check
					std::cout << "x_C, x_D = " << x_C << ", " << x_D << std::endl;
					std::cout << "-l_AB = " << -l_AB << std::endl;
					assert(false);
				}
			}
			else { // Non-collinear case : intersection exists if non-parallel and the ray is not collinear with a neighbor edge
				
				// The following check should only be used for non-collinear cases
				if (AB[1] > 0) {
					if (C.y < r_i->O.y - 0.11 && D.y < r_i->O.y - 0.11) {
						continue;
					}
				}
				else if (AB[1] < 0) {
					if (C.y > r_i->O.y + 0.11 && D.y > r_i->O.y + 0.11) {
						continue;
					}
				}
				if (AB[0] > 0) {
					if (C.x < r_i->O.x - 0.11 && D.x < r_i->O.x - 0.11) {
						continue;
					}
				}
				else if (AB[0] < 0) {
					if (C.x > r_i->O.x + 0.11 && D.x > r_i->O.x + 0.11) {
						continue;
					}
				}

				Vertex* v1 = (*it_h)->v1, *v2 = (*it_h)->v2;
				bool collinear_with_neighbor_1 = false, collinear_with_neighbor_2 = false;
				for (int path = 0; path < v1->directions.size(); ++path) {
					Outer_Edge* e_neighbor = static_cast<Outer_Edge*>(v1->directions[path].second->e);
					if (r_i->collinear_boundaries.count(e_neighbor->id_edge) > 0 && e_neighbor->init_boundary_location.first == (*it_h)->init_boundary_location.first) {
						collinear_with_neighbor_1 = true;
						break;
					}
				}
				for (int path = 0; path < v2->directions.size(); ++path) {
					Outer_Edge* e_neighbor = static_cast<Outer_Edge*>(v2->directions[path].second->e);
					if (r_i->collinear_boundaries.count(e_neighbor->id_edge) > 0 && e_neighbor->init_boundary_location.first == (*it_h)->init_boundary_location.first) {
						collinear_with_neighbor_2 = true;
						break;
					}
				}
				if (collinear_with_neighbor_1 && collinear_with_neighbor_2) {
					continue;
				}

				double determinant = AB[0] * CD[1] - AB[1] * CD[0];
				Vec2d BD = D - B;
				double x_0 = (BD[0] * CD[1] - BD[1] * CD[0]) / determinant;
				double x_i, y_i;
				bool intersect = false;
				   // Let B(x) = B + x * AB, D(y) = D + y * CD
				   // We aim at finding x_0 >= -l_AB and 0 > y >= -l_CD such that B(x_0) = D(y_0)
				   // First heck distance from edge ends to the ray.
				double x_D = AB.dot(BD);
				double x_C = AB.dot(C - B);
				double dist_to_v = 1e4;
				if (x_C >= -l_AB && x_D >= -l_AB) {
					Point2d C_proj(B.x + x_C*AB[0], B.y + x_C*AB[1]);
					double dist_to_C = norm(C_proj - C);
					Point2d D_proj(B.x + x_D*AB[0], B.y + x_D*AB[1]);
					double dist_to_D = norm(D_proj - D);
					dist_to_v = MIN(dist_to_C, dist_to_D);
					if (dist_to_C <= corner_eps && dist_to_D <= corner_eps) {
						if (x_C < x_D) {
							x_i = x_C;
							y_i = -l_CD;
						}
						else {
							x_i = x_D;
							y_i = 0;
						}
					}
					else if (dist_to_C < corner_eps) {
						x_i = x_C;
						y_i = -l_CD;
					}
					else if (dist_to_D < corner_eps) {
						x_i = x_D;
						y_i = 0;
					}
				}
				else if (x_C >= -l_AB) {
					Point2d C_proj(B.x + x_C*AB[0], B.y + x_C*AB[1]);
					dist_to_v = fabs(-AB[1] * (C_proj.x - C.x) + AB[0] * (C_proj.y - C.y));
					x_i = x_C;
					y_i = -l_CD;
				}
				else if (x_D >= -l_AB) {
					Point2d D_proj(B.x + x_D*AB[0], B.y + x_D*AB[1]);
					dist_to_v = fabs(-AB[1] * (D_proj.x - D.x) + AB[0] * (D_proj.y - D.y));
					x_i = x_D;
					y_i = 0;
				}
				intersect = (dist_to_v <= corner_eps);

				if (!intersect) { // If not intersecting the corners
					if (fabs(determinant) >= 1e-6) {// If not parallel, directly compute intersection.
						if (x_0 >= -l_AB) {
							// First check the intersection with the edge
							x_i = x_0;
							y_i = -(AB[0] * BD[1] - AB[1] * BD[0]) / determinant;
							if (y_i < 0 && y_i >= -l_CD) {
								intersect = true;
							}
						}
					}
				}

				if (intersect) {
					int contour_id = (*it_h)->init_boundary_location.first;
					if (contour_id == 0) { // record nearest intersection with the outer contour
						if (x_i < t_i[0]) {
							t_i[0] = x_i;
							t_j[0] = y_i;
							r_j[0] = POLYGON_OUTER;
							outer_edges_index[0] = i;
							++intersection;
						}
					}
					else { // record all intersections with inner contours closer than the computed outer contour intersection
						if (x_i < t_i[0]) {
							t_i.push_back(x_i);
							t_j.push_back(y_i);
							r_j.push_back(POLYGON_INNER);
							outer_edges_index.push_back(i);
						}
					}
				}
			}

		}
	}
	assert(intersection > 0 && r_j[0] == POLYGON_OUTER);
}

void Geometry::intersection_boundary(Ray* r_i, int rows, int cols, double & t_i, Image_Boundary & r_j, double & t_j)
{
	Vec2d dir = r_i->OA;
	t_i = FLT_MAX;

	if (dir[1] > 0) {
		// Collision with the line y = rows-0.5
		double _t_i = (rows - 0.5 - r_i->A.y) / dir[1];
		if (_t_i < t_i) {
			t_i = _t_i;
			t_j = r_i->A.x + t_i * dir[0];
			r_j = TOP_IMAGE;
		}
	} else if (dir[1] < 0) {
		// Collision with the line y = -0.5
		double _t_i = (-0.5 - r_i->A.y) / dir[1];
		if (_t_i < t_i) {
			t_i = _t_i;
			t_j = r_i->A.x + t_i * dir[0];
			r_j = BOTTOM_IMAGE;
		}
	}

	if (dir[0] < 0) {
		// Collision with the line x = -0.5
		double _t_i = (-0.5 - r_i->A.x) / dir[0];
		if (_t_i < t_i) {
			t_i = _t_i;
			t_j = r_i->A.y + t_i * dir[1];
			r_j = LEFT_IMAGE;
		}
	} else if (dir[0] > 0) {
		// Collision with the line x = cols-0.5
		double _t_i = (cols - 0.5 - r_i->A.x) / dir[0];
		if (_t_i < t_i) {
			t_i = _t_i;
			t_j = r_i->A.y + t_i * dir[1];
			r_j = RIGHT_IMAGE;
		}
	}
}


void Geometry::direct_intersection(Ray* r_i, Image_Boundary b_j, double rows, double cols, double & t_i, double & t_j)
{
	switch (b_j) {
	case TOP_IMAGE: // Collision with the line y = rows-0.5
		t_i = (rows - 0.5 - r_i->A.y) / r_i->OA[1];
		t_j = r_i->A.x + t_i * r_i->OA[0];
		return;
	case BOTTOM_IMAGE: // Collision with the line y = -0.5
		t_i = (-0.5 - r_i->A.y) / r_i->OA[1];
		t_j = r_i->A.x + t_i * r_i->OA[0];
		return;
	case LEFT_IMAGE: // Collision with the line x = -0.5
		t_i = (-0.5 - r_i->A.x) / r_i->OA[0];
		t_j = r_i->A.y + t_i * r_i->OA[1];
		return;
	case RIGHT_IMAGE : // Collision with the line x = cols-0.5
		t_i = (cols - 0.5 - r_i->A.x) / r_i->OA[0];
		t_j = r_i->A.y + t_i * r_i->OA[1];
		return;
	}
}


void Geometry::direct_intersection(Ray* r_i, Ray* r_j, double & t_i, double & t_j) 
{
	const Point2d &B = r_i->A, &D = r_j->A;
	const Vec2d &AB = r_i->OA;
	Vec2d DC = -r_j->OA;
	double det = AB[0] * DC[1] - AB[1] * DC[0];
	if (fabs(det) < 1e-6) {
		t_i = t_j = FLT_MAX;
	} else {
		Vec2d BD = D - B;
		t_i = (BD[0] * DC[1] - BD[1] * DC[0]) / det;
		t_j = (AB[0] * BD[1] - AB[1] * BD[0]) / det;
	}
}


bool Geometry::intersection(Ray *s, Ray *t, double max_time_s0, double max_time_s1, 
	unsigned int & s_index, unsigned int & t_index, double & s_time, double & t_time) {
	const Point2d & B = s->A;
	const Point2d & D = t->A;
	const Vec2d & AB = s->OA;
	Vec2d DC = -t->OA;
	// Let B(x) = B + x * AB, D(y) = D + y * CD
	// We aim at finding x_0 >= -l_AB and y_0 >= -l_CD such that B(x_0) = C(y_0)
	double determinant = AB[0] * DC[1] - AB[1] * DC[0];
	// If segments are parallel, return false
	if (fabs(determinant) < 1e-6) {
		return false;
	} else {
		Vec2d BD = D - B;
		double l_AB = s->initial_length;
		double l_CD = t->initial_length;
		//double l_AB = sqrt((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y));
		//double l_CD = sqrt((D.x - C.x) * (D.x - C.x) + (D.y - C.y) * (D.y - C.y));
		double x_i = (BD[0] * DC[1] - BD[1] * DC[0]) / determinant;
		double y_i = (AB[0] * BD[1] - AB[1] * BD[0]) / determinant;
		if (x_i >= -l_AB) {
			s_index = s->index;
			s_time = x_i;
			if (s_time > max_time_s0) return false;
			if (s->t_collinear_boundaries.size() > 0) {
				for (list<pair<double, double>>::iterator it = s->t_collinear_boundaries.begin(); it != s->t_collinear_boundaries.end(); ++it) {
					double t1 = it->first, t2 = it->second;
					if (s_time >= t1 && s_time <= t2) return false;
				}
			}
		} else {
			Ray* s_opp = s->opposite();
			if (s_opp == NULL) return false;
			s_index = s->index + 1;
			if (s_opp->O == s->O) {
				s_time = -2 * l_AB - x_i;
			}
			else {
				s_time = -s->parent->length - x_i;
				if (s_time < -l_AB) return false;
			}

			if (s_time > max_time_s1) return false;
			if (s_opp->t_collinear_boundaries.size() > 0) {
				for (list<pair<double, double>>::iterator it = s_opp->t_collinear_boundaries.begin(); it != s_opp->t_collinear_boundaries.end(); ++it) {
					double t1 = it->first, t2 = it->second;
					if (s_time >= t1 && s_time <= t2) return false;
				}
			}
		}
		if (y_i >= -l_CD) {
			t_index = t->index;
			t_time = y_i;
			if (t->t_collinear_boundaries.size() > 0) {
				for (list<pair<double, double>>::iterator it = t->t_collinear_boundaries.begin(); it != t->t_collinear_boundaries.end(); ++it) {
					double t1 = it->first, t2 = it->second;
					if (t_time >= t1 && t_time <= t2) return false;
				}
			}
		} else {
			Ray* t_opp = t->opposite();
			if (t_opp == NULL) return false;
			t_index = t->index + 1;
			if (t_opp->O == t->O) {
				t_time = -2 * l_CD - y_i;
			}
			else {
				t_time = -t->parent->length - y_i;
				if (t_time < -l_CD) return false;
			}

			if (t_opp->t_collinear_boundaries.size() > 0) {
				for (list<pair<double, double>>::iterator it = t_opp->t_collinear_boundaries.begin(); it != t_opp->t_collinear_boundaries.end(); ++it) {
					double t1 = it->first, t2 = it->second;
					if (t_time >= t1 && t_time <= t2) return false;
				}
			}
		}
		return true;
	}
}


bool Geometry::intersection_colinear(Ray* s, Ray* t, double & s_time, double & t_time)
{
	Point2d & B = s->A;
	Vec2d & AB = s->OA;
	Point2d & D = t->A;
	Vec2d & CD = t->OA;
	// As s and t are colinear, either AB = CD or AB = -CD
	// If AB = CD, the equation has no solution
	assert(fabs(AB[0] * CD[1] - AB[1] * CD[0]) < 1e-5);
	if (fabs(AB[0] - CD[0]) < 1e-5 && fabs(AB[1] - CD[1]) < 1e-5) {
		return false;
	} else {
		Point2d I = (B + D) / 2;
		// We determine s_0 such that I = B + s_0 * AB
		double s_0 = (fabs(AB[0]) > 1e-6 ? (I.x - B.x) / AB[0] : (I.y - B.y) / AB[1]);
		// We determine t_0 such that I = D + t_0 * CD
		double t_0 = (fabs(CD[0]) > 1e-6 ? (I.x - D.x) / CD[0] : (I.y - D.y) / CD[1]);
		// We assert s_0 and t_0 are valid values
		if (s_0 > -s->initial_length && t_0 > -t->initial_length) {
			s_time = s_0;
			t_time = t_0;
			return true;
		} else {
			return false;
		}
	}
}


bool Geometry::are_overlapping(Segment *s, Segment *t, Point2d & P, Point2d & Q)
{
	// If s and t are not parallel, they can't be overlapped
    Vec2d & AB = s->finalDirection;
    Vec2d & CD = t->finalDirection;
	if (fabs(AB[0] * CD[1] - AB[1] * CD[0]) > 1e-5) {
		return false;
	} else {
		// We determine if C or D belongs to [AB]
        Point2d & A = s->finalEnd1;
        Point2d & B = s->finalEnd2;
        Point2d & C = t->finalEnd1;
        Point2d & D = t->finalEnd2;
		double l_AB = s->length;
		double t_C, t_D;
		if (fabs(AB[0]) < 1e-6) {
			if (fabs(C.x - A.x) < 1e-2) {
				t_C = (C.y - A.y) / AB[1];
				t_D = (D.y - A.y) / AB[1];
			} else {
				return false;
			}
		} else if (fabs(AB[1]) < 1e-6) {
			if (fabs(C.y - A.y) < 1e-2) {
				t_C = (C.x - A.x) / AB[0];
				t_D = (D.x - A.x) / AB[0];
			} else {
				return false;
			}
		} else {
			double t_Cx = (C.x - A.x) / AB[0];
			double t_Cy = (C.y - A.y) / AB[1];
			double t_Dx = (D.x - A.x) / AB[0];
			double t_Dy = (D.y - A.y) / AB[1];
			if ((fabs(t_Cx - t_Cy) < 1e-2) && (fabs(t_Dx - t_Dy) < 1e-2)) {				
				t_C = t_Cx;
				t_D = t_Dx;
			} else {
				return false;
			}
		}
		if ((t_C < 0 && t_D < 0) || t_C > l_AB && t_D > l_AB) {
			return false;
		} else {
			// Now, we know are [AB] and [CD] are overlapping
			// A, B, C and D are assigned to abscissas on the line [AB]
			// The ends of the merged segment correspond to the argmin and the argmax of such abscissas
			vector<Point2d> points;
			vector<double> abscissas;
			points.push_back(A); points.push_back(B); 
			points.push_back(C); points.push_back(D);
			abscissas.push_back(0); abscissas.push_back(l_AB);
			abscissas.push_back(t_C); abscissas.push_back(t_D);
			double minimum = FLT_MAX, maximum = -FLT_MAX;
			int argmin = -1, argmax = -1;
			for (int i = 0 ; i < 4 ; i++) {
				if (abscissas[i] < minimum) {
					minimum = abscissas[i];
					argmin = i;
				}
				if (abscissas[i] > maximum) {
					maximum = abscissas[i];
					argmax = i;
				}
			}
			P = points[argmin];
			Q = points[argmax];
			return true;
		}
	}
}

/*----------------------------------------------------------------------------*/
/** Build a region of pixels that share the same angle, up to a
tolerance 'prec', starting at point (x,y).
*/
void Geometry::region_grow(Point2i seed, Matrix<double> & I_m, Matrix<double> & I_t, vector<Point2i> & reg, double & grad_angle, 
	vector<set<Point2i, pointcomp>> & bins, double value_thresh, double bin_width, double prec)
{
	assert(reg.empty());
	int m_rows = static_cast<int>(I_t.rows);
	int m_cols = static_cast<int>(I_t.cols);
	// initialize region's gradient angle
	int i = int(I_t.rows) - 1 - seed.y; // rows
	int j = seed.x; // columns
	grad_angle = I_t(i, j);
	// add seed to container
	reg.push_back(seed);

	double sumdx = I_m(i, j)*cos(I_t(i, j)), sumdy = I_m(i, j)*sin(I_t(i, j));
	/*region growing*/
	for (int k = 0; k < reg.size(); ++k) 
		for (int xx = reg[k].x - 1; xx <= reg[k].x + 1; ++xx) 
			for (int yy = reg[k].y - 1; yy <= reg[k].y + 1; ++yy) 
				if (xx >= 0 && yy >= 0 && xx < m_cols && yy < m_rows) {
					i = m_rows - 1 - yy;
					j = xx;
					int id_bin = (I_m(i, j) > value_thresh &&  I_t(i, j) != NOTDEF) ? floor((I_m(i, j) - value_thresh) / bin_width) : -1;
					if (id_bin == -1) continue;
					double d_angle = Geometry::angle_diff_abs(grad_angle, I_t(i, j));
					if (d_angle > prec) continue;
					id_bin = id_bin >= bins.size() ? bins.size() - 1 : id_bin;
					Point2i neighbor(xx, yy);
					int erased = bins[id_bin].erase(neighbor);
					if (erased > 0) {
						/* add point */
						reg.push_back(neighbor);
						/* update region's angle */
						sumdx += I_m(i, j)*cos(I_t(i, j));
						sumdy += I_m(i, j)*sin(I_t(i, j));
						grad_angle = atan2(sumdy, sumdx);
					}
				}
}

/*----------------------------------------------------------------------------*/
/** Computes a rectangle that covers a region of points.
*/
void Geometry::region2rect(vector<Point2i> & reg,
	Matrix<double> & I_m, double reg_angle,
	double prec, double p, rect * rec)
{
	// reg_angle is the level line angle of the region, which is perpendicular to the gradient direction
	double x, y, dx, dy, l, w, theta, weight, sum, l_min, l_max, w_min, w_max;

	int m_rows = static_cast<int>(I_m.rows);
	int reg_size = reg.size();
	/* check parameters */
	assert(reg_size >= 1 && rec != NULL);

	/* center of the region:

	It is computed as the weighted sum of the coordinates
	of all the pixels in the region. The norm of the gradient
	is used as the weight of a pixel. The sum is as follows:
	cx = \sum_i G(i).x_i
	cy = \sum_i G(i).y_i
	where G(i) is the norm of the gradient of pixel i
	and x_i,y_i are its coordinates.
	*/
	x = y = sum = 0.0;
	for (int k = 0; k<reg_size; ++k)
	{
		int i = m_rows - 1 - reg[k].y;
		int j = reg[k].x;
		weight = I_m(i, j);
		x += (double)reg[k].x * weight;
		y += (double)reg[k].y * weight;
		sum += weight;
	}
	assert(sum > 0.0);
	x /= sum;
	y /= sum;

	/* theta */
	theta = get_rec_theta(reg, x, y, I_m, reg_angle, prec);

	/* length and width:

	'l' and 'w' are computed as the distance from the center of the
	region to pixel i, projected along the rectangle axis (dx,dy) and
	to the orthogonal axis (-dy,dx), respectively.

	The length of the rectangle goes from l_min to l_max, where l_min
	and l_max are the minimum and maximum values of l in the region.
	Analogously, the width is selected from w_min to w_max, where
	w_min and w_max are the minimum and maximum of w for the pixels
	in the region.
	*/
	dx = cos(theta);
	dy = sin(theta);
	l_min = l_max = w_min = w_max = 0.0;
	for (int k = 0; k<reg_size; k++)
	{
		l = ((double)reg[k].x - x) * dx + ((double)reg[k].y - y) * dy;
		w = -((double)reg[k].x - x) * dy + ((double)reg[k].y - y) * dx;

		if (l > l_max) l_max = l;
		if (l < l_min) l_min = l;
		if (w > w_max) w_max = w;
		if (w < w_min) w_min = w;
	}

	/* store values */
	rec->x1 = x + l_min * dx;
	rec->y1 = y + l_min * dy;
	rec->x2 = x + l_max * dx;
	rec->y2 = y + l_max * dy;
	rec->width = w_max - w_min;
	rec->w_max = w_max;
	rec->w_min = w_min;
	rec->x = x;
	rec->y = y;
	rec->theta = theta;
	rec->dx = dx;
	rec->dy = dy;
	rec->prec = prec;
	rec->p = p;

	/* we impose a minimal width of one pixel

	A sharp horizontal or vertical step would produce a perfectly
	horizontal or vertical region. The width computed would be
	zero. But that corresponds to a one pixels width transition in
	the image.
	*/
	if (rec->width < 1.0) {
		double expand = 1.0 - rec->width;
		rec->w_max += 0.5*expand;
		rec->w_min -= 0.5*expand;
		rec->width = 1.0;
	}
}

double Geometry::get_rec_theta(vector<Point2i> & reg, double x, double y,
	Matrix<double> & I_m, double reg_angle, double prec)
{
	double lambda, theta, weight;
	double Ixx = 0.0;
	double Iyy = 0.0;
	double Ixy = 0.0;
	int k;

	int m_rows = static_cast<int>(I_m.rows);
	int reg_size = reg.size();
	/* compute inertia matrix */
	for (k = 0; k<reg_size; ++k)
	{
		int i = m_rows - 1 - reg[k].y;
		int j = reg[k].x;
		weight = I_m(i,j);
		Ixx += ((double)reg[k].y - y) * ((double)reg[k].y - y) * weight;
		Iyy += ((double)reg[k].x - x) * ((double)reg[k].x - x) * weight;
		Ixy -= ((double)reg[k].x - x) * ((double)reg[k].y - y) * weight;
	}

	/* compute smallest eigenvalue */
	lambda = 0.5 * (Ixx + Iyy - sqrt((Ixx - Iyy)*(Ixx - Iyy) + 4.0*Ixy*Ixy));

	/* compute angle */
	if (fabs(Ixx) < 1e-6 && fabs(Iyy) < 1e-6 && fabs(Ixy) < 1e-6) theta = reg_angle;
	else theta = fabs(Ixx)>fabs(Iyy) ? atan2(lambda - Ixx, Ixy) : atan2(Ixy, lambda - Iyy);

	/* The previous procedure doesn't cares about orientation,
	so it could be wrong by 90, 180 degrees. Here is corrected if necessary. */
	while (angle_diff_signed(reg_angle, theta) > M_1_4_PI) theta += M_1_2_PI;
	while (angle_diff_signed(reg_angle, theta) < -M_1_4_PI) theta -= M_1_2_PI;

	// allow +-5 degree difference between theta and reg_angle
	return angle_diff_abs(reg_angle, theta) < 0.087266 ? theta : reg_angle;
}

/*----------------------------------------------------------------------------*/
/** Refine a rectangle.

If a rectangle is not with the right density of region points,
'reduce_region_radius' is called to try to satisfy this condition.
*/
bool Geometry::refine(vector<Point2i> & reg, int * reg_size, 
	Matrix<double> & I_m, vector<set<Point2i, pointcomp>> & bins, double value_thresh, double bin_width,
	double & reg_angle, double & grad_angle, double prec, double p, rect * rec, Matrix<double> & I_t, double density_th)
{
	/* compute region points density */
	double density = (double)*reg_size /
		(dist(rec->x1, rec->y1, rec->x2, rec->y2) * rec->width);
	/*------ First try: reduce rectangle width ------*/
	if (density < density_th)
		reduce_rect_width(reg, reg_size, I_m, bins, value_thresh, bin_width, 
			reg_angle, grad_angle, prec, p, rec, I_t, density_th);

	/*------ Second try: reduce region radius ------*/
	if (density < density_th)
		return reduce_region_radius(reg, reg_size, I_m, bins, value_thresh, bin_width,
			reg_angle, grad_angle, prec, p, rec, I_t, density_th);

	/* if this point is reached, the density criterion is satisfied */
	return TRUE;
}

bool Geometry::reduce_rect_width(vector<Point2i> & reg, int * reg_size,
	Matrix<double> & I_m, vector<set<Point2i, pointcomp>> & bins, double value_thresh, double bin_width,
	double & reg_angle, double & grad_angle, double prec, double p, rect * rec,
	Matrix<double> & I_t, double density_th)
{
	double density;
	vector<Point2i> new_reg;
	vector<Point2i> removed_from_reg;
	removed_from_reg.reserve(*reg_size);
	rect r;
	double delta = 0.5;
	double delta_2 = delta / 2.0;
	/* check parameters */
	assert(*reg_size > 0 && prec > 0.0 && rec != NULL);

	/* compute region points density */
	density = (double)*reg_size /
		(dist(rec->x1, rec->y1, rec->x2, rec->y2) * rec->width);

	/* if the density criterion is satisfied there is nothing to do */
	if (density >= density_th) return TRUE;

	/* while the density criterion is not satisfied, remove farther pixels */
	/* try to reduce width */
	new_reg = reg;
	int new_reg_size = *reg_size;
	rect_copy(rec, &r);
	for (int n = 0; n<6; n++)
	{
		if ((r.width - delta) >= 0.5)
		{
			if (r.w_max >= fabs(r.w_min)) {
				r.w_max -= delta;
				if (r.w_max < 0) r.w_max = 1e-5;
			}
			else {
				r.w_min += delta;
				if (r.w_min > 0) r.w_min = -1e-5;
			}
			r.width -= delta;
			for (int k = 0; k < new_reg_size; k++) {
				double w = -((double)new_reg[k].x - r.x) * r.dy + ((double)new_reg[k].y - r.y) * r.dx;
				if (w < r.w_min || w > r.w_max)
				{
					removed_from_reg.push_back(new_reg[k]);
					/* remove point from the region by replacing with the last element*/
					new_reg[k].x = new_reg[new_reg_size - 1].x; /* if i==*reg_size-1 copy itself */
					new_reg[k].y = new_reg[new_reg_size - 1].y;
					--new_reg_size;
					--k; /* to avoid skipping one point */
				}
			}
			new_reg.resize(new_reg_size);
			/* reject if the region is too small.
			2 is the minimal region size for 'region2rect' to work. */
			if (new_reg_size < 2) return FALSE;

			/* re-compute region points density */
			double density_new = (double)new_reg_size /
				(dist(r.x1, r.y1, r.x2, r.y2) * r.width);
			if (density_new > density)
			{
				density = density_new;
				/* save region */
				reg = new_reg;
				*reg_size = new_reg_size;
				/* re-compute region angle*/
				double sumdx = 0, sumdy = 0;
				for (int k = 0; k < *reg_size; k++) {
					uint i = uint(int(I_t.rows) - 1 - reg[k].y);
					uint j = uint(reg[k].x);
					sumdx += I_m(i, j)*cos(I_t(i, j));
					sumdy += I_m(i, j)*sin(I_t(i, j));
				}
				grad_angle = atan2(sumdy, sumdx);
				reg_angle = grad_angle + PI / 2; // level line angle
				if (reg_angle > PI) reg_angle -= M_2__PI;
				/* re-compute rectangle */
				region2rect(reg, I_m, reg_angle, prec, p, &r);
				/* save rectangle */
				rect_copy(&r, rec);
				/* add back removed points to bins */
				for (int k = 0; k < removed_from_reg.size(); ++k) {
					uint i = uint(int(I_t.rows) - 1 - removed_from_reg[k].y);
					uint j = uint(removed_from_reg[k].x);
					int id_bin = jclamp( 0, floor((I_m(i, j) - value_thresh) / bin_width), bins.size()-1);
					bins[id_bin].insert(removed_from_reg[k]);
				}
				if (density >= density_th) return TRUE;
			}
		}
	}

	return (density >= density_th);
}

/*----------------------------------------------------------------------------*/
/** Reduce the region size, by elimination the points far from the
starting point, until that leads to rectangle with the right
density of region points or to discard the region if too small.
*/
bool Geometry::reduce_region_radius(vector<Point2i> & reg, int * reg_size, 
	Matrix<double> & I_m, vector<set<Point2i, pointcomp>> & bins, double value_thresh, double bin_width,
	double & reg_angle, double & grad_angle, double prec, double p, rect * rec, Matrix<double> & I_t, double density_th)
{
	double density, rad1, rad2, rad, xc, yc;
	/* check parameters */
	assert(*reg_size > 0 &&  prec > 0.0 && rec != NULL);

	/* compute region points density */
	density = (double)*reg_size /
		(dist(rec->x1, rec->y1, rec->x2, rec->y2) * rec->width);

	/* if the density criterion is satisfied there is nothing to do */
	if (density >= density_th) return TRUE;

	/* compute region's radius */
	xc = (double)reg[0].x;
	yc = (double)reg[0].y;
	rad1 = dist(xc, yc, rec->x1, rec->y1);
	rad2 = dist(xc, yc, rec->x2, rec->y2);
	rad = rad1 > rad2 ? rad1 : rad2;

	/* while the density criterion is not satisfied, remove farther pixels */
	while (density < density_th)
	{
		rad *= 0.75; /* reduce region's radius to 75% of its value */

					 /* remove points from the region and update 'used' map */
		for (int k = 0; k<*reg_size; k++)
			if (dist(xc, yc, (double)reg[k].x, (double)reg[k].y) > rad)
			{
				uint i = uint(int(I_m.rows) - 1 - reg[k].y);
				uint j = uint(reg[k].x);
				int id_bin = jclamp(0, floor((I_m(i, j) - value_thresh) / bin_width), bins.size() - 1);
				bins[id_bin].insert(reg[k]);
				/* remove point from the region by replacing with the last element*/
				reg[k].x = reg[*reg_size - 1].x; /* if i==*reg_size-1 copy itself */
				reg[k].y = reg[*reg_size - 1].y;
				--(*reg_size);
				--k; /* to avoid skipping one point */
			}
		reg.resize(*reg_size);
		/* reject if the region is too small.
		2 is the minimal region size for 'region2rect' to work. */
		if (*reg_size < 2) return FALSE;
		/* re-compute region angle*/
		double sumdx = 0, sumdy = 0;
		for (int k = 0; k < *reg_size; k++) {
			uint i = uint(int(I_t.rows) - 1 - reg[k].y);
			uint j = uint(reg[k].x);
			sumdx += I_m(i, j)*cos(I_t(i, j));
			sumdy += I_m(i, j)*sin(I_t(i, j));
		}
		grad_angle = atan2(sumdy, sumdx);
		reg_angle = grad_angle + PI / 2; // level line angle
		if (reg_angle > PI) reg_angle -= M_2__PI;
		/* re-compute rectangle */
		region2rect(reg, I_m, reg_angle, prec, p, rec);

		/* re-compute region points density */
		density = (double)*reg_size /
			(dist(rec->x1, rec->y1, rec->x2, rec->y2) * rec->width);
	}

	/* if this point is reached, the density criterion is satisfied */
	return TRUE;
}

/*----------------------------------------------------------------------------*/
/** Signed angle difference from b to a in [-PI, PI]
*/
double Geometry::angle_diff_signed(double & a, double & b) {
	double diff = a - b;
	while (diff <= -M_PI) diff += M_2__PI;
	while (diff >   M_PI) diff -= M_2__PI;
	return diff;
}

/*----------------------------------------------------------------------------*/
/** Absolute value angle difference in [0, PI]
*/
double Geometry::angle_diff_abs(double & a, double & b) {
	double diff = a - b;
	while (diff <= -M_PI) diff += M_2__PI;
	while (diff >   M_PI) diff -= M_2__PI;
	if (diff < 0.0) diff = -diff;
	return diff;
}

double Geometry::dist(double x1, double y1, double x2, double y2) {
	return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
}

void Geometry::rect_copy(struct rect * in, struct rect * out)
{
	/* check parameters */
	assert(in != NULL || out != NULL);

	/* copy values */
	out->x1 = in->x1;
	out->y1 = in->y1;
	out->x2 = in->x2;
	out->y2 = in->y2;
	out->width = in->width;
	out->w_max = in->w_max;
	out->w_min = in->w_min;
	out->x = in->x;
	out->y = in->y;
	out->theta = in->theta;
	out->dx = in->dx;
	out->dy = in->dy;
	out->prec = in->prec;
	out->p = in->p;
}