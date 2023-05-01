#include "refine_partition.h"
#include "defs.h"
#include "trace.h"
#include <ctime>
#include "svg.h"

void Refine_Partition::preprocess_graph(Kinetic_Model *model)
{
	std::cout << "Preprocessing ...... ";

	model->graph->estimate_cdf_grad_m(model->I_grad_m, true);
	model->graph->estimate_cdf_grad_m(model->I_lab_grad_m, false);
	model->graph->weight_all_edges(model->used, model->I_grad_m, model->I_grad_t);
	model->graph->find_pixels_inside_all_facets();
	model->graph->find_active_gradients_all_facets(model->I_grad_m, model->I_grad_t, true); // always call after weighting edges
	model->graph->find_active_gradients_all_facets(model->I_lab_grad_m, model->I_lab_grad_t, false);
	model->graph->compute_feature_all_facets(model->I_prob, model->I);
	model->graph->naive_labeling();
#if NOT_MEASURING_PERFORMANCES
	delete model->initial_graph;
	model->initial_graph = new Partition(*(model->graph));
#endif
	model->graph->detect_line_segments_all_facets(model->I, model->I_grad_m, model->I_grad_t, true);
	model->graph->detect_line_segments_all_facets(model->I, model->I_lab_grad_m, model->I_lab_grad_t, false);
	model->graph->extend_detected_segments_all_facets(model->I, model->I_prob, model->used, model->I_grad_m, model->I_grad_t, model->reg_angles, model->bin_to_angle_index);
	model->graph->build_dual_graph(model->I, model->I_grad_m, model->I_grad_t);
	model->graph->compute_energy();
	model->graph->mark_edges_next_operation();

	make_layer(model);
	std::cout << "DONE. " << std::endl;
	std::cout << "[#Vertices, #Edges] = " << "[" << model->graph->get_id_vertices() << ", " << model->graph->id_edges << "]" << std::endl;
}

void Refine_Partition::grouping(Kinetic_Model *model)
{
	clock_t t_begin, t_end;
	t_begin = clock();

	if (model->graph->get_n_iterations() == 0) {
		preprocess_graph(model);
	}

	merge_and_split(model, model->params->max_iters);
	model->graph->mark_edges_next_operation();

	make_layer(model);

	t_end = clock();
	trace(model->params->verbose_level, 5, "** Output " + std::to_string(model->graph->faces.size()) + " facets in " + std::to_string(float(t_end - t_begin) / CLOCKS_PER_SEC) + " s.");
}

void Refine_Partition::merge_and_split(Kinetic_Model *model, int max_iters)
{	
	Partition* P = model->graph;
	std::string save_location = model->params->path_output_directory + "\\" + model->basename + "_iterations";

	for (int i = 0; i < max_iters; ++i) {
		trace(model->params->verbose_level, 5, "i = " + std::to_string(i));

		// call next operation
		Operator_Type op = P->realize_next_operation(model->I, model->I_prob, model->used, model->I_grad_m, model->I_grad_t, model->I_lab_grad_m, model->I_lab_grad_t, model->reg_angles, model->bin_to_angle_index);
		if (op == STOP) break;
	}	
}


void Refine_Partition::make_layer(Kinetic_Model *model)
{
	Point2d shift(0.5, 0.5);

	Matrix<uchar> & I = model->I;

	uchar color_0[3] = { 100, 240, 255 };
	uchar color_1[3] = { 0, 240, 240 };
	uchar color_2[3] = { 35, 120, 180 };
	uchar color_3[3] = { 105, 20, 105 };

	model->clear_line_items(model->L_prop);
	Edge* e = model->graph->edges_head;
	while (e != NULL) {
		Point2d pt_1 = e->v1->pt + shift;
		Point2d pt_2 = e->v2->pt + shift;
		Point2d pc_1 = Point2d(jclamp(0, pt_1.x, I.cols - 1), jclamp(0, I.rows - pt_1.y, I.rows - 1));
		Point2d pc_2 = Point2d(jclamp(0, pt_2.x, I.cols - 1), jclamp(0, I.rows - pt_2.y, I.rows - 1));
		model->add_line(model->L_prop, pc_1.x, pc_1.y, pc_2.x, pc_2.y, 255, 0, 0);
		e = e->e_next;
	}

	model->clear_line_items(model->L_next_operation);
	for (int i = 0; i < model->graph->edges_next_operation.size(); ++i) {
		const auto& next_op_edges = model->graph->edges_next_operation[i];
		for (const auto& point_pair : next_op_edges) {
			Point2d pt_1 = point_pair.first + shift;
			Point2d pt_2 = point_pair.second + shift;
			Point2d pc_1 = Point2d(jclamp(0, pt_1.x, I.cols - 1), jclamp(0, I.rows - pt_1.y, I.rows - 1));
			Point2d pc_2 = Point2d(jclamp(0, pt_2.x, I.cols - 1), jclamp(0, I.rows - pt_2.y, I.rows - 1));
			switch (i) {
				case 0: 
					model->add_line(model->L_next_operation, pc_1.x, pc_1.y, pc_2.x, pc_2.y, 0, 255, 0);
					break;
				case 1:
					model->add_line(model->L_next_operation, pc_1.x, pc_1.y, pc_2.x, pc_2.y, 0, 0, 255);
					break;
				case 2: 
					model->add_line(model->L_next_operation, pc_1.x, pc_1.y, pc_2.x, pc_2.y, 0, 255, 255);
					break;
			}
		}
	}

	model->clear_line_items(model->L_split_segments);
	for (Face* f: model->graph->faces) {
		int n_segs = 0;
		for (auto it = f->detected_segments.begin(); it != f->detected_segments.end(); ++it) {
			Point2d pt_1 = (*it)->end1 + shift;
			Point2d pt_2 = (*it)->end2 + shift;
			Point2d pc_1 = Point2d(jclamp(0, pt_1.x, I.cols - 1), jclamp(0, I.rows - pt_1.y, I.rows - 1));
			Point2d pc_2 = Point2d(jclamp(0, pt_2.x, I.cols - 1), jclamp(0, I.rows - pt_2.y, I.rows - 1));
			model->add_line(model->L_split_segments, pc_1.x, pc_1.y, pc_2.x, pc_2.y, 105, 20, 105);
			++n_segs;
			if (n_segs >= 10) break;
		}
	}
}
