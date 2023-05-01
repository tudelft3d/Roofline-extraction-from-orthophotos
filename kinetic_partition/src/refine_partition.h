#pragma once
#include "kinetic_model.h"

namespace Refine_Partition {

	void preprocess_graph(Kinetic_Model *model);

	void grouping(Kinetic_Model *model);

	void merge_and_split(Kinetic_Model *model, int max_iters);

	void make_layer(Kinetic_Model *model);

	void svg_animation(Kinetic_Model *model, std::string & directory, int stepsize);

}