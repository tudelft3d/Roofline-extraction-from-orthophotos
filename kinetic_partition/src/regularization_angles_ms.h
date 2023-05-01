#ifndef REGULARIZATION_ANGLES_RTREE_H
#define REGULARIZATION_ANGLES_RTREE_H

#include "defs.h"
#include "kinetic_model.h"
#include "regularization_angles.h"
#include "r_angle.h"

#include <vector>
#include <random>

using std::vector;
using std::default_random_engine;
using std::uniform_int_distribution;

class Regularization_Angles_Mean_Shift : public Regularization_Angles
{
public:
    Regularization_Angles_Mean_Shift();

    virtual ~Regularization_Angles_Mean_Shift();

    void regularize(Kinetic_Model* model);


    void build_rtree(vector<R_Angle *> & regularities, Boost_RTree & rtree_regularities);

    void search_neighborhood(Segment* s, double D, Boost_RTree & rtree_regularities, vector<R_Angle *> & regularities, list<R_Angle *> & neighborhood);

    void search_neighborhood(R_Angle* r, double D, Boost_RTree & rtree_regularities, vector<R_Angle *> & regularities, list<R_Angle *> & neighborhood);


    void identify_regularities(Kinetic_Model* model, Boost_RTree & rtree_segments, vector<R_Angle *> & regularities, bool include_artificial);

    bool propose_regularity_by_mean_shift(Kinetic_Model* model, Boost_RTree & rtree_segments, Segment* seed, double & value, bool include_artificial);

    void smooth_regularities(Kinetic_Model* model, Boost_RTree & rtree_regularities, vector<R_Angle *> & all_regularities);

    bool smooth_by_mean_shift(Kinetic_Model* model, R_Angle* seed, Boost_RTree & rtree_regularities, vector<R_Angle *> &all_regularities, int min_terms, double & z_0);

    void assign_segments_to_local_regularities(Kinetic_Model* model, Boost_RTree & rtree_regularities, vector<R_Angle *> & all_regularities);

private:
    bool fabs_mod_90(double a, double b, double eps) {
        for (int k = -1; k <= 1; k++) {
            if (fabs(a + 90 * k - b) < eps) return true;
        }
        return false;
    }

    bool inside_interval_mod_90(double a, double b, double eps) {
        for (int k = -1; k <= 1; k++) {
            if (a - eps <= b + 90 * k && b + 90 * k <= a + eps) return true;
        }
        return false;
    }

    bool fabs_mod_180(double a, double b, double eps) {
        for (int k = -1; k <= 1; k++) {
            if (fabs(a + 180 * k - b) < eps) return true;
        }
        return false;
    }

    bool inside_interval_mod_180(double a, double b, double eps) {
        for (int k = -1; k <= 1; k++) {
            if (a - eps <= b + 180 * k && b + 180 * k <= a + eps) return true;
        }
        return false;
    }

    double kernel(double k, double x) {
        return exp(-k * x * x);
    }
};

#endif // REGULARIZATION_ANGLES_RTREE_H
