#ifndef R_ORDINATE_H
#define R_ORDINATE_H

#include "segment_ray.h"

class R_Ordinate
{
public:
    R_Ordinate(uint _index, Segment* _segment, double _y_observed, double _weight);

    ~R_Ordinate();

    void smooth(double _y_smoothed);

    uint index;
    Segment* segment;
    double y_observed;
    double weight;

    double amplitude;
    double y_smoothed;

    bool smoothed;
    bool queued;
};

#endif // R_ORDINATE_H
