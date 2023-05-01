#pragma once

#ifndef REGULARITY_ANGLE_LINE_H
#define REGULARITY_ANGLE_LINE_H
#include "segment_ray.h"

class R_Angle
{
public:
    R_Angle(uint _index, Segment* _segment, double _angle_observed, double _weight);

    ~R_Angle();

    void smooth(double _angle_smoothed);

    uint index;
    Segment* segment;
    double angle_observed;
    double weight;

    double amplitude;
    double angle_smoothed;

    bool smoothed;
    bool queued;
};

#endif // REGULARITY_ANGLE_LINE_H
