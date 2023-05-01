#include "r_angle.h"


R_Angle::R_Angle(uint _index, Segment *_segment, double _angle_observed, double _weight)
{
    index = _index;
    segment = _segment;
    angle_observed = _angle_observed;
    weight = _weight;

    amplitude = FLT_MAX;
    angle_smoothed = FLT_MAX;

    smoothed = false;
    queued = false;
}


R_Angle::~R_Angle()
{

}


void R_Angle::smooth(double _angle_smoothed)
{
    smoothed = true;
    angle_smoothed = _angle_smoothed;
}
