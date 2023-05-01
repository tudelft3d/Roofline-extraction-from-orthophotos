#include "r_ordinate.h"


R_Ordinate::R_Ordinate(uint _index, Segment* _segment, double _y_observed, double _weight)
{
    index = _index;
    segment = _segment;
    y_observed = _y_observed;
    weight = _weight;

    amplitude = FLT_MAX;
    y_smoothed = FLT_MAX;
    smoothed = false;
    queued = false;
}


R_Ordinate::~R_Ordinate()
{

}


void R_Ordinate::smooth(double _y_smoothed)
{
    y_smoothed = _y_smoothed;
    smoothed = true;
}
