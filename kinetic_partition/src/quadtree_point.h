#ifndef QUADTREE_POINT_H
#define QUADTREE_POINT_H

#include <opencv2/core.hpp>

using cv::Point2d;


class Quadtree_Point
{
public:
    Quadtree_Point();

    virtual ~Quadtree_Point();

    virtual Point2d keypoint() = 0;
};

#endif // QUADTREE_POINT_H
