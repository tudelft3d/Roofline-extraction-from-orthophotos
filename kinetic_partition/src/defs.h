#pragma once

#ifndef Q_MOC_RUN
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#endif

const char separator =
#ifdef _WIN32
    '\\';
#else
    '/';
#endif


#ifndef jclamp
#define jclamp(a, x, b) ((x) < (a) ? (a) : ((x) > (b) ? (b) : (x)))
#endif
#ifndef SIGN
#define SIGN(a) ((a > 0) - (a < 0))
#endif

#ifndef PI
#define PI 3.141592653589783238462643383279
#endif
/** 1/4 pi */
#ifndef M_1_4_PI
#define M_1_4_PI 0.78539816339
#endif
/** 1/2 pi */
#ifndef M_1_2_PI
#define M_1_2_PI 1.57079632679
#endif
/** 3/2 pi */
#ifndef M_3_2_PI
#define M_3_2_PI 4.71238898038
#endif
/** 2 pi */
#ifndef M_2__PI
#define M_2__PI  6.28318530718
#endif
#ifndef NOTDEF
#define NOTDEF -1024.0
#endif
#ifndef INVALID_LABEL
#define INVALID_LABEL 255
#endif
typedef unsigned char uchar;
typedef unsigned int uint;

struct float2 { float x, y; };
struct float3 { float x, y, z; };

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 2, bg::cs::cartesian> Boost_Point;
typedef bg::model::box<Boost_Point> Boost_Box;
typedef std::pair<Boost_Box, uint> Boost_Value;
typedef bgi::rtree<Boost_Value, bgi::quadratic<16> > Boost_RTree;
