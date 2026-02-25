/**
 * ContourGeometry.h
 *
 * !! 此文件仅用于独立编译测试 !!
 * !! 用户项目中请替换为自己的几何类型头文件 !!
 *
 * 作者：ContourDev
 * 创建时间：2026-02-18
 */

#ifndef CONTOURGEOMETRY_H
#define CONTOURGEOMETRY_H

#include <vector>
#include <cmath>

#ifndef EPSILON
#define EPSILON 1e-8
#endif

class MPoint
{
public:
    float X;
    float Y;

    MPoint() : X(0), Y(0) {}
    MPoint(const float x, const float y) : X(x), Y(y) {}
    MPoint(const MPoint& pt) : X(pt.X), Y(pt.Y) {}

    bool operator==(const MPoint& pt) const
    {
        return fabs(pt.X - X) < EPSILON && fabs(pt.Y - Y) < EPSILON;
    }

    bool operator!=(const MPoint& pt) const
    {
        return fabs(pt.X - X) > EPSILON || fabs(pt.Y - Y) > EPSILON;
    }

    MPoint& operator=(const MPoint& pt)
    {
        if (this == &pt) return *this;
        X = pt.X; Y = pt.Y;
        return *this;
    }
};

class MCurve
{
public:
    std::vector<MPoint> Points;
    MCurve() {}
    virtual ~MCurve() {}
};

class MRing : public MCurve
{
public:
    MRing() {}
};

class MPolygon
{
public:
    std::vector<MRing> OuterRings;
    std::vector<MRing> InnerRings;

    void Clear() { OuterRings.clear(); InnerRings.clear(); }
};

struct MContourLine
{
    MCurve Curve;
    float Level;
};

struct MContourFill
{
    MPolygon Polygon;
    float LevelMax, LevelMin;
    int ColorIndex;
};

#endif // CONTOURGEOMETRY_H
