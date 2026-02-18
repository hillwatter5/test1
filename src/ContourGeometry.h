/**
 * ContourGeometry.h
 * 功能说明：定义等值线与填色算法所需的基础几何类型，
 *          包括点、矩形、曲线、环和多边形。
 * 作者：ContourDev
 * 创建时间：2026-02-18
 */

#ifndef CONTOURGEOMETRY_H
#define CONTOURGEOMETRY_H

#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>

#ifndef EPSILON
#define EPSILON 1e-8
#endif

namespace slib {
    inline double Distance(float x1, float y1, float x2, float y2)
    {
        double dx = x2 - x1;
        double dy = y2 - y1;
        return std::sqrt(dx * dx + dy * dy);
    }
}

/**
 * 形状类型枚举
 */
enum MG_Shape
{
    MG_Shape_Point = 0,
    MG_Shape_Curve,
    MG_Shape_Ring,
    MG_Shape_Polygon
};

#define MG_BLOCKS_SHAPE_DEFINE_ID() int ShapeID = 0

// ============================================================
//  MPoint - 二维点
// ============================================================
class MPoint
{
public:
    float X;
    float Y;

public:
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
        X = pt.X;
        Y = pt.Y;
        return *this;
    }

    bool Contain(const MPoint& pt) const
    {
        return (fabs(X - pt.X) < EPSILON && fabs(Y - pt.Y) < EPSILON);
    }
};

// ============================================================
//  MRect - 矩形包围盒
// ============================================================
class MRect
{
public:
    float MinX;
    float MinY;
    float MaxX;
    float MaxY;

public:
    MRect() : MinX(0), MinY(0), MaxX(0), MaxY(0) {}

    MRect(float minX, float minY, float maxX, float maxY)
        : MinX(minX), MinY(minY), MaxX(maxX), MaxY(maxY) {}

    bool IsEmpty() const
    {
        return (fabs(MaxX - MinX) < EPSILON && fabs(MaxY - MinY) < EPSILON);
    }

    void Clear()
    {
        MinX = MinY = MaxX = MaxY = 0;
    }

    /**
     * 将另一个矩形合并到当前矩形中（取并集包围盒）
     */
    void Join(const MRect& other)
    {
        if (IsEmpty()) {
            MinX = other.MinX;
            MinY = other.MinY;
            MaxX = other.MaxX;
            MaxY = other.MaxY;
        } else {
            if (other.MinX < MinX) MinX = other.MinX;
            if (other.MinY < MinY) MinY = other.MinY;
            if (other.MaxX > MaxX) MaxX = other.MaxX;
            if (other.MaxY > MaxY) MaxY = other.MaxY;
        }
    }

    bool ContainPoint(const MPoint& pt) const
    {
        return pt.X >= MinX && pt.X <= MaxX &&
               pt.Y >= MinY && pt.Y <= MaxY;
    }
};

// ============================================================
//  MCurve - 多段线（开放曲线）
// ============================================================
class MCurve
{
public:
    std::vector<MPoint> Points;

public:
    MCurve() {}

    virtual ~MCurve() {}

    virtual void Clear()
    {
        Points.clear();
    }

    void AddPoint(const MPoint& pt)
    {
        Points.push_back(pt);
    }

    virtual bool Contain(const MPoint& pt) const
    {
        for (size_t i = 0; i < Points.size(); i++) {
            if (Points[i].Contain(pt))
                return true;
        }
        return false;
    }

    virtual void UpdateMBR()
    {
        _MBR.Clear();
        if (Points.empty()) return;
        _MBR.MinX = _MBR.MaxX = Points[0].X;
        _MBR.MinY = _MBR.MaxY = Points[0].Y;
        for (size_t i = 1; i < Points.size(); i++) {
            if (Points[i].X < _MBR.MinX) _MBR.MinX = Points[i].X;
            if (Points[i].X > _MBR.MaxX) _MBR.MaxX = Points[i].X;
            if (Points[i].Y < _MBR.MinY) _MBR.MinY = Points[i].Y;
            if (Points[i].Y > _MBR.MaxY) _MBR.MaxY = Points[i].Y;
        }
    }

    virtual MRect GetMBR()
    {
        if (_MBR.IsEmpty())
            UpdateMBR();
        return _MBR;
    }

    virtual double Perimeter() const
    {
        int n = (int)Points.size();
        if (n < 2) return 0;
        double perimeter = 0;
        for (int i = 0; i < n - 1; i++) {
            perimeter += slib::Distance(Points[i].X, Points[i].Y,
                                        Points[i + 1].X, Points[i + 1].Y);
        }
        return perimeter;
    }

protected:
    MRect _MBR;
};

// ============================================================
//  MRing - 闭合环（多边形边界）
// ============================================================
class MRing : public MCurve
{
public:
    MRing() {}

    /**
     * 射线法判断点是否在环内
     */
    bool Contain(const MPoint& pt) const override
    {
        int n = (int)Points.size();
        if (n < 3) return false;
        bool inside = false;
        for (int i = 0, j = n - 1; i < n; j = i++) {
            if (((Points[i].Y > pt.Y) != (Points[j].Y > pt.Y)) &&
                (pt.X < (Points[j].X - Points[i].X) *
                 (pt.Y - Points[i].Y) / (Points[j].Y - Points[i].Y) +
                 Points[i].X))
            {
                inside = !inside;
            }
        }
        return inside;
    }

    /**
     * 判断环是否包含另一个矩形
     */
    bool Contain(const MRect& rect) const
    {
        MPoint corners[4] = {
            MPoint(rect.MinX, rect.MinY),
            MPoint(rect.MaxX, rect.MinY),
            MPoint(rect.MaxX, rect.MaxY),
            MPoint(rect.MinX, rect.MaxY)
        };
        for (int i = 0; i < 4; i++) {
            if (!Contain(corners[i]))
                return false;
        }
        return true;
    }

    /**
     * 判断环是否包含另一个环
     */
    bool Contain(const MRing* pRing) const
    {
        if (!pRing || pRing->Points.empty()) return false;
        for (size_t i = 0; i < pRing->Points.size(); i++) {
            if (!Contain(pRing->Points[i]))
                return false;
        }
        return true;
    }

    /**
     * 判断与另一个环是否相交（简化：包围盒相交判断）
     */
    bool IsIntersectWith(MRing* pRing)
    {
        MRect mbrA = GetMBR();
        MRect mbrB = pRing->GetMBR();
        return !(mbrA.MaxX < mbrB.MinX || mbrA.MinX > mbrB.MaxX ||
                 mbrA.MaxY < mbrB.MinY || mbrA.MinY > mbrB.MaxY);
    }

    /**
     * 判断与矩形是否相交
     */
    bool IsIntersectWith(MRect* pRect)
    {
        MRect mbr = GetMBR();
        return !(mbr.MaxX < pRect->MinX || mbr.MinX > pRect->MaxX ||
                 mbr.MaxY < pRect->MinY || mbr.MinY > pRect->MaxY);
    }

    /**
     * 计算有符号面积（Shoelace 公式）
     * 正值 = 逆时针（外环），负值 = 顺时针（内环/孔洞）
     */
    double SignedArea() const
    {
        int n = (int)Points.size();
        if (n < 3) return 0;
        double area = 0;
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            area += (double)Points[i].X * Points[j].Y;
            area -= (double)Points[j].X * Points[i].Y;
        }
        return area / 2.0;
    }

    double Perimeter() const override
    {
        int n = (int)Points.size();
        if (n < 2) return 0;
        double perimeter = 0;
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            perimeter += slib::Distance(Points[i].X, Points[i].Y,
                                        Points[j].X, Points[j].Y);
        }
        return perimeter;
    }
};

// ============================================================
//  MPolygon - 多边形（含外环与内孔）
// ============================================================
class MPolygon
{
public:
    std::vector<MRing> OuterRings;
    std::vector<MRing> InnerRings;

public:
    MPolygon() {}

    MG_Shape GetType() const { return MG_Shape_Polygon; }

    MRect GetMBR()
    {
        if (_MBR.IsEmpty())
            UpdateMBR();
        return _MBR;
    }

    void UpdateMBR()
    {
        _MBR.Clear();
        int nRingNum = (int)OuterRings.size();
        if (nRingNum > 0) {
            for (int i = 0; i < nRingNum; i++) {
                OuterRings[i].UpdateMBR();
                _MBR.Join(OuterRings[i].GetMBR());
            }
        } else {
            nRingNum = (int)InnerRings.size();
            for (int i = 0; i < nRingNum; i++) {
                InnerRings[i].UpdateMBR();
                _MBR.Join(InnerRings[i].GetMBR());
            }
        }
    }

    int GetTotalPointCount() const
    {
        int count = 0;
        for (size_t i = 0; i < OuterRings.size(); i++)
            count += (int)OuterRings[i].Points.size();
        for (size_t i = 0; i < InnerRings.size(); i++)
            count += (int)InnerRings[i].Points.size();
        return (int)(count + OuterRings.size() + InnerRings.size());
    }

    bool Contain(const MPoint& pt) const
    {
        for (size_t i = 0; i < OuterRings.size(); i++) {
            if (OuterRings[i].Contain(pt))
                return true;
        }
        return false;
    }

    void Clear()
    {
        OuterRings.clear();
        InnerRings.clear();
        _MBR.Clear();
    }

private:
    MRect _MBR;
    MG_BLOCKS_SHAPE_DEFINE_ID();
};

#endif // CONTOURGEOMETRY_H
