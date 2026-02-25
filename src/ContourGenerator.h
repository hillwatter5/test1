/**
 * ContourGenerator.h
 * 功能说明：等值线与填色多边形的核心计算类。
 *          基于 Marching Squares 算法，复现 matplotlib.contour/contourf 的底层逻辑。
 *          计算结果用 MContourLine 和 MContourFill 存储。
 * 作者：ContourDev
 * 创建时间：2026-02-18
 */

#ifndef CONTOURGENERATOR_H
#define CONTOURGENERATOR_H

// !! 请替换为您项目中提供 MPoint/MCurve/MRing/MPolygon/MContourLine/MContourFill 的头文件 !!
#include "ContourGeometry.h"

#include <vector>
#include <map>
#include <utility>
#include <cfloat>
#include <cmath>

/**
 * 合并计算结果：等值线 + 填色多边形
 */
struct CMergedResult
{
    std::vector<MContourLine> ContourLines;
    std::vector<MContourFill> FilledContours;
};

// ---- 算法内部辅助类型 ----

typedef std::pair<MPoint, MPoint> CSegPair;

struct CPointKey
{
    long long Qx, Qy;
    CPointKey() : Qx(0), Qy(0) {}
    explicit CPointKey(const MPoint& pt)
    {
        Qx = (long long)std::round(pt.X * 10000000.0);
        Qy = (long long)std::round(pt.Y * 10000000.0);
    }
    bool operator<(const CPointKey& o) const { return Qx < o.Qx || (Qx == o.Qx && Qy < o.Qy); }
    bool operator==(const CPointKey& o) const { return Qx == o.Qx && Qy == o.Qy; }
};

struct CEdgeKey
{
    CPointKey Start, End;
    CEdgeKey() {}
    CEdgeKey(const CPointKey& s, const CPointKey& e) : Start(s), End(e) {}
    bool operator<(const CEdgeKey& o) const
    {
        if (Start < o.Start) return true;
        if (o.Start < Start) return false;
        return End < o.End;
    }
};

struct CCachedIsect { bool Has; MPoint Pt; CCachedIsect() : Has(false) {} };

/** Shoelace 有符号面积，正=逆时针（外环），负=顺时针（内环） */
double ComputeSignedArea(const std::vector<MPoint>& pts);

/**
 * 等值线与填色多边形核心计算引擎
 *
 * 用法示例：
 *   CContourGenerator gen;
 *   gen.Initialize(xCoords, cols, yCoords, rows, values);
 *   float levels[] = {1.0f, 2.5f, 3.0f};
 *   CMergedResult result = gen.ComputeContourAndFill(levels, 3);
 *   // result.ContourLines  → vector<MContourLine>
 *   // result.FilledContours → vector<MContourFill>
 */
class CContourGenerator
{
public:
    CContourGenerator();
    ~CContourGenerator();

    /**
     * 初始化网格数据
     * @param xCoords 列坐标数组，长度 cols
     * @param cols    列数
     * @param yCoords 行坐标数组，长度 rows
     * @param rows    行数
     * @param values  网格值，长度 rows*cols，行优先
     */
    void Initialize(const float* xCoords, int cols,
                    const float* yCoords, int rows,
                    const float* values);

    /**
     * 仅计算等值线
     * N 个 level → 每个 level 可能产生多条 MContourLine
     */
    std::vector<MContourLine> ComputeContourLines(
        const float* levels, int numLevels);

    /**
     * 仅计算填色多边形
     * N 个 level → N+1 个区间，每个区间可能多个 MContourFill
     *   ColorIndex=0   : value < levels[0]
     *   ColorIndex=k   : levels[k-1] <= value < levels[k]
     *   ColorIndex=N   : value >= levels[N-1]
     */
    std::vector<MContourFill> ComputeFilledContours(
        const float* levels, int numLevels);

    /**
     * 合并计算（性能优化：单次遍历 + 交点缓存）
     */
    CMergedResult ComputeContourAndFill(
        const float* levels, int numLevels);

private:
    int _rows, _cols;
    std::vector<float> _xCoords, _yCoords, _values;

    float GetValue(int row, int col) const;
    MPoint GetGridPoint(int row, int col) const;

    std::vector<MCurve> ChainSegments(std::vector<CSegPair>& segments) const;

    void CancelOpposingEdges(
        std::vector<CSegPair>& edges,
        std::vector<CSegPair>& remaining) const;

    std::vector<MRing> TraceRings(std::vector<CSegPair>& edges) const;

    std::vector<MPolygon> GroupRingsToPolygons(std::vector<MRing>& rings) const;

    std::vector<MPolygon> AssemblePolygonsFromEdges(
        std::vector<CSegPair>& edges) const;
};

#endif // CONTOURGENERATOR_H
