/**
 * ContourGenerator.h
 * 功能说明：等值线与填色多边形的核心计算类声明。
 *          基于 Marching Squares 算法，复现 matplotlib.contour/contourf 的底层逻辑。
 * 作者：ContourDev
 * 创建时间：2026-02-18
 */

#ifndef CONTOURGENERATOR_H
#define CONTOURGENERATOR_H

#include "ContourGeometry.h"
#include <vector>
#include <map>
#include <utility>
#include <cfloat>

/**
 * 等值线计算结果：一个等值级别对应的所有曲线
 */
struct CContourLineResult
{
    float Level;
    std::vector<MCurve> Curves;
};

/**
 * 填色计算结果：一个等值区间对应的所有多边形
 *
 * 对 N 个 levels，产生 N+1 个区间：
 *   ColorIndex=0 : value < levels[0]
 *   ColorIndex=k : levels[k-1] <= value < levels[k]  (1 <= k < N)
 *   ColorIndex=N : value >= levels[N-1]
 */
struct CFilledContourResult
{
    float LevelLow;
    float LevelHigh;
    int ColorIndex;
    std::vector<MPolygon> Polygons;

    CFilledContourResult() : LevelLow(0), LevelHigh(0), ColorIndex(-1) {}
};

/**
 * 合并计算结果：等值线 + 填色多边形
 */
struct CMergedResult
{
    /** N 个 level 对应 N 条等值线结果 */
    std::vector<CContourLineResult> ContourLines;

    /** N 个 level 对应 N+1 个填色区间结果，每个含 ColorIndex */
    std::vector<CFilledContourResult> FilledContours;
};

/**
 * 量化后的点坐标，用于浮点比较时做哈希键
 */
struct CPointKey
{
    long long Qx;
    long long Qy;

    CPointKey() : Qx(0), Qy(0) {}

    explicit CPointKey(const MPoint& pt)
    {
        Qx = (long long)std::round(pt.X * 10000000.0);
        Qy = (long long)std::round(pt.Y * 10000000.0);
    }

    bool operator<(const CPointKey& other) const
    {
        return (Qx < other.Qx) || (Qx == other.Qx && Qy < other.Qy);
    }

    bool operator==(const CPointKey& other) const
    {
        return Qx == other.Qx && Qy == other.Qy;
    }
};

/**
 * 有向边的量化键，用于边对消
 */
struct CEdgeKey
{
    CPointKey Start;
    CPointKey End;

    CEdgeKey() {}
    CEdgeKey(const CPointKey& s, const CPointKey& e) : Start(s), End(e) {}

    bool operator<(const CEdgeKey& other) const
    {
        if (Start < other.Start) return true;
        if (other.Start < Start) return false;
        return End < other.End;
    }

    bool operator==(const CEdgeKey& other) const
    {
        return Start == other.Start && End == other.End;
    }
};

/**
 * 边交点缓存条目：一条网格边在某个 level 上的交点
 */
struct CCachedIsect
{
    bool Has;
    MPoint Pt;

    CCachedIsect() : Has(false) {}
};

typedef std::pair<MPoint, MPoint> CSegPair;

/**
 * 等值线和填色区间的核心计算引擎
 *
 * 基于 Marching Squares 算法：
 * - 等值线：扫描所有网格单元，用线性插值计算等值线与网格边的交点，
 *           生成线段后链接为连续多段线（MCurve）。
 * - 填色：对每个单元格生成等值区间内的有向边，然后对相邻单元共享边对消，
 *         追踪剩余边组装为闭合环（MRing），最终按拓扑归组为 MPolygon。
 */
class CContourGenerator
{
public:
    CContourGenerator();
    ~CContourGenerator();

    /**
     * 初始化网格数据
     * @param xCoords 列坐标数组，长度为 cols
     * @param cols    列数（x 方向网格点数）
     * @param yCoords 行坐标数组，长度为 rows
     * @param rows    行数（y 方向网格点数）
     * @param values  网格值数组，长度为 rows*cols，行优先存储
     */
    void Initialize(const float* xCoords, int cols,
                    const float* yCoords, int rows,
                    const float* values);

    /**
     * 仅计算等值线（N 个 level → N 组曲线）
     */
    std::vector<CContourLineResult> ComputeContourLines(
        const float* levels, int numLevels);

    /**
     * 仅计算填色多边形（N 个 level → N+1 个区间）
     * 区间规则：
     *   [0] value < levels[0]        (ColorIndex=0)
     *   [k] levels[k-1] <= value < levels[k]  (ColorIndex=k)
     *   [N] value >= levels[N-1]     (ColorIndex=N)
     */
    std::vector<CFilledContourResult> ComputeFilledContours(
        const float* levels, int numLevels);

    /**
     * 合并计算等值线与填色多边形（性能优化：单次网格遍历 + 边交点缓存）
     *
     * 相比分别调用 ComputeContourLines + ComputeFilledContours：
     *   - 每个单元格的 4 角值和坐标只读取一次（减少 ~50% 内存访问）
     *   - 同一网格边的交点计算一次、等值线和填色共享（减少重复插值）
     */
    CMergedResult ComputeContourAndFill(
        const float* levels, int numLevels);

private:
    int _rows;
    int _cols;
    std::vector<float> _xCoords;
    std::vector<float> _yCoords;
    std::vector<float> _values;

    float GetValue(int row, int col) const;
    MPoint GetGridPoint(int row, int col) const;
    MPoint Interpolate(const MPoint& p1, float v1,
                       const MPoint& p2, float v2,
                       float level) const;

    void ComputeCellSegments(
        int row, int col, float level,
        std::vector<CSegPair>& segments) const;

    std::vector<MCurve> ChainSegments(
        std::vector<CSegPair>& segments) const;

    void ComputeCellFillEdges(
        int row, int col,
        float levelLow, float levelHigh,
        std::vector<CSegPair>& edges) const;

    void CancelOpposingEdges(
        std::vector<CSegPair>& edges,
        std::vector<CSegPair>& remaining) const;

    std::vector<MRing> TraceRings(
        std::vector<CSegPair>& edges) const;

    std::vector<MPolygon> GroupRingsToPolygons(
        std::vector<MRing>& rings) const;

    /** 对一个区间的有向边做后处理：对消→追踪→归组 */
    std::vector<MPolygon> AssemblePolygonsFromEdges(
        std::vector<CSegPair>& edges) const;
};

#endif // CONTOURGENERATOR_H
