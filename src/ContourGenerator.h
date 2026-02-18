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
 */
struct CFilledContourResult
{
    float LevelLow;
    float LevelHigh;
    std::vector<MPolygon> Polygons;
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
     * 计算等值线
     * @param levels   等值级别数组
     * @param numLevels 级别数
     * @return 每个级别对应的等值线结果
     */
    std::vector<CContourLineResult> ComputeContourLines(
        const float* levels, int numLevels);

    /**
     * 计算填色多边形
     * @param levels   等值级别数组（从小到大排列）
     * @param numLevels 级别数，产生 numLevels-1 个填色区间
     * @return 每个区间对应的填色结果
     */
    std::vector<CFilledContourResult> ComputeFilledContours(
        const float* levels, int numLevels);

private:
    int _rows;
    int _cols;
    std::vector<float> _xCoords;
    std::vector<float> _yCoords;
    std::vector<float> _values;

    /** 获取网格值 value[row][col] */
    float GetValue(int row, int col) const;

    /** 获取网格点坐标 */
    MPoint GetGridPoint(int row, int col) const;

    /** 线性插值：在两点之间找到等值点 */
    MPoint Interpolate(const MPoint& p1, float v1,
                       const MPoint& p2, float v2,
                       float level) const;

    /**
     * Marching Squares 核心：计算单元格内的等值线段
     * 一个单元格最多产生 2 条线段（鞍点时）
     */
    void ComputeCellSegments(
        int row, int col, float level,
        std::vector<std::pair<MPoint, MPoint>>& segments) const;

    /** 将散乱线段链接为连续多段线 */
    std::vector<MCurve> ChainSegments(
        std::vector<std::pair<MPoint, MPoint>>& segments) const;

    /**
     * 计算单元格内填色区间 [levelLow, levelHigh) 的有向多边形边
     * 沿单元格边界行走，收集处于"区间内"的顶点和插值交点
     */
    void ComputeCellFillEdges(
        int row, int col,
        float levelLow, float levelHigh,
        std::vector<std::pair<MPoint, MPoint>>& edges) const;

    /** 对消相邻单元格共享的反向边 */
    void CancelOpposingEdges(
        std::vector<std::pair<MPoint, MPoint>>& edges,
        std::vector<std::pair<MPoint, MPoint>>& remaining) const;

    /** 将剩余有向边追踪为闭合环 */
    std::vector<MRing> TraceRings(
        std::vector<std::pair<MPoint, MPoint>>& edges) const;

    /** 按有符号面积分类外环/内环，并归组为 MPolygon */
    std::vector<MPolygon> GroupRingsToPolygons(
        std::vector<MRing>& rings) const;
};

#endif // CONTOURGENERATOR_H
