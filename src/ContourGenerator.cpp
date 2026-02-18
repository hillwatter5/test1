/**
 * ContourGenerator.cpp
 * 功能说明：等值线与填色多边形核心算法的完整实现。
 *          算法对齐 matplotlib._contour 的 Marching Squares 逻辑。
 * 作者：ContourDev
 * 创建时间：2026-02-18
 */

#include "ContourGenerator.h"
#include <cmath>
#include <algorithm>
#include <map>
#include <set>

// ============================================================
//  Marching Squares 查找表
//
//  单元格角点编号（逆时针）：
//      c3(col, row+1) ---- c2(col+1, row+1)
//           |                    |
//      c0(col, row)   ---- c1(col+1, row)
//
//  边编号：
//      e0: c0→c1（底边）  e1: c1→c2（右边）
//      e2: c2→c3（顶边）  e3: c3→c0（左边）
//
//  Case Index = (c0>=level?1:0) | (c1>=level?2:0)
//             | (c2>=level?4:0) | (c3>=level?8:0)
// ============================================================

struct CSegmentDef
{
    int Edge1;
    int Edge2;
};

static const int SEGMENT_COUNT[16] = {
    0, 1, 1, 1, 1, 2, 1, 1,
    1, 1, 2, 1, 1, 1, 1, 0
};

static const int EDGE_CORNERS[4][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}
};

/**
 * 默认查找表。鞍点 case 5/10 给出"分离"连接方式，
 * 当中心值 >= level 时翻转为"连通"方式。
 */
static const CSegmentDef CASE_TABLE[16][2] = {
    {{-1,-1},{-1,-1}}, {{ 0, 3},{-1,-1}}, {{ 0, 1},{-1,-1}}, {{ 3, 1},{-1,-1}},
    {{ 1, 2},{-1,-1}}, {{ 0, 3},{ 1, 2}}, {{ 0, 2},{-1,-1}}, {{ 3, 2},{-1,-1}},
    {{ 2, 3},{-1,-1}}, {{ 0, 2},{-1,-1}}, {{ 0, 1},{ 2, 3}}, {{ 1, 2},{-1,-1}},
    {{ 1, 3},{-1,-1}}, {{ 0, 1},{-1,-1}}, {{ 0, 3},{-1,-1}}, {{-1,-1},{-1,-1}},
};

static const CSegmentDef SADDLE_FLIP[2][2] = {
    {{ 0, 1},{ 2, 3}},   // case  5 翻转
    {{ 0, 3},{ 1, 2}},   // case 10 翻转
};

// ============================================================
//  静态辅助函数（必须在调用者之前定义）
// ============================================================

/**
 * 规范化线性插值。
 * 为了确保同一条网格边从两个相邻单元格计算时结果一致，
 * 统一按坐标较小端 → 较大端方向插值，避免浮点运算路径不同导致结果偏差。
 */
static MPoint CanonicalInterpolate(
    const MPoint& p1, float v1,
    const MPoint& p2, float v2,
    float level)
{
    if (fabs(v2 - v1) < (float)EPSILON)
        return MPoint((p1.X + p2.X) * 0.5f, (p1.Y + p2.Y) * 0.5f);

    const MPoint* pa = &p1;
    const MPoint* pb = &p2;
    float va = v1, vb = v2;

    bool needSwap = (p1.X > p2.X + (float)EPSILON) ||
                    (fabs(p1.X - p2.X) < (float)EPSILON && p1.Y > p2.Y + (float)EPSILON);
    if (needSwap) {
        pa = &p2; pb = &p1;
        va = v2;  vb = v1;
    }

    float t = (level - va) / (vb - va);
    t = std::max(0.0f, std::min(1.0f, t));

    if (t <= 0.0f) return *pa;
    if (t >= 1.0f) return *pb;

    return MPoint(pa->X + t * (pb->X - pa->X),
                  pa->Y + t * (pb->Y - pa->Y));
}

/**
 * 从多段线尾部向前延伸
 */
static void ExtendChainForward(
    std::vector<MPoint>& chain,
    const std::vector<std::pair<MPoint, MPoint>>& segments,
    std::map<CPointKey, std::vector<int>>& endpointMap,
    std::vector<bool>& visited)
{
    bool extended = true;
    while (extended) {
        extended = false;
        CPointKey tailKey(chain.back());
        auto it = endpointMap.find(tailKey);
        if (it == endpointMap.end()) break;

        for (size_t k = 0; k < it->second.size(); k++) {
            int idx = it->second[k];
            if (visited[idx]) continue;

            CPointKey sk(segments[idx].first);
            CPointKey ek(segments[idx].second);

            if (sk == tailKey) {
                visited[idx] = true;
                chain.push_back(segments[idx].second);
                extended = true;
                break;
            } else if (ek == tailKey) {
                visited[idx] = true;
                chain.push_back(segments[idx].first);
                extended = true;
                break;
            }
        }
    }
}

/**
 * 从多段线头部向后延伸
 */
static void ExtendChainBackward(
    std::vector<MPoint>& chain,
    const std::vector<std::pair<MPoint, MPoint>>& segments,
    std::map<CPointKey, std::vector<int>>& endpointMap,
    std::vector<bool>& visited)
{
    bool extended = true;
    while (extended) {
        extended = false;
        CPointKey headKey(chain.front());
        auto it = endpointMap.find(headKey);
        if (it == endpointMap.end()) break;

        for (size_t k = 0; k < it->second.size(); k++) {
            int idx = it->second[k];
            if (visited[idx]) continue;

            CPointKey sk(segments[idx].first);
            CPointKey ek(segments[idx].second);

            if (sk == headKey) {
                visited[idx] = true;
                chain.insert(chain.begin(), segments[idx].second);
                extended = true;
                break;
            } else if (ek == headKey) {
                visited[idx] = true;
                chain.insert(chain.begin(), segments[idx].first);
                extended = true;
                break;
            }
        }
    }
}

/**
 * 沿单元格边界逆时针行走，收集 [levelLow, levelHigh) 区间内的顶点。
 * 行走顺序：c0→c1（底边）→c2（右边）→c3（顶边）→c0（左边）
 */
static void CollectFillVertices(
    const MPoint corners[4], const float vals[4], const int states[4],
    float levelLow, float levelHigh,
    std::vector<MPoint>& polyVerts)
{
    polyVerts.clear();

    auto interpolateEdge = [&](int cA, int cB, float level) -> MPoint {
        return CanonicalInterpolate(corners[cA], vals[cA],
                                    corners[cB], vals[cB], level);
    };

    const int edgeStart[4] = {0, 1, 2, 3};
    const int edgeEnd[4]   = {1, 2, 3, 0};

    for (int e = 0; e < 4; e++) {
        int cA = edgeStart[e];
        int cB = edgeEnd[e];
        float vA = vals[cA], vB = vals[cB];

        if (states[cA] == 1)
            polyVerts.push_back(corners[cA]);

        // 收集本边上与 levelLow/levelHigh 的交点，按参数 t 排序
        struct Intersection { float T; MPoint Pt; };
        std::vector<Intersection> isects;

        if ((vA < levelLow) != (vB < levelLow)) {
            float t = (levelLow - vA) / (vB - vA);
            isects.push_back({t, interpolateEdge(cA, cB, levelLow)});
        }
        if ((vA < levelHigh) != (vB < levelHigh)) {
            float t = (levelHigh - vA) / (vB - vA);
            isects.push_back({t, interpolateEdge(cA, cB, levelHigh)});
        }

        std::sort(isects.begin(), isects.end(),
                  [](const Intersection& a, const Intersection& b) {
                      return a.T < b.T;
                  });

        for (size_t k = 0; k < isects.size(); k++)
            polyVerts.push_back(isects[k].Pt);
    }
}

/**
 * 判断是否需要对鞍点单元格做拆分。
 * 对角线角点同为 between、另一对角非 between、中心也非 between 时为真。
 */
static bool NeedSaddleSplit(const int states[4], int centerState)
{
    bool diag02 = (states[0] == 1 && states[2] == 1);
    bool diag13 = (states[1] == 1 && states[3] == 1);
    bool other02 = (states[0] != 1 && states[2] != 1);
    bool other13 = (states[1] != 1 && states[3] != 1);

    if (diag02 && other13 && centerState != 1) return true;
    if (diag13 && other02 && centerState != 1) return true;
    return false;
}

/**
 * 将鞍点产生的自交多边形拆分为两个独立多边形并输出有向边。
 */
static void SplitSaddlePolygon(
    const std::vector<MPoint>& polyVerts,
    std::vector<std::pair<MPoint, MPoint>>& edges)
{
    int n = (int)polyVerts.size();
    if (n < 6) {
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            edges.push_back(std::make_pair(polyVerts[i], polyVerts[j]));
        }
        return;
    }

    int half = n / 2;
    std::vector<MPoint> g1(polyVerts.begin(), polyVerts.begin() + half);
    std::vector<MPoint> g2(polyVerts.begin() + half, polyVerts.end());

    for (size_t i = 0; i < g1.size(); i++) {
        size_t j = (i + 1) % g1.size();
        edges.push_back(std::make_pair(g1[i], g1[j]));
    }
    for (size_t i = 0; i < g2.size(); i++) {
        size_t j = (i + 1) % g2.size();
        edges.push_back(std::make_pair(g2[i], g2[j]));
    }
}

// ============================================================
//  构造与析构
// ============================================================

CContourGenerator::CContourGenerator()
    : _rows(0), _cols(0)
{
}

CContourGenerator::~CContourGenerator()
{
}

// ============================================================
//  初始化网格数据
// ============================================================

void CContourGenerator::Initialize(
    const float* xCoords, int cols,
    const float* yCoords, int rows,
    const float* values)
{
    _rows = rows;
    _cols = cols;
    _xCoords.assign(xCoords, xCoords + cols);
    _yCoords.assign(yCoords, yCoords + rows);
    _values.assign(values, values + rows * cols);
}

// ============================================================
//  基础工具方法
// ============================================================

float CContourGenerator::GetValue(int row, int col) const
{
    return _values[row * _cols + col];
}

MPoint CContourGenerator::GetGridPoint(int row, int col) const
{
    return MPoint(_xCoords[col], _yCoords[row]);
}

MPoint CContourGenerator::Interpolate(
    const MPoint& p1, float v1,
    const MPoint& p2, float v2,
    float level) const
{
    return CanonicalInterpolate(p1, v1, p2, v2, level);
}

// ============================================================
//  等值线段生成（Marching Squares 核心）
// ============================================================

/**
 * 计算单个单元格内的等值线段。
 *
 * 步骤：
 * 1. 读取四角值 → 4-bit caseIndex
 * 2. 查表确定被穿过的边对
 * 3. 鞍点（case 5/10）用中心值修正连接方式
 * 4. 对每条被穿过的边线性插值得交点
 * 5. 交点对组成线段输出
 */
void CContourGenerator::ComputeCellSegments(
    int row, int col, float level,
    std::vector<std::pair<MPoint, MPoint>>& segments) const
{
    float vals[4];
    vals[0] = GetValue(row, col);
    vals[1] = GetValue(row, col + 1);
    vals[2] = GetValue(row + 1, col + 1);
    vals[3] = GetValue(row + 1, col);

    MPoint corners[4];
    corners[0] = GetGridPoint(row, col);
    corners[1] = GetGridPoint(row, col + 1);
    corners[2] = GetGridPoint(row + 1, col + 1);
    corners[3] = GetGridPoint(row + 1, col);

    int caseIndex = 0;
    if (vals[0] >= level) caseIndex |= 1;
    if (vals[1] >= level) caseIndex |= 2;
    if (vals[2] >= level) caseIndex |= 4;
    if (vals[3] >= level) caseIndex |= 8;

    int numSegs = SEGMENT_COUNT[caseIndex];
    if (numSegs == 0) return;

    const CSegmentDef* segs = CASE_TABLE[caseIndex];

    // 鞍点修正
    if (caseIndex == 5) {
        float center = (vals[0] + vals[1] + vals[2] + vals[3]) * 0.25f;
        if (center >= level)
            segs = SADDLE_FLIP[0];
    } else if (caseIndex == 10) {
        float center = (vals[0] + vals[1] + vals[2] + vals[3]) * 0.25f;
        if (center >= level)
            segs = SADDLE_FLIP[1];
    }

    for (int s = 0; s < numSegs; s++) {
        int e1 = segs[s].Edge1;
        int e2 = segs[s].Edge2;
        if (e1 < 0 || e2 < 0) continue;

        int cA1 = EDGE_CORNERS[e1][0], cB1 = EDGE_CORNERS[e1][1];
        int cA2 = EDGE_CORNERS[e2][0], cB2 = EDGE_CORNERS[e2][1];

        MPoint pt1 = Interpolate(corners[cA1], vals[cA1],
                                 corners[cB1], vals[cB1], level);
        MPoint pt2 = Interpolate(corners[cA2], vals[cA2],
                                 corners[cB2], vals[cB2], level);

        // 过滤退化线段（角点值恰好等于等值级别时可能产生零长度线段）
        if (pt1 == pt2) continue;

        segments.push_back(std::make_pair(pt1, pt2));
    }
}

// ============================================================
//  线段链接为连续多段线
// ============================================================

/**
 * 尝试将共享端点的碎片链合并。
 * 解决角点恰好位于等值面上时产生的链断裂问题。
 */
static void MergeChainFragments(std::vector<std::vector<MPoint>>& chains)
{
    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t i = 0; i < chains.size() && !merged; i++) {
            if (chains[i].empty()) continue;
            CPointKey tailKey(chains[i].back());
            CPointKey headKey(chains[i].front());

            for (size_t j = i + 1; j < chains.size() && !merged; j++) {
                if (chains[j].empty()) continue;
                CPointKey jHead(chains[j].front());
                CPointKey jTail(chains[j].back());

                if (tailKey == jHead) {
                    // i 的尾部接 j 的头部
                    chains[i].insert(chains[i].end(),
                                     chains[j].begin() + 1, chains[j].end());
                    chains[j].clear();
                    merged = true;
                } else if (headKey == jTail) {
                    // j 的尾部接 i 的头部
                    chains[j].insert(chains[j].end(),
                                     chains[i].begin() + 1, chains[i].end());
                    chains[i].swap(chains[j]);
                    chains[j].clear();
                    merged = true;
                } else if (tailKey == jTail) {
                    // 翻转 j 后拼到 i 尾部
                    std::reverse(chains[j].begin(), chains[j].end());
                    chains[i].insert(chains[i].end(),
                                     chains[j].begin() + 1, chains[j].end());
                    chains[j].clear();
                    merged = true;
                } else if (headKey == jHead) {
                    // 翻转 i 后拼到 j 头部
                    std::reverse(chains[i].begin(), chains[i].end());
                    chains[i].insert(chains[i].end(),
                                     chains[j].begin() + 1, chains[j].end());
                    chains[j].clear();
                    merged = true;
                }
            }
        }
    }

    // 清除空链
    std::vector<std::vector<MPoint>> result;
    for (size_t i = 0; i < chains.size(); i++) {
        if (!chains[i].empty())
            result.push_back(chains[i]);
    }
    chains.swap(result);
}

/**
 * 用量化键判等去除多段线中的连续重复点（比浮点 EPSILON 更可靠）
 */
static void RemoveDuplicatePoints(std::vector<MPoint>& pts)
{
    if (pts.size() < 2) return;
    std::vector<MPoint> cleaned;
    cleaned.push_back(pts[0]);
    CPointKey prevKey(pts[0]);
    for (size_t i = 1; i < pts.size(); i++) {
        CPointKey curKey(pts[i]);
        if (!(curKey == prevKey)) {
            cleaned.push_back(pts[i]);
            prevKey = curKey;
        }
    }
    pts.swap(cleaned);
}

/**
 * 建立端点→线段索引映射，从任一未访问线段出发向两端延伸，
 * 直到无法继续或形成闭环，输出一条完整多段线。
 * 最后对碎片链尝试再次合并（处理角点恰好在等值面上的情况）。
 */
std::vector<MCurve> CContourGenerator::ChainSegments(
    std::vector<std::pair<MPoint, MPoint>>& segments) const
{
    int n = (int)segments.size();
    if (n == 0) return {};

    std::map<CPointKey, std::vector<int>> endpointMap;
    for (int i = 0; i < n; i++) {
        endpointMap[CPointKey(segments[i].first)].push_back(i);
        endpointMap[CPointKey(segments[i].second)].push_back(i);
    }

    std::vector<bool> visited(n, false);
    std::vector<std::vector<MPoint>> chains;

    for (int i = 0; i < n; i++) {
        if (visited[i]) continue;
        visited[i] = true;

        std::vector<MPoint> chain;
        chain.push_back(segments[i].first);
        chain.push_back(segments[i].second);

        ExtendChainForward(chain, segments, endpointMap, visited);
        ExtendChainBackward(chain, segments, endpointMap, visited);
        RemoveDuplicatePoints(chain);

        if (chain.size() >= 2)
            chains.push_back(chain);
    }

    // 二次合并：尝试将共享端点的碎片链首尾拼接
    MergeChainFragments(chains);

    // 合并后再次去重
    std::vector<MCurve> curves;
    for (size_t i = 0; i < chains.size(); i++) {
        RemoveDuplicatePoints(chains[i]);
        if (chains[i].size() >= 2) {
            MCurve curve;
            curve.Points = chains[i];
            curves.push_back(curve);
        }
    }
    return curves;
}

// ============================================================
//  等值线计算入口
// ============================================================

std::vector<CContourLineResult> CContourGenerator::ComputeContourLines(
    const float* levels, int numLevels)
{
    std::vector<CContourLineResult> results;

    for (int lv = 0; lv < numLevels; lv++) {
        float level = levels[lv];
        std::vector<std::pair<MPoint, MPoint>> allSegments;

        for (int row = 0; row < _rows - 1; row++)
            for (int col = 0; col < _cols - 1; col++)
                ComputeCellSegments(row, col, level, allSegments);

        CContourLineResult result;
        result.Level = level;
        result.Curves = ChainSegments(allSegments);
        results.push_back(result);
    }

    return results;
}

// ============================================================
//  填色多边形：单元格有向边生成
// ============================================================

/**
 * 对单元格角点按值分为 below(0)/between(1)/above(2) 三种状态，
 * 沿边界行走收集 between 区域的顶点（角点 + 等值交点），
 * 输出相邻顶点构成的有向边。鞍点时做拆分。
 */
void CContourGenerator::ComputeCellFillEdges(
    int row, int col,
    float levelLow, float levelHigh,
    std::vector<std::pair<MPoint, MPoint>>& edges) const
{
    float vals[4];
    vals[0] = GetValue(row, col);
    vals[1] = GetValue(row, col + 1);
    vals[2] = GetValue(row + 1, col + 1);
    vals[3] = GetValue(row + 1, col);

    MPoint corners[4];
    corners[0] = GetGridPoint(row, col);
    corners[1] = GetGridPoint(row, col + 1);
    corners[2] = GetGridPoint(row + 1, col + 1);
    corners[3] = GetGridPoint(row + 1, col);

    int states[4];
    for (int i = 0; i < 4; i++) {
        if (vals[i] < levelLow)        states[i] = 0;
        else if (vals[i] < levelHigh)  states[i] = 1;
        else                           states[i] = 2;
    }

    // 快速排除：无 between 角点且无边穿越
    bool anyBetween = false;
    bool hasCross = false;
    for (int i = 0; i < 4; i++) {
        if (states[i] == 1) anyBetween = true;
        int j = (i + 1) % 4;
        float vA = vals[i], vB = vals[j];
        if ((vA < levelLow) != (vB < levelLow)) hasCross = true;
        if ((vA < levelHigh) != (vB < levelHigh)) hasCross = true;
    }
    if (!anyBetween && !hasCross) return;

    std::vector<MPoint> polyVerts;
    CollectFillVertices(corners, vals, states, levelLow, levelHigh, polyVerts);
    RemoveDuplicatePoints(polyVerts);
    if (polyVerts.size() < 3) return;

    float center = (vals[0] + vals[1] + vals[2] + vals[3]) * 0.25f;
    int centerState = (center < levelLow) ? 0 : (center < levelHigh) ? 1 : 2;

    if (NeedSaddleSplit(states, centerState) && polyVerts.size() >= 6) {
        SplitSaddlePolygon(polyVerts, edges);
    } else {
        for (size_t i = 0; i < polyVerts.size(); i++) {
            size_t j = (i + 1) % polyVerts.size();
            edges.push_back(std::make_pair(polyVerts[i], polyVerts[j]));
        }
    }
}

// ============================================================
//  填色计算入口
// ============================================================

std::vector<CFilledContourResult> CContourGenerator::ComputeFilledContours(
    const float* levels, int numLevels)
{
    std::vector<CFilledContourResult> results;
    if (numLevels < 2) return results;

    for (int lv = 0; lv < numLevels - 1; lv++) {
        float levelLow = levels[lv];
        float levelHigh = levels[lv + 1];
        std::vector<std::pair<MPoint, MPoint>> allEdges;

        for (int row = 0; row < _rows - 1; row++)
            for (int col = 0; col < _cols - 1; col++)
                ComputeCellFillEdges(row, col, levelLow, levelHigh, allEdges);

        std::vector<std::pair<MPoint, MPoint>> remaining;
        CancelOpposingEdges(allEdges, remaining);
        std::vector<MRing> rings = TraceRings(remaining);

        CFilledContourResult result;
        result.LevelLow = levelLow;
        result.LevelHigh = levelHigh;
        result.Polygons = GroupRingsToPolygons(rings);
        results.push_back(result);
    }

    return results;
}

// ============================================================
//  对消相邻单元格的反向共享边
// ============================================================

/**
 * 两个相邻单元格共享网格边时，各自填色多边形产生方向相反的有向边。
 * 对消后剩余的边即为整个填色区域的外边界。
 */
void CContourGenerator::CancelOpposingEdges(
    std::vector<std::pair<MPoint, MPoint>>& edges,
    std::vector<std::pair<MPoint, MPoint>>& remaining) const
{
    std::map<CEdgeKey, int> edgeCount;
    std::map<CEdgeKey, std::pair<MPoint, MPoint>> edgeData;

    for (size_t i = 0; i < edges.size(); i++) {
        CPointKey ks(edges[i].first);
        CPointKey ke(edges[i].second);
        CEdgeKey fwd(ks, ke);
        CEdgeKey rev(ke, ks);

        if (edgeCount.count(rev) && edgeCount[rev] > 0) {
            edgeCount[rev]--;
        } else {
            edgeCount[fwd]++;
            edgeData[fwd] = edges[i];
        }
    }

    remaining.clear();
    for (auto& kv : edgeCount) {
        for (int c = 0; c < kv.second; c++)
            remaining.push_back(edgeData[kv.first]);
    }
}

// ============================================================
//  追踪闭合环
// ============================================================

/**
 * 建立起点→终点邻接表，逐条追踪直到回到起点形成闭合环。
 */
std::vector<MRing> CContourGenerator::TraceRings(
    std::vector<std::pair<MPoint, MPoint>>& edges) const
{
    std::vector<MRing> rings;
    if (edges.empty()) return rings;

    std::map<CPointKey, std::vector<std::pair<CPointKey, MPoint>>> adjacency;
    std::map<CPointKey, MPoint> keyToPoint;

    for (size_t i = 0; i < edges.size(); i++) {
        CPointKey ks(edges[i].first);
        CPointKey ke(edges[i].second);
        adjacency[ks].push_back(std::make_pair(ke, edges[i].second));
        keyToPoint[ks] = edges[i].first;
        keyToPoint[ke] = edges[i].second;
    }

    int maxSteps = (int)edges.size() + 1;

    while (!adjacency.empty()) {
        CPointKey startKey = adjacency.begin()->first;
        MPoint startPt = keyToPoint[startKey];

        MRing ring;
        ring.AddPoint(startPt);
        CPointKey currentKey = startKey;
        bool ok = true;

        for (int step = 0; step < maxSteps; step++) {
            auto adj = adjacency.find(currentKey);
            if (adj == adjacency.end() || adj->second.empty()) {
                ok = false;
                break;
            }

            CPointKey nextKey = adj->second.back().first;
            MPoint nextPt = adj->second.back().second;
            adj->second.pop_back();
            if (adj->second.empty()) adjacency.erase(adj);

            if (nextKey == startKey) break;

            ring.AddPoint(nextPt);
            currentKey = nextKey;
        }

        if (ok && ring.Points.size() >= 3) {
            RemoveDuplicatePoints(ring.Points);
            if (ring.Points.size() >= 3)
                rings.push_back(ring);
        }
    }

    return rings;
}

// ============================================================
//  按拓扑归组环为 MPolygon
// ============================================================

/**
 * 正面积 → 逆时针 → 外环，负面积 → 顺时针 → 内环。
 * 内环分配给包含它的面积最小的外环。
 */
std::vector<MPolygon> CContourGenerator::GroupRingsToPolygons(
    std::vector<MRing>& rings) const
{
    if (rings.empty()) return {};

    std::vector<int> outerIdx, innerIdx;
    for (size_t i = 0; i < rings.size(); i++) {
        if (rings[i].SignedArea() > 0)
            outerIdx.push_back((int)i);
        else
            innerIdx.push_back((int)i);
    }

    if (outerIdx.empty()) {
        for (size_t i = 0; i < rings.size(); i++) {
            std::reverse(rings[i].Points.begin(), rings[i].Points.end());
            outerIdx.push_back((int)i);
        }
        innerIdx.clear();
    }

    std::vector<MPolygon> polyVec(outerIdx.size());
    for (size_t i = 0; i < outerIdx.size(); i++)
        polyVec[i].OuterRings.push_back(rings[outerIdx[i]]);

    for (size_t i = 0; i < innerIdx.size(); i++) {
        MPoint testPt = rings[innerIdx[i]].Points[0];
        int bestOuter = -1;
        double bestArea = 1e30;
        for (size_t j = 0; j < outerIdx.size(); j++) {
            if (rings[outerIdx[j]].Contain(testPt)) {
                double a = fabs(rings[outerIdx[j]].SignedArea());
                if (a < bestArea) { bestArea = a; bestOuter = (int)j; }
            }
        }
        if (bestOuter >= 0)
            polyVec[bestOuter].InnerRings.push_back(rings[innerIdx[i]]);
    }

    return polyVec;
}
