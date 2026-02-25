/**
 * ContourGenerator.cpp
 * 功能说明：等值线与填色多边形核心算法的完整实现。
 *          算法对齐 matplotlib._contour 的 Marching Squares 逻辑。
 *          支持单次网格遍历 + 边交点缓存的合并计算模式。
 * 作者：ContourDev
 * 创建时间：2026-02-18
 */

#include "ContourGenerator.h"
#include <cmath>
#include <algorithm>
#include <map>

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

struct CSegmentDef { int Edge1; int Edge2; };

static const int SEGMENT_COUNT[16] = {
    0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 0
};

static const int EDGE_CORNERS[4][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}
};

static const CSegmentDef CASE_TABLE[16][2] = {
    {{-1,-1},{-1,-1}}, {{ 0, 3},{-1,-1}}, {{ 0, 1},{-1,-1}}, {{ 3, 1},{-1,-1}},
    {{ 1, 2},{-1,-1}}, {{ 0, 3},{ 1, 2}}, {{ 0, 2},{-1,-1}}, {{ 3, 2},{-1,-1}},
    {{ 2, 3},{-1,-1}}, {{ 0, 2},{-1,-1}}, {{ 0, 1},{ 2, 3}}, {{ 1, 2},{-1,-1}},
    {{ 1, 3},{-1,-1}}, {{ 0, 1},{-1,-1}}, {{ 0, 3},{-1,-1}}, {{-1,-1},{-1,-1}},
};

static const CSegmentDef SADDLE_FLIP[2][2] = {
    {{ 0, 1},{ 2, 3}}, {{ 0, 3},{ 1, 2}},
};

// ============================================================
//  公共静态辅助函数
// ============================================================

static MPoint CanonicalInterpolate(
    const MPoint& p1, float v1, const MPoint& p2, float v2, float level)
{
    if (fabs(v2 - v1) < (float)EPSILON)
        return MPoint((p1.X + p2.X) * 0.5f, (p1.Y + p2.Y) * 0.5f);

    const MPoint* pa = &p1;
    const MPoint* pb = &p2;
    float va = v1, vb = v2;

    bool needSwap = (p1.X > p2.X + (float)EPSILON) ||
                    (fabs(p1.X - p2.X) < (float)EPSILON && p1.Y > p2.Y + (float)EPSILON);
    if (needSwap) { pa = &p2; pb = &p1; va = v2; vb = v1; }

    float t = (level - va) / (vb - va);
    t = std::max(0.0f, std::min(1.0f, t));
    if (t <= 0.0f) return *pa;
    if (t >= 1.0f) return *pb;
    return MPoint(pa->X + t * (pb->X - pa->X), pa->Y + t * (pb->Y - pa->Y));
}

static void RemoveDuplicatePoints(std::vector<MPoint>& pts)
{
    if (pts.size() < 2) return;
    std::vector<MPoint> cleaned;
    cleaned.reserve(pts.size());
    cleaned.push_back(pts[0]);
    CPointKey prevKey(pts[0]);
    for (size_t i = 1; i < pts.size(); i++) {
        CPointKey curKey(pts[i]);
        if (!(curKey == prevKey)) { cleaned.push_back(pts[i]); prevKey = curKey; }
    }
    pts.swap(cleaned);
}

static void ExtendChainForward(
    std::vector<MPoint>& chain,
    const std::vector<CSegPair>& segments,
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
            CPointKey sk(segments[idx].first), ek(segments[idx].second);
            if (sk == tailKey) {
                visited[idx] = true; chain.push_back(segments[idx].second); extended = true; break;
            } else if (ek == tailKey) {
                visited[idx] = true; chain.push_back(segments[idx].first); extended = true; break;
            }
        }
    }
}

static void ExtendChainBackward(
    std::vector<MPoint>& chain,
    const std::vector<CSegPair>& segments,
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
            CPointKey sk(segments[idx].first), ek(segments[idx].second);
            if (sk == headKey) {
                visited[idx] = true; chain.insert(chain.begin(), segments[idx].second); extended = true; break;
            } else if (ek == headKey) {
                visited[idx] = true; chain.insert(chain.begin(), segments[idx].first); extended = true; break;
            }
        }
    }
}

static void MergeChainFragments(std::vector<std::vector<MPoint>>& chains)
{
    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t i = 0; i < chains.size() && !merged; i++) {
            if (chains[i].empty()) continue;
            CPointKey tailI(chains[i].back()), headI(chains[i].front());
            for (size_t j = i + 1; j < chains.size() && !merged; j++) {
                if (chains[j].empty()) continue;
                CPointKey headJ(chains[j].front()), tailJ(chains[j].back());
                if (tailI == headJ) {
                    chains[i].insert(chains[i].end(), chains[j].begin() + 1, chains[j].end());
                    chains[j].clear(); merged = true;
                } else if (headI == tailJ) {
                    chains[j].insert(chains[j].end(), chains[i].begin() + 1, chains[i].end());
                    chains[i].swap(chains[j]); chains[j].clear(); merged = true;
                } else if (tailI == tailJ) {
                    std::reverse(chains[j].begin(), chains[j].end());
                    chains[i].insert(chains[i].end(), chains[j].begin() + 1, chains[j].end());
                    chains[j].clear(); merged = true;
                } else if (headI == headJ) {
                    std::reverse(chains[i].begin(), chains[i].end());
                    chains[i].insert(chains[i].end(), chains[j].begin() + 1, chains[j].end());
                    chains[j].clear(); merged = true;
                }
            }
        }
    }
    std::vector<std::vector<MPoint>> result;
    for (size_t i = 0; i < chains.size(); i++)
        if (!chains[i].empty()) result.push_back(chains[i]);
    chains.swap(result);
}

static bool NeedSaddleSplit(const int states[4], int centerState)
{
    bool d02 = (states[0] == 1 && states[2] == 1);
    bool d13 = (states[1] == 1 && states[3] == 1);
    bool o02 = (states[0] != 1 && states[2] != 1);
    bool o13 = (states[1] != 1 && states[3] != 1);
    return (d02 && o13 && centerState != 1) || (d13 && o02 && centerState != 1);
}

static void SplitSaddlePolygon(
    const std::vector<MPoint>& polyVerts, std::vector<CSegPair>& edges)
{
    int n = (int)polyVerts.size();
    if (n < 6) {
        for (int i = 0; i < n; i++)
            edges.push_back({polyVerts[i], polyVerts[(i + 1) % n]});
        return;
    }
    int half = n / 2;
    std::vector<MPoint> g1(polyVerts.begin(), polyVerts.begin() + half);
    std::vector<MPoint> g2(polyVerts.begin() + half, polyVerts.end());
    for (size_t i = 0; i < g1.size(); i++)
        edges.push_back({g1[i], g1[(i + 1) % g1.size()]});
    for (size_t i = 0; i < g2.size(); i++)
        edges.push_back({g2[i], g2[(i + 1) % g2.size()]});
}

// ============================================================
//  核心静态函数：基于预读单元格数据生成等值线段
// ============================================================

/**
 * 从预读的角点数据生成等值线段。
 * 可选地使用缓存的边交点（cache 非空时），避免重复插值。
 */
static void GenerateCellContour(
    const float vals[4], const MPoint corners[4], float level,
    int levelIdx, int numLevels,
    const std::vector<CCachedIsect>* cache,
    std::vector<CSegPair>& segments)
{
    int caseIndex = 0;
    if (vals[0] >= level) caseIndex |= 1;
    if (vals[1] >= level) caseIndex |= 2;
    if (vals[2] >= level) caseIndex |= 4;
    if (vals[3] >= level) caseIndex |= 8;

    int numSegs = SEGMENT_COUNT[caseIndex];
    if (numSegs == 0) return;

    const CSegmentDef* segs = CASE_TABLE[caseIndex];
    if (caseIndex == 5 || caseIndex == 10) {
        float center = (vals[0] + vals[1] + vals[2] + vals[3]) * 0.25f;
        if (center >= level)
            segs = (caseIndex == 5) ? SADDLE_FLIP[0] : SADDLE_FLIP[1];
    }

    for (int s = 0; s < numSegs; s++) {
        int e1 = segs[s].Edge1, e2 = segs[s].Edge2;
        if (e1 < 0 || e2 < 0) continue;

        MPoint pt1, pt2;

        // 从缓存查找或现场计算交点
        if (cache && levelIdx >= 0) {
            pt1 = (*cache)[e1 * numLevels + levelIdx].Pt;
            pt2 = (*cache)[e2 * numLevels + levelIdx].Pt;
        } else {
            int cA1 = EDGE_CORNERS[e1][0], cB1 = EDGE_CORNERS[e1][1];
            int cA2 = EDGE_CORNERS[e2][0], cB2 = EDGE_CORNERS[e2][1];
            pt1 = CanonicalInterpolate(corners[cA1], vals[cA1], corners[cB1], vals[cB1], level);
            pt2 = CanonicalInterpolate(corners[cA2], vals[cA2], corners[cB2], vals[cB2], level);
        }

        if (!(pt1 == pt2))
            segments.push_back({pt1, pt2});
    }
}

// ============================================================
//  核心静态函数：基于预读单元格数据生成填色有向边
// ============================================================

/**
 * 沿单元格边界行走收集填色区间 [levelLow, levelHigh) 的多边形顶点。
 * 支持缓存查找：lowIdx/highIdx 为交点在缓存中的 level 下标，-1 表示无交点。
 */
static void CollectFillVerticesEx(
    const MPoint corners[4], const float vals[4], const int states[4],
    float levelLow, float levelHigh,
    int lowIdx, int highIdx, int numLevels,
    const std::vector<CCachedIsect>* cache,
    std::vector<MPoint>& polyVerts)
{
    polyVerts.clear();
    const int edgeStart[4] = {0, 1, 2, 3};
    const int edgeEnd[4]   = {1, 2, 3, 0};

    for (int e = 0; e < 4; e++) {
        int cA = edgeStart[e], cB = edgeEnd[e];
        float vA = vals[cA], vB = vals[cB];

        if (states[cA] == 1)
            polyVerts.push_back(corners[cA]);

        struct Isect { float T; MPoint Pt; };
        Isect isects[2];
        int numIsects = 0;

        // 与 levelLow 的交叉
        if ((vA < levelLow) != (vB < levelLow)) {
            float t = (levelLow - vA) / (vB - vA);
            MPoint pt;
            if (cache && lowIdx >= 0)
                pt = (*cache)[e * numLevels + lowIdx].Pt;
            else
                pt = CanonicalInterpolate(corners[cA], vA, corners[cB], vB, levelLow);
            isects[numIsects++] = {t, pt};
        }
        // 与 levelHigh 的交叉
        if ((vA < levelHigh) != (vB < levelHigh)) {
            float t = (levelHigh - vA) / (vB - vA);
            MPoint pt;
            if (cache && highIdx >= 0 && highIdx < numLevels)
                pt = (*cache)[e * numLevels + highIdx].Pt;
            else
                pt = CanonicalInterpolate(corners[cA], vA, corners[cB], vB, levelHigh);
            isects[numIsects++] = {t, pt};
        }

        if (numIsects == 2 && isects[0].T > isects[1].T)
            std::swap(isects[0], isects[1]);

        for (int k = 0; k < numIsects; k++)
            polyVerts.push_back(isects[k].Pt);
    }
}

/**
 * 从预读单元格数据生成填色区间 [levelLow, levelHigh) 的有向边。
 * polyVertsBuf 是可复用的临时缓冲，避免每个单元格重新分配。
 */
static void GenerateCellFill(
    const float vals[4], const MPoint corners[4],
    float levelLow, float levelHigh,
    int lowIdx, int highIdx, int numLevels,
    const std::vector<CCachedIsect>* cache,
    std::vector<MPoint>& polyVertsBuf,
    std::vector<CSegPair>& edges)
{
    int states[4];
    for (int i = 0; i < 4; i++) {
        if (vals[i] < levelLow)        states[i] = 0;
        else if (vals[i] < levelHigh)  states[i] = 1;
        else                           states[i] = 2;
    }

    bool anyBetween = false, hasCross = false;
    for (int i = 0; i < 4; i++) {
        if (states[i] == 1) anyBetween = true;
        int j = (i + 1) % 4;
        if ((vals[i] < levelLow) != (vals[j] < levelLow)) hasCross = true;
        if ((vals[i] < levelHigh) != (vals[j] < levelHigh)) hasCross = true;
    }
    if (!anyBetween && !hasCross) return;

    CollectFillVerticesEx(corners, vals, states, levelLow, levelHigh,
                          lowIdx, highIdx, numLevels, cache, polyVertsBuf);
    RemoveDuplicatePoints(polyVertsBuf);
    if (polyVertsBuf.size() < 3) return;

    float center = (vals[0] + vals[1] + vals[2] + vals[3]) * 0.25f;
    int centerState = (center < levelLow) ? 0 : (center < levelHigh) ? 1 : 2;

    if (NeedSaddleSplit(states, centerState) && polyVertsBuf.size() >= 6) {
        SplitSaddlePolygon(polyVertsBuf, edges);
    } else {
        for (size_t i = 0; i < polyVertsBuf.size(); i++)
            edges.push_back({polyVertsBuf[i], polyVertsBuf[(i + 1) % polyVertsBuf.size()]});
    }
}

// ============================================================
//  构造、析构、初始化、工具方法
// ============================================================

CContourGenerator::CContourGenerator() : _rows(0), _cols(0) {}
CContourGenerator::~CContourGenerator() {}

void CContourGenerator::Initialize(
    const float* xCoords, int cols,
    const float* yCoords, int rows,
    const float* values)
{
    _rows = rows; _cols = cols;
    _xCoords.assign(xCoords, xCoords + cols);
    _yCoords.assign(yCoords, yCoords + rows);
    _values.assign(values, values + rows * cols);
}

float CContourGenerator::GetValue(int row, int col) const
{
    return _values[row * _cols + col];
}

MPoint CContourGenerator::GetGridPoint(int row, int col) const
{
    return MPoint(_xCoords[col], _yCoords[row]);
}

MPoint CContourGenerator::Interpolate(
    const MPoint& p1, float v1, const MPoint& p2, float v2, float level) const
{
    return CanonicalInterpolate(p1, v1, p2, v2, level);
}

// ============================================================
//  兼容旧接口的成员方法（委托到静态函数）
// ============================================================

void CContourGenerator::ComputeCellSegments(
    int row, int col, float level, std::vector<CSegPair>& segments) const
{
    float vals[4] = { GetValue(row,col), GetValue(row,col+1),
                      GetValue(row+1,col+1), GetValue(row+1,col) };
    MPoint corners[4] = { GetGridPoint(row,col), GetGridPoint(row,col+1),
                          GetGridPoint(row+1,col+1), GetGridPoint(row+1,col) };
    GenerateCellContour(vals, corners, level, -1, 0, nullptr, segments);
}

void CContourGenerator::ComputeCellFillEdges(
    int row, int col, float levelLow, float levelHigh,
    std::vector<CSegPair>& edges) const
{
    float vals[4] = { GetValue(row,col), GetValue(row,col+1),
                      GetValue(row+1,col+1), GetValue(row+1,col) };
    MPoint corners[4] = { GetGridPoint(row,col), GetGridPoint(row,col+1),
                          GetGridPoint(row+1,col+1), GetGridPoint(row+1,col) };
    std::vector<MPoint> buf;
    GenerateCellFill(vals, corners, levelLow, levelHigh,
                     -1, -1, 0, nullptr, buf, edges);
}

// ============================================================
//  线段链接
// ============================================================

std::vector<MCurve> CContourGenerator::ChainSegments(
    std::vector<CSegPair>& segments) const
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
        if (chain.size() >= 2) chains.push_back(chain);
    }

    MergeChainFragments(chains);

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
//  边对消、闭合环追踪、拓扑归组
// ============================================================

void CContourGenerator::CancelOpposingEdges(
    std::vector<CSegPair>& edges, std::vector<CSegPair>& remaining) const
{
    std::map<CEdgeKey, int> edgeCount;
    std::map<CEdgeKey, CSegPair> edgeData;

    for (size_t i = 0; i < edges.size(); i++) {
        CEdgeKey fwd(CPointKey(edges[i].first), CPointKey(edges[i].second));
        CEdgeKey rev(CPointKey(edges[i].second), CPointKey(edges[i].first));
        if (edgeCount.count(rev) && edgeCount[rev] > 0) {
            edgeCount[rev]--;
        } else {
            edgeCount[fwd]++;
            edgeData[fwd] = edges[i];
        }
    }

    remaining.clear();
    for (auto& kv : edgeCount)
        for (int c = 0; c < kv.second; c++)
            remaining.push_back(edgeData[kv.first]);
}

std::vector<MRing> CContourGenerator::TraceRings(
    std::vector<CSegPair>& edges) const
{
    std::vector<MRing> rings;
    if (edges.empty()) return rings;

    std::map<CPointKey, std::vector<std::pair<CPointKey, MPoint>>> adjacency;
    std::map<CPointKey, MPoint> keyToPoint;
    for (size_t i = 0; i < edges.size(); i++) {
        CPointKey ks(edges[i].first), ke(edges[i].second);
        adjacency[ks].push_back({ke, edges[i].second});
        keyToPoint[ks] = edges[i].first;
        keyToPoint[ke] = edges[i].second;
    }

    int maxSteps = (int)edges.size() + 1;
    while (!adjacency.empty()) {
        CPointKey startKey = adjacency.begin()->first;
        MRing ring;
        ring.AddPoint(keyToPoint[startKey]);
        CPointKey cur = startKey;
        bool ok = true;
        for (int step = 0; step < maxSteps; step++) {
            auto adj = adjacency.find(cur);
            if (adj == adjacency.end() || adj->second.empty()) { ok = false; break; }
            auto next = adj->second.back();
            adj->second.pop_back();
            if (adj->second.empty()) adjacency.erase(adj);
            if (next.first == startKey) break;
            ring.AddPoint(next.second);
            cur = next.first;
        }
        if (ok && ring.Points.size() >= 3) {
            RemoveDuplicatePoints(ring.Points);
            if (ring.Points.size() >= 3) rings.push_back(ring);
        }
    }
    return rings;
}

std::vector<MPolygon> CContourGenerator::GroupRingsToPolygons(
    std::vector<MRing>& rings) const
{
    if (rings.empty()) return {};
    std::vector<int> outerIdx, innerIdx;
    for (size_t i = 0; i < rings.size(); i++) {
        if (rings[i].SignedArea() > 0) outerIdx.push_back((int)i);
        else                           innerIdx.push_back((int)i);
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
        int best = -1; double bestA = 1e30;
        for (size_t j = 0; j < outerIdx.size(); j++) {
            if (rings[outerIdx[j]].Contain(testPt)) {
                double a = fabs(rings[outerIdx[j]].SignedArea());
                if (a < bestA) { bestA = a; best = (int)j; }
            }
        }
        if (best >= 0) polyVec[best].InnerRings.push_back(rings[innerIdx[i]]);
    }
    return polyVec;
}

std::vector<MPolygon> CContourGenerator::AssemblePolygonsFromEdges(
    std::vector<CSegPair>& edges) const
{
    std::vector<CSegPair> remaining;
    CancelOpposingEdges(edges, remaining);
    std::vector<MRing> rings = TraceRings(remaining);
    return GroupRingsToPolygons(rings);
}

// ============================================================
//  独立接口：仅计算等值线
// ============================================================

std::vector<CContourLineResult> CContourGenerator::ComputeContourLines(
    const float* levels, int numLevels)
{
    std::vector<CContourLineResult> results(numLevels);
    for (int lv = 0; lv < numLevels; lv++) {
        results[lv].Level = levels[lv];
        std::vector<CSegPair> allSegs;
        for (int r = 0; r < _rows - 1; r++)
            for (int c = 0; c < _cols - 1; c++)
                ComputeCellSegments(r, c, levels[lv], allSegs);
        results[lv].Curves = ChainSegments(allSegs);
    }
    return results;
}

// ============================================================
//  独立接口：仅计算填色（已修正为 N+1 个区间）
// ============================================================

std::vector<CFilledContourResult> CContourGenerator::ComputeFilledContours(
    const float* levels, int numLevels)
{
    if (numLevels < 1) return {};
    int numIntervals = numLevels + 1;
    std::vector<CFilledContourResult> results(numIntervals);

    for (int i = 0; i < numIntervals; i++) {
        float lo = (i == 0) ? -FLT_MAX : levels[i - 1];
        float hi = (i == numLevels) ? FLT_MAX : levels[i];
        results[i].LevelLow = lo;
        results[i].LevelHigh = hi;
        results[i].ColorIndex = i;

        std::vector<CSegPair> allEdges;
        for (int r = 0; r < _rows - 1; r++)
            for (int c = 0; c < _cols - 1; c++)
                ComputeCellFillEdges(r, c, lo, hi, allEdges);

        results[i].Polygons = AssemblePolygonsFromEdges(allEdges);
    }
    return results;
}

// ============================================================
//  合并计算：单次遍历 + 边交点缓存（性能优化核心）
// ============================================================

/**
 * 合并计算等值线与填色多边形。
 *
 * 性能优化策略：
 * 1. 单元格的 4 角值和坐标只读取一次（N+1 次降为 1 次）
 * 2. 预计算 4 条边 × N 个 level 的交点缓存表（edgeCache）
 * 3. 等值线和填色共享同一缓存表，避免重复插值
 *
 * 对比分别调用：
 *   分别调用 → 读取 (N + N+1) × Cells 次单元格数据
 *   合并调用 → 读取 1 × Cells 次单元格数据 + 计算 4N 个缓存交点
 */
CMergedResult CContourGenerator::ComputeContourAndFill(
    const float* levels, int numLevels)
{
    CMergedResult result;
    if (numLevels < 1) return result;

    int numIntervals = numLevels + 1;

    // 分配每个 level/区间的线段和边存储
    std::vector<std::vector<CSegPair>> segsByLevel(numLevels);
    std::vector<std::vector<CSegPair>> edgesByInterval(numIntervals);

    // 构建扩展边界数组：[-FLT_MAX, levels[0], ..., levels[N-1], FLT_MAX]
    std::vector<float> bounds(numIntervals + 1);
    bounds[0] = -FLT_MAX;
    for (int i = 0; i < numLevels; i++) bounds[i + 1] = levels[i];
    bounds[numIntervals] = FLT_MAX;

    // 预分配可复用缓冲
    std::vector<CCachedIsect> edgeCache(4 * numLevels);
    std::vector<MPoint> polyVertsBuf;
    polyVertsBuf.reserve(16);

    // ---- 单次网格遍历 ----
    for (int row = 0; row < _rows - 1; row++) {
        for (int col = 0; col < _cols - 1; col++) {
            // 步骤 1：只读取一次单元格数据
            float vals[4] = {
                GetValue(row, col),     GetValue(row, col + 1),
                GetValue(row + 1, col + 1), GetValue(row + 1, col)
            };
            MPoint corners[4] = {
                GetGridPoint(row, col),     GetGridPoint(row, col + 1),
                GetGridPoint(row + 1, col + 1), GetGridPoint(row + 1, col)
            };

            // 步骤 2：预计算 4 边 × N 级别的交点缓存
            for (int e = 0; e < 4; e++) {
                int cA = EDGE_CORNERS[e][0], cB = EDGE_CORNERS[e][1];
                for (int lv = 0; lv < numLevels; lv++) {
                    int idx = e * numLevels + lv;
                    edgeCache[idx].Has =
                        (vals[cA] < levels[lv]) != (vals[cB] < levels[lv]);
                    if (edgeCache[idx].Has) {
                        edgeCache[idx].Pt = CanonicalInterpolate(
                            corners[cA], vals[cA], corners[cB], vals[cB], levels[lv]);
                    }
                }
            }

            // 步骤 3：生成各 level 的等值线段（从缓存查找交点）
            for (int lv = 0; lv < numLevels; lv++) {
                GenerateCellContour(vals, corners, levels[lv],
                                    lv, numLevels, &edgeCache,
                                    segsByLevel[lv]);
            }

            // 步骤 4：生成各区间的填色有向边（从缓存查找交点）
            for (int i = 0; i < numIntervals; i++) {
                int lowIdx  = i - 1;                  // -1 表示 -FLT_MAX，无缓存
                int highIdx = (i < numLevels) ? i : -1; // -1 表示 FLT_MAX，无缓存
                GenerateCellFill(vals, corners, bounds[i], bounds[i + 1],
                                 lowIdx, highIdx, numLevels, &edgeCache,
                                 polyVertsBuf, edgesByInterval[i]);
            }
        }
    }

    // ---- 后处理 ----
    result.ContourLines.resize(numLevels);
    for (int lv = 0; lv < numLevels; lv++) {
        result.ContourLines[lv].Level = levels[lv];
        result.ContourLines[lv].Curves = ChainSegments(segsByLevel[lv]);
    }

    result.FilledContours.resize(numIntervals);
    for (int i = 0; i < numIntervals; i++) {
        result.FilledContours[i].LevelLow = bounds[i];
        result.FilledContours[i].LevelHigh = bounds[i + 1];
        result.FilledContours[i].ColorIndex = i;
        result.FilledContours[i].Polygons = AssemblePolygonsFromEdges(edgesByInterval[i]);
    }

    return result;
}
