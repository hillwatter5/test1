/**
 * ContourGenerator.cpp
 * 功能说明：等值线与填色多边形核心算法完整实现。
 *          算法对齐 matplotlib._contour 的 Marching Squares 逻辑。
 *          不依赖 MRing/MPolygon 上的任何扩展方法，
 *          环面积和包含判断均由本文件内静态函数独立实现。
 * 作者：ContourDev
 * 创建时间：2026-02-18
 */

#include "ContourGenerator.h"
#include <algorithm>

// ============================================================
//  Marching Squares 查找表
//
//  单元格角点编号（逆时针）：
//      c3(col, row+1) ---- c2(col+1, row+1)
//           |                    |
//      c0(col, row)   ---- c1(col+1, row)
//
//  边编号：e0=底(c0→c1) e1=右(c1→c2) e2=顶(c2→c3) e3=左(c3→c0)
//
//  caseIndex = (c0>=level?1:0) | (c1>=level?2:0)
//            | (c2>=level?4:0) | (c3>=level?8:0)
// ============================================================

struct CSegDef { int E1, E2; };

static const int SEG_COUNT[16] = {0,1,1,1,1,2,1,1,1,1,2,1,1,1,1,0};

static const int EC[4][2] = {{0,1},{1,2},{2,3},{3,0}};

static const CSegDef CTAB[16][2] = {
    {{-1,-1},{-1,-1}}, {{ 0, 3},{-1,-1}}, {{ 0, 1},{-1,-1}}, {{ 3, 1},{-1,-1}},
    {{ 1, 2},{-1,-1}}, {{ 0, 3},{ 1, 2}}, {{ 0, 2},{-1,-1}}, {{ 3, 2},{-1,-1}},
    {{ 2, 3},{-1,-1}}, {{ 0, 2},{-1,-1}}, {{ 0, 1},{ 2, 3}}, {{ 1, 2},{-1,-1}},
    {{ 1, 3},{-1,-1}}, {{ 0, 1},{-1,-1}}, {{ 0, 3},{-1,-1}}, {{-1,-1},{-1,-1}},
};

static const CSegDef SFLIP[2][2] = {
    {{ 0, 1},{ 2, 3}},
    {{ 0, 3},{ 1, 2}},
};

// ============================================================
//  独立几何辅助函数（不依赖 MRing/MPolygon 上的方法）
// ============================================================

/**
 * Shoelace 公式计算环的有符号面积。正=逆时针（外环），负=顺时针（内环）。
 */
double ComputeSignedArea(const std::vector<MPoint>& pts)
{
    int n = (int)pts.size();
    if (n < 3) return 0;
    double area = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += (double)pts[i].X * pts[j].Y;
        area -= (double)pts[j].X * pts[i].Y;
    }
    return area * 0.5;
}

/**
 * 射线法判断点是否在多边形环内
 */
static bool PointInRing(const std::vector<MPoint>& pts, const MPoint& pt)
{
    int n = (int)pts.size();
    if (n < 3) return false;
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        if (((pts[i].Y > pt.Y) != (pts[j].Y > pt.Y)) &&
            (pt.X < (pts[j].X - pts[i].X) * (pt.Y - pts[i].Y) /
             (pts[j].Y - pts[i].Y) + pts[i].X))
            inside = !inside;
    }
    return inside;
}

// ============================================================
//  插值和去重
// ============================================================

/**
 * 规范化线性插值：统一从坐标较小端向较大端方向插值，
 * 确保同一网格边从两个相邻单元格计算时结果完全一致。
 */
static MPoint CanonicalInterpolate(
    const MPoint& p1, float v1, const MPoint& p2, float v2, float level)
{
    if (fabs(v2 - v1) < (float)EPSILON)
        return MPoint((p1.X + p2.X) * 0.5f, (p1.Y + p2.Y) * 0.5f);

    const MPoint* pa = &p1; const MPoint* pb = &p2;
    float va = v1, vb = v2;
    bool sw = (p1.X > p2.X + (float)EPSILON) ||
              (fabs(p1.X - p2.X) < (float)EPSILON && p1.Y > p2.Y + (float)EPSILON);
    if (sw) { pa = &p2; pb = &p1; va = v2; vb = v1; }

    float t = (level - va) / (vb - va);
    t = std::max(0.0f, std::min(1.0f, t));
    if (t <= 0.0f) return *pa;
    if (t >= 1.0f) return *pb;
    return MPoint(pa->X + t * (pb->X - pa->X), pa->Y + t * (pb->Y - pa->Y));
}

static void Dedup(std::vector<MPoint>& pts)
{
    if (pts.size() < 2) return;
    std::vector<MPoint> out;
    out.reserve(pts.size());
    out.push_back(pts[0]);
    CPointKey prev(pts[0]);
    for (size_t i = 1; i < pts.size(); i++) {
        CPointKey cur(pts[i]);
        if (!(cur == prev)) { out.push_back(pts[i]); prev = cur; }
    }
    pts.swap(out);
}

// ============================================================
//  线段链接辅助
// ============================================================

static void ExtendFwd(
    std::vector<MPoint>& ch, const std::vector<CSegPair>& segs,
    std::map<CPointKey, std::vector<int>>& ep, std::vector<bool>& vis)
{
    bool ext = true;
    while (ext) {
        ext = false;
        CPointKey tk(ch.back());
        auto it = ep.find(tk);
        if (it == ep.end()) break;
        for (size_t k = 0; k < it->second.size(); k++) {
            int idx = it->second[k];
            if (vis[idx]) continue;
            CPointKey sk(segs[idx].first), ek(segs[idx].second);
            if (sk == tk) { vis[idx]=true; ch.push_back(segs[idx].second); ext=true; break; }
            if (ek == tk) { vis[idx]=true; ch.push_back(segs[idx].first);  ext=true; break; }
        }
    }
}

static void ExtendBwd(
    std::vector<MPoint>& ch, const std::vector<CSegPair>& segs,
    std::map<CPointKey, std::vector<int>>& ep, std::vector<bool>& vis)
{
    bool ext = true;
    while (ext) {
        ext = false;
        CPointKey hk(ch.front());
        auto it = ep.find(hk);
        if (it == ep.end()) break;
        for (size_t k = 0; k < it->second.size(); k++) {
            int idx = it->second[k];
            if (vis[idx]) continue;
            CPointKey sk(segs[idx].first), ek(segs[idx].second);
            if (sk == hk) { vis[idx]=true; ch.insert(ch.begin(),segs[idx].second); ext=true; break; }
            if (ek == hk) { vis[idx]=true; ch.insert(ch.begin(),segs[idx].first);  ext=true; break; }
        }
    }
}

static void MergeFragments(std::vector<std::vector<MPoint>>& chs)
{
    bool mg = true;
    while (mg) {
        mg = false;
        for (size_t i = 0; i < chs.size() && !mg; i++) {
            if (chs[i].empty()) continue;
            CPointKey ti(chs[i].back()), hi(chs[i].front());
            for (size_t j = i+1; j < chs.size() && !mg; j++) {
                if (chs[j].empty()) continue;
                CPointKey hj(chs[j].front()), tj(chs[j].back());
                if (ti == hj) {
                    chs[i].insert(chs[i].end(), chs[j].begin()+1, chs[j].end());
                    chs[j].clear(); mg=true;
                } else if (hi == tj) {
                    chs[j].insert(chs[j].end(), chs[i].begin()+1, chs[i].end());
                    chs[i].swap(chs[j]); chs[j].clear(); mg=true;
                } else if (ti == tj) {
                    std::reverse(chs[j].begin(), chs[j].end());
                    chs[i].insert(chs[i].end(), chs[j].begin()+1, chs[j].end());
                    chs[j].clear(); mg=true;
                } else if (hi == hj) {
                    std::reverse(chs[i].begin(), chs[i].end());
                    chs[i].insert(chs[i].end(), chs[j].begin()+1, chs[j].end());
                    chs[j].clear(); mg=true;
                }
            }
        }
    }
    std::vector<std::vector<MPoint>> r;
    for (auto& c : chs) if (!c.empty()) r.push_back(c);
    chs.swap(r);
}

// ============================================================
//  单元格级别：等值线段生成
// ============================================================

static void GenCellContour(
    const float v[4], const MPoint c[4], float level,
    int lvIdx, int nLev, const std::vector<CCachedIsect>* cache,
    std::vector<CSegPair>& out)
{
    int ci = 0;
    if (v[0] >= level) ci |= 1;
    if (v[1] >= level) ci |= 2;
    if (v[2] >= level) ci |= 4;
    if (v[3] >= level) ci |= 8;

    int ns = SEG_COUNT[ci];
    if (ns == 0) return;

    const CSegDef* sd = CTAB[ci];
    if (ci == 5 || ci == 10) {
        float ctr = (v[0]+v[1]+v[2]+v[3]) * 0.25f;
        if (ctr >= level) sd = (ci == 5) ? SFLIP[0] : SFLIP[1];
    }

    for (int s = 0; s < ns; s++) {
        int e1 = sd[s].E1, e2 = sd[s].E2;
        if (e1 < 0) continue;
        MPoint p1, p2;
        if (cache && lvIdx >= 0) {
            p1 = (*cache)[e1*nLev+lvIdx].Pt;
            p2 = (*cache)[e2*nLev+lvIdx].Pt;
        } else {
            p1 = CanonicalInterpolate(c[EC[e1][0]], v[EC[e1][0]], c[EC[e1][1]], v[EC[e1][1]], level);
            p2 = CanonicalInterpolate(c[EC[e2][0]], v[EC[e2][0]], c[EC[e2][1]], v[EC[e2][1]], level);
        }
        if (!(p1 == p2)) out.push_back({p1, p2});
    }
}

// ============================================================
//  单元格级别：填色多边形边生成
// ============================================================

static bool NeedSplit(const int st[4], int cs)
{
    bool d02=(st[0]==1&&st[2]==1), d13=(st[1]==1&&st[3]==1);
    bool o02=(st[0]!=1&&st[2]!=1), o13=(st[1]!=1&&st[3]!=1);
    return (d02&&o13&&cs!=1)||(d13&&o02&&cs!=1);
}

static void SplitSaddle(const std::vector<MPoint>& pv, std::vector<CSegPair>& out)
{
    int n=(int)pv.size();
    if (n<6) { for(int i=0;i<n;i++) out.push_back({pv[i],pv[(i+1)%n]}); return; }
    int h=n/2;
    for(int i=0;i<h;i++) out.push_back({pv[i],pv[(i+1)%h]});
    for(int i=h;i<n;i++) out.push_back({pv[i],pv[h+((i-h+1)%(n-h))]});
}

static void GenCellFill(
    const float v[4], const MPoint c[4],
    float lo, float hi, int loIdx, int hiIdx, int nLev,
    const std::vector<CCachedIsect>* cache,
    std::vector<MPoint>& buf, std::vector<CSegPair>& out)
{
    int st[4];
    for (int i=0;i<4;i++) { st[i]=(v[i]<lo)?0:(v[i]<hi)?1:2; }

    bool any=false, cross=false;
    for (int i=0;i<4;i++) {
        if (st[i]==1) any=true;
        int j=(i+1)%4;
        if ((v[i]<lo)!=(v[j]<lo)) cross=true;
        if ((v[i]<hi)!=(v[j]<hi)) cross=true;
    }
    if (!any && !cross) return;

    // 沿边界逆时针行走收集 between 区域顶点
    buf.clear();
    const int eS[4]={0,1,2,3}, eE[4]={1,2,3,0};
    for (int e=0;e<4;e++) {
        int cA=eS[e], cB=eE[e];
        float vA=v[cA], vB=v[cB];
        if (st[cA]==1) buf.push_back(c[cA]);

        struct Is { float T; MPoint P; };
        Is is[2]; int ni=0;
        if ((vA<lo)!=(vB<lo)) {
            float t=(lo-vA)/(vB-vA);
            MPoint p = (cache&&loIdx>=0) ? (*cache)[e*nLev+loIdx].Pt
                       : CanonicalInterpolate(c[cA],vA,c[cB],vB,lo);
            is[ni++] = {t, p};
        }
        if ((vA<hi)!=(vB<hi)) {
            float t=(hi-vA)/(vB-vA);
            MPoint p = (cache&&hiIdx>=0&&hiIdx<nLev) ? (*cache)[e*nLev+hiIdx].Pt
                       : CanonicalInterpolate(c[cA],vA,c[cB],vB,hi);
            is[ni++] = {t, p};
        }
        if (ni==2 && is[0].T>is[1].T) std::swap(is[0],is[1]);
        for (int k=0;k<ni;k++) buf.push_back(is[k].P);
    }

    Dedup(buf);
    if (buf.size()<3) return;

    float ctr=(v[0]+v[1]+v[2]+v[3])*0.25f;
    int cs=(ctr<lo)?0:(ctr<hi)?1:2;

    if (NeedSplit(st,cs) && buf.size()>=6) {
        SplitSaddle(buf, out);
    } else {
        for (size_t i=0;i<buf.size();i++)
            out.push_back({buf[i], buf[(i+1)%buf.size()]});
    }
}

// ============================================================
//  CContourGenerator 实现
// ============================================================

CContourGenerator::CContourGenerator() : _rows(0), _cols(0) {}
CContourGenerator::~CContourGenerator() {}

void CContourGenerator::Initialize(
    const float* xCoords, int cols,
    const float* yCoords, int rows,
    const float* values)
{
    _rows=rows; _cols=cols;
    _xCoords.assign(xCoords, xCoords+cols);
    _yCoords.assign(yCoords, yCoords+rows);
    _values.assign(values, values+rows*cols);
}

float CContourGenerator::GetValue(int row, int col) const
{
    return _values[row*_cols+col];
}

MPoint CContourGenerator::GetGridPoint(int row, int col) const
{
    return MPoint(_xCoords[col], _yCoords[row]);
}

// ============================================================
//  线段链接为 MCurve
// ============================================================

std::vector<MCurve> CContourGenerator::ChainSegments(
    std::vector<CSegPair>& segs) const
{
    int n=(int)segs.size();
    if (n==0) return {};

    std::map<CPointKey,std::vector<int>> ep;
    for (int i=0;i<n;i++) {
        ep[CPointKey(segs[i].first)].push_back(i);
        ep[CPointKey(segs[i].second)].push_back(i);
    }

    std::vector<bool> vis(n,false);
    std::vector<std::vector<MPoint>> chs;

    for (int i=0;i<n;i++) {
        if (vis[i]) continue;
        vis[i]=true;
        std::vector<MPoint> ch;
        ch.push_back(segs[i].first);
        ch.push_back(segs[i].second);
        ExtendFwd(ch, segs, ep, vis);
        ExtendBwd(ch, segs, ep, vis);
        Dedup(ch);
        if (ch.size()>=2) chs.push_back(ch);
    }
    MergeFragments(chs);

    std::vector<MCurve> curves;
    for (auto& ch : chs) {
        Dedup(ch);
        if (ch.size()>=2) {
            MCurve cv;
            cv.Points = ch;
            curves.push_back(cv);
        }
    }
    return curves;
}

// ============================================================
//  边对消、环追踪、归组为 MPolygon
// ============================================================

void CContourGenerator::CancelOpposingEdges(
    std::vector<CSegPair>& edges, std::vector<CSegPair>& rem) const
{
    std::map<CEdgeKey,int> cnt;
    std::map<CEdgeKey,CSegPair> dat;
    for (auto& e : edges) {
        CEdgeKey f(CPointKey(e.first),CPointKey(e.second));
        CEdgeKey r(CPointKey(e.second),CPointKey(e.first));
        if (cnt.count(r)&&cnt[r]>0) cnt[r]--;
        else { cnt[f]++; dat[f]=e; }
    }
    rem.clear();
    for (auto& kv : cnt)
        for (int i=0;i<kv.second;i++) rem.push_back(dat[kv.first]);
}

std::vector<MRing> CContourGenerator::TraceRings(
    std::vector<CSegPair>& edges) const
{
    std::vector<MRing> rings;
    if (edges.empty()) return rings;

    std::map<CPointKey,std::vector<std::pair<CPointKey,MPoint>>> adj;
    std::map<CPointKey,MPoint> k2p;
    for (auto& e : edges) {
        CPointKey ks(e.first), ke(e.second);
        adj[ks].push_back({ke, e.second});
        k2p[ks]=e.first; k2p[ke]=e.second;
    }

    int mx=(int)edges.size()+1;
    while (!adj.empty()) {
        CPointKey sk=adj.begin()->first;
        MRing ring;
        ring.Points.push_back(k2p[sk]);
        CPointKey cur=sk;
        bool ok=true;
        for (int step=0;step<mx;step++) {
            auto it=adj.find(cur);
            if (it==adj.end()||it->second.empty()) {ok=false; break;}
            auto nx=it->second.back();
            it->second.pop_back();
            if (it->second.empty()) adj.erase(it);
            if (nx.first==sk) break;
            ring.Points.push_back(nx.second);
            cur=nx.first;
        }
        if (ok && ring.Points.size()>=3) {
            Dedup(ring.Points);
            if (ring.Points.size()>=3) rings.push_back(ring);
        }
    }
    return rings;
}

std::vector<MPolygon> CContourGenerator::GroupRingsToPolygons(
    std::vector<MRing>& rings) const
{
    if (rings.empty()) return {};

    std::vector<int> outerIdx, innerIdx;
    for (size_t i=0;i<rings.size();i++) {
        double a = ComputeSignedArea(rings[i].Points);
        if (a > 0) outerIdx.push_back((int)i);
        else        innerIdx.push_back((int)i);
    }
    if (outerIdx.empty()) {
        for (size_t i=0;i<rings.size();i++) {
            std::reverse(rings[i].Points.begin(), rings[i].Points.end());
            outerIdx.push_back((int)i);
        }
        innerIdx.clear();
    }

    std::vector<MPolygon> pv(outerIdx.size());
    for (size_t i=0;i<outerIdx.size();i++)
        pv[i].OuterRings.push_back(rings[outerIdx[i]]);

    for (size_t i=0;i<innerIdx.size();i++) {
        MPoint tp = rings[innerIdx[i]].Points[0];
        int best=-1; double bestA=1e30;
        for (size_t j=0;j<outerIdx.size();j++) {
            if (PointInRing(rings[outerIdx[j]].Points, tp)) {
                double a=fabs(ComputeSignedArea(rings[outerIdx[j]].Points));
                if (a<bestA) { bestA=a; best=(int)j; }
            }
        }
        if (best>=0) pv[best].InnerRings.push_back(rings[innerIdx[i]]);
    }
    return pv;
}

std::vector<MPolygon> CContourGenerator::AssemblePolygonsFromEdges(
    std::vector<CSegPair>& edges) const
{
    std::vector<CSegPair> rem;
    CancelOpposingEdges(edges, rem);
    std::vector<MRing> rings = TraceRings(rem);
    return GroupRingsToPolygons(rings);
}

// ============================================================
//  公共接口：仅计算等值线
// ============================================================

std::vector<MContourLine> CContourGenerator::ComputeContourLines(
    const float* levels, int numLevels)
{
    std::vector<MContourLine> result;

    for (int lv=0; lv<numLevels; lv++) {
        std::vector<CSegPair> allSegs;
        for (int r=0; r<_rows-1; r++)
            for (int c=0; c<_cols-1; c++) {
                float v[4]={GetValue(r,c),GetValue(r,c+1),GetValue(r+1,c+1),GetValue(r+1,c)};
                MPoint cn[4]={GetGridPoint(r,c),GetGridPoint(r,c+1),GetGridPoint(r+1,c+1),GetGridPoint(r+1,c)};
                GenCellContour(v, cn, levels[lv], -1, 0, nullptr, allSegs);
            }
        std::vector<MCurve> curves = ChainSegments(allSegs);
        for (auto& cv : curves) {
            MContourLine cl;
            cl.Level = levels[lv];
            cl.Curve = cv;
            result.push_back(cl);
        }
    }
    return result;
}

// ============================================================
//  公共接口：仅计算填色
// ============================================================

std::vector<MContourFill> CContourGenerator::ComputeFilledContours(
    const float* levels, int numLevels)
{
    if (numLevels<1) return {};
    int nI = numLevels + 1;
    std::vector<MContourFill> result;

    for (int i=0; i<nI; i++) {
        float lo = (i==0) ? -FLT_MAX : levels[i-1];
        float hi = (i==numLevels) ? FLT_MAX : levels[i];

        std::vector<CSegPair> allEdges;
        std::vector<MPoint> buf;
        for (int r=0; r<_rows-1; r++)
            for (int c=0; c<_cols-1; c++) {
                float v[4]={GetValue(r,c),GetValue(r,c+1),GetValue(r+1,c+1),GetValue(r+1,c)};
                MPoint cn[4]={GetGridPoint(r,c),GetGridPoint(r,c+1),GetGridPoint(r+1,c+1),GetGridPoint(r+1,c)};
                GenCellFill(v, cn, lo, hi, -1, -1, 0, nullptr, buf, allEdges);
            }

        std::vector<MPolygon> polys = AssemblePolygonsFromEdges(allEdges);

        float lMin, lMax;
        if (i == 0)           { lMin = levels[0];          lMax = levels[0]; }
        else if (i==numLevels){ lMin = levels[numLevels-1]; lMax = levels[numLevels-1]; }
        else                  { lMin = levels[i-1];         lMax = levels[i]; }

        for (auto& p : polys) {
            MContourFill cf;
            cf.Polygon = p;
            cf.LevelMin = lMin;
            cf.LevelMax = lMax;
            cf.ColorIndex = i;
            result.push_back(cf);
        }
    }
    return result;
}

// ============================================================
//  公共接口：合并计算（性能优化核心）
// ============================================================

CMergedResult CContourGenerator::ComputeContourAndFill(
    const float* levels, int numLevels)
{
    CMergedResult result;
    if (numLevels < 1) return result;

    int nI = numLevels + 1;

    // 扩展边界：[-FLT_MAX, levels[0..N-1], FLT_MAX]
    std::vector<float> bounds(nI + 1);
    bounds[0] = -FLT_MAX;
    for (int i=0;i<numLevels;i++) bounds[i+1] = levels[i];
    bounds[nI] = FLT_MAX;

    // 每个 level / 区间的线段和边
    std::vector<std::vector<CSegPair>> segByLv(numLevels);
    std::vector<std::vector<CSegPair>> edgeByI(nI);

    // 可复用缓冲
    std::vector<CCachedIsect> cache(4 * numLevels);
    std::vector<MPoint> buf;
    buf.reserve(16);

    // ---- 单次网格遍历 ----
    for (int row=0; row<_rows-1; row++) {
        for (int col=0; col<_cols-1; col++) {
            float v[4] = {
                GetValue(row,col), GetValue(row,col+1),
                GetValue(row+1,col+1), GetValue(row+1,col)
            };
            MPoint cn[4] = {
                GetGridPoint(row,col), GetGridPoint(row,col+1),
                GetGridPoint(row+1,col+1), GetGridPoint(row+1,col)
            };

            // 预计算 4边×N级别 交点缓存
            for (int e=0;e<4;e++) {
                int cA=EC[e][0], cB=EC[e][1];
                for (int lv=0;lv<numLevels;lv++) {
                    int idx=e*numLevels+lv;
                    cache[idx].Has = (v[cA]<levels[lv])!=(v[cB]<levels[lv]);
                    if (cache[idx].Has)
                        cache[idx].Pt = CanonicalInterpolate(
                            cn[cA], v[cA], cn[cB], v[cB], levels[lv]);
                }
            }

            // 等值线段
            for (int lv=0;lv<numLevels;lv++)
                GenCellContour(v, cn, levels[lv], lv, numLevels, &cache, segByLv[lv]);

            // 填色有向边
            for (int i=0;i<nI;i++) {
                int loIdx=i-1, hiIdx=(i<numLevels)?i:-1;
                GenCellFill(v, cn, bounds[i], bounds[i+1],
                            loIdx, hiIdx, numLevels, &cache, buf, edgeByI[i]);
            }
        }
    }

    // ---- 后处理：等值线 ----
    for (int lv=0;lv<numLevels;lv++) {
        std::vector<MCurve> curves = ChainSegments(segByLv[lv]);
        for (auto& cv : curves) {
            MContourLine cl;
            cl.Level = levels[lv];
            cl.Curve = cv;
            result.ContourLines.push_back(cl);
        }
    }

    // ---- 后处理：填色 ----
    for (int i=0;i<nI;i++) {
        std::vector<MPolygon> polys = AssemblePolygonsFromEdges(edgeByI[i]);

        float lMin, lMax;
        if (i==0)           { lMin=levels[0];          lMax=levels[0]; }
        else if (i==numLevels){ lMin=levels[numLevels-1]; lMax=levels[numLevels-1]; }
        else                  { lMin=levels[i-1];         lMax=levels[i]; }

        for (auto& p : polys) {
            MContourFill cf;
            cf.Polygon = p;
            cf.LevelMin = lMin;
            cf.LevelMax = lMax;
            cf.ColorIndex = i;
            result.FilledContours.push_back(cf);
        }
    }

    return result;
}
