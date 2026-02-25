/**
 * main.cpp
 * 功能说明：等值线与填色算法测试。验证 MContourLine/MContourFill 输出格式。
 * 作者：ContourDev
 * 创建时间：2026-02-18
 */

#include "ContourGenerator.h"
#include <iostream>
#include <iomanip>

static void BuildTestGrid(
    int n, float lo, float hi,
    std::vector<float>& x, std::vector<float>& y, std::vector<float>& v)
{
    float s = (hi - lo) / (n - 1);
    x.resize(n); y.resize(n); v.resize(n * n);
    for (int i = 0; i < n; i++) { x[i] = lo + i * s; y[i] = lo + i * s; }
    for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++)
            v[r * n + c] = x[c] * x[c] + y[r] * y[r];
}

int main()
{
    const int N = 31;
    std::vector<float> xc, yc, vals;
    BuildTestGrid(N, -2.5f, 2.5f, xc, yc, vals);

    CContourGenerator gen;
    gen.Initialize(xc.data(), N, yc.data(), N, vals.data());

    // ---- 测试: levels = [1, 2.5, 3] ----
    std::cout << "==========================================" << std::endl;
    std::cout << "  f(x,y) = x^2 + y^2,  levels = [1, 2.5, 3]" << std::endl;
    std::cout << "  期望: 3 条等值线, 4 个填色区间 (ColorIndex 0~3)" << std::endl;
    std::cout << "==========================================" << std::endl;

    float levels[] = {1.0f, 2.5f, 3.0f};
    CMergedResult res = gen.ComputeContourAndFill(levels, 3);

    std::cout << "\n--- MContourLine 列表 (" << res.ContourLines.size() << " 条) ---" << std::endl;
    for (size_t i = 0; i < res.ContourLines.size(); i++) {
        const MContourLine& cl = res.ContourLines[i];
        std::cout << std::fixed << std::setprecision(2)
                  << "  [" << i << "] Level=" << cl.Level
                  << "  点数=" << cl.Curve.Points.size() << std::endl;
    }

    std::cout << "\n--- MContourFill 列表 (" << res.FilledContours.size() << " 条) ---" << std::endl;
    for (size_t i = 0; i < res.FilledContours.size(); i++) {
        const MContourFill& cf = res.FilledContours[i];
        std::cout << std::fixed << std::setprecision(2)
                  << "  [" << i << "] ColorIndex=" << cf.ColorIndex
                  << "  LevelMin=" << cf.LevelMin << "  LevelMax=" << cf.LevelMax
                  << "  外环=" << cf.Polygon.OuterRings.size()
                  << "  内环=" << cf.Polygon.InnerRings.size();

        double netArea = 0;
        for (auto& r : cf.Polygon.OuterRings)
            netArea += ComputeSignedArea(r.Points);
        for (auto& r : cf.Polygon.InnerRings)
            netArea += ComputeSignedArea(r.Points);
        std::cout << "  净面积=" << std::setprecision(3) << netArea << std::endl;
    }

    std::cout << "\n--- 理论对比 ---" << std::endl;
    std::cout << "  ColorIndex=0 (v<1)       : 面积≈π       = 3.142" << std::endl;
    std::cout << "  ColorIndex=1 (1≤v<2.5)   : 面积≈π*1.5   = 4.712" << std::endl;
    std::cout << "  ColorIndex=2 (2.5≤v<3)   : 面积≈π*0.5   = 1.571" << std::endl;
    std::cout << "  ColorIndex=3 (v≥3)       : 面积≈25-3π   = 15.575" << std::endl;

    std::cout << "\n==========================================" << std::endl;
    std::cout << "  测试完成" << std::endl;
    std::cout << "==========================================" << std::endl;
    return 0;
}
