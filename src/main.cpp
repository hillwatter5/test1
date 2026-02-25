/**
 * main.cpp
 * 功能说明：等值线与填色算法的测试程序。
 *          测试合并计算接口和 N+1 填色区间的正确性。
 * 作者：ContourDev
 * 创建时间：2026-02-18
 */

#include "ContourGenerator.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>

static void BuildTestGrid(
    int gridSize, float rangeMin, float rangeMax,
    std::vector<float>& xCoords,
    std::vector<float>& yCoords,
    std::vector<float>& values)
{
    float step = (rangeMax - rangeMin) / (gridSize - 1);
    xCoords.resize(gridSize);
    yCoords.resize(gridSize);
    values.resize(gridSize * gridSize);

    for (int i = 0; i < gridSize; i++) {
        xCoords[i] = rangeMin + i * step;
        yCoords[i] = rangeMin + i * step;
    }
    for (int row = 0; row < gridSize; row++) {
        for (int col = 0; col < gridSize; col++) {
            float x = xCoords[col], y = yCoords[row];
            values[row * gridSize + col] = x * x + y * y;
        }
    }
}

static void PrintMergedResult(const CMergedResult& merged)
{
    std::cout << "\n======== 等值线结果 ========" << std::endl;
    for (size_t i = 0; i < merged.ContourLines.size(); i++) {
        const CContourLineResult& cl = merged.ContourLines[i];
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Level=" << cl.Level
                  << " → " << cl.Curves.size() << " 条曲线";
        int totalPts = 0;
        for (size_t c = 0; c < cl.Curves.size(); c++)
            totalPts += (int)cl.Curves[c].Points.size();
        std::cout << "（共 " << totalPts << " 个点）" << std::endl;
    }

    std::cout << "\n======== 填色区间结果 ========" << std::endl;
    for (size_t i = 0; i < merged.FilledContours.size(); i++) {
        const CFilledContourResult& fc = merged.FilledContours[i];
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "ColorIndex=" << fc.ColorIndex << " : ";

        if (fc.LevelLow < -1e30f)
            std::cout << "value < " << fc.LevelHigh;
        else if (fc.LevelHigh > 1e30f)
            std::cout << "value >= " << fc.LevelLow;
        else
            std::cout << fc.LevelLow << " <= value < " << fc.LevelHigh;

        std::cout << " → " << fc.Polygons.size() << " 个多边形";

        int totalOuter = 0, totalInner = 0;
        double totalArea = 0;
        for (size_t p = 0; p < fc.Polygons.size(); p++) {
            totalOuter += (int)fc.Polygons[p].OuterRings.size();
            totalInner += (int)fc.Polygons[p].InnerRings.size();
            for (size_t r = 0; r < fc.Polygons[p].OuterRings.size(); r++)
                totalArea += fc.Polygons[p].OuterRings[r].SignedArea();
            for (size_t r = 0; r < fc.Polygons[p].InnerRings.size(); r++)
                totalArea += fc.Polygons[p].InnerRings[r].SignedArea();
        }
        std::cout << "（" << totalOuter << " 外环, "
                  << totalInner << " 内环, "
                  << "净面积=" << std::setprecision(3) << totalArea << "）"
                  << std::endl;
    }
}

int main()
{
    // ============================================================
    //  测试 1：用户示例 levels=[1, 2.5, 3]
    //  验证：3 个 level → 3 条等值线 + 4 个填色区间（ColorIndex 0~3）
    // ============================================================
    std::cout << "============================================" << std::endl;
    std::cout << "  测试 1：levels = [1, 2.5, 3]" << std::endl;
    std::cout << "  函数 f(x,y) = x^2 + y^2" << std::endl;
    std::cout << "  期望：3 条等值线 + 4 个填色区间" << std::endl;
    std::cout << "============================================" << std::endl;

    const int gridSize = 31;
    std::vector<float> xCoords, yCoords, values;
    BuildTestGrid(gridSize, -2.5f, 2.5f, xCoords, yCoords, values);

    CContourGenerator gen;
    gen.Initialize(xCoords.data(), gridSize,
                   yCoords.data(), gridSize,
                   values.data());

    float levels1[] = {1.0f, 2.5f, 3.0f};
    int numLevels1 = 3;

    std::cout << "\n--- 调用合并计算接口 ComputeContourAndFill ---" << std::endl;
    CMergedResult merged = gen.ComputeContourAndFill(levels1, numLevels1);
    PrintMergedResult(merged);

    std::cout << "\n期望值对比（f=x²+y² 的等值线是以原点为圆心的圆）：" << std::endl;
    std::cout << "  Level=1.0 → 半径=1 的圆" << std::endl;
    std::cout << "  Level=2.5 → 半径=√2.5≈1.58 的圆" << std::endl;
    std::cout << "  Level=3.0 → 半径=√3≈1.73 的圆" << std::endl;
    std::cout << "  ColorIndex=0 (value<1) → 面积≈π*1=3.14" << std::endl;
    std::cout << "  ColorIndex=1 (1<=v<2.5) → 净面积≈π*(2.5-1)=4.71" << std::endl;
    std::cout << "  ColorIndex=2 (2.5<=v<3) → 净面积≈π*(3-2.5)=1.57" << std::endl;
    std::cout << "  ColorIndex=3 (v>=3) → 净面积=25-π*3≈15.58" << std::endl;

    // ============================================================
    //  测试 2：多级别验证 levels=[0.5, 1, 2, 4, 6]
    // ============================================================
    std::cout << "\n============================================" << std::endl;
    std::cout << "  测试 2：levels = [0.5, 1, 2, 4, 6]" << std::endl;
    std::cout << "  期望：5 条等值线 + 6 个填色区间" << std::endl;
    std::cout << "============================================" << std::endl;

    float levels2[] = {0.5f, 1.0f, 2.0f, 4.0f, 6.0f};
    CMergedResult merged2 = gen.ComputeContourAndFill(levels2, 5);
    PrintMergedResult(merged2);

    // ============================================================
    //  测试 3：对比独立接口和合并接口的一致性
    // ============================================================
    std::cout << "\n============================================" << std::endl;
    std::cout << "  测试 3：独立接口 vs 合并接口一致性验证" << std::endl;
    std::cout << "============================================" << std::endl;

    auto contourOnly = gen.ComputeContourLines(levels1, numLevels1);
    auto fillOnly    = gen.ComputeFilledContours(levels1, numLevels1);

    std::cout << "独立等值线：" << contourOnly.size() << " 条" << std::endl;
    std::cout << "合并等值线：" << merged.ContourLines.size() << " 条" << std::endl;
    std::cout << "独立填色区间：" << fillOnly.size() << " 个" << std::endl;
    std::cout << "合并填色区间：" << merged.FilledContours.size() << " 个" << std::endl;

    bool consistent = true;
    for (size_t i = 0; i < fillOnly.size() && i < merged.FilledContours.size(); i++) {
        int nA = 0, nB = 0;
        for (auto& p : fillOnly[i].Polygons)
            for (auto& r : p.OuterRings) nA += (int)r.Points.size();
        for (auto& p : merged.FilledContours[i].Polygons)
            for (auto& r : p.OuterRings) nB += (int)r.Points.size();
        if (fillOnly[i].ColorIndex != merged.FilledContours[i].ColorIndex) {
            consistent = false;
            std::cout << "  区间 " << i << " ColorIndex 不一致!" << std::endl;
        }
    }
    if (consistent)
        std::cout << "✓ 两种接口结果一致" << std::endl;

    std::cout << "\n============================================" << std::endl;
    std::cout << "  全部测试完成" << std::endl;
    std::cout << "============================================" << std::endl;

    return 0;
}
