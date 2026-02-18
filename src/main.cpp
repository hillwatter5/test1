/**
 * main.cpp
 * 功能说明：等值线与填色算法的测试程序。
 *          构造 X²+Y² 的二维网格数据，调用核心函数并输出结果。
 * 作者：ContourDev
 * 创建时间：2026-02-18
 */

#include "ContourGenerator.h"
#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * 构造测试网格数据：f(x,y) = x² + y²
 * 网格范围 [-2, 2] × [-2, 2]，分辨率由 gridSize 控制
 */
static void BuildTestGrid(
    int gridSize,
    std::vector<float>& xCoords,
    std::vector<float>& yCoords,
    std::vector<float>& values)
{
    float minVal = -2.0f, maxVal = 2.0f;
    float step = (maxVal - minVal) / (gridSize - 1);

    xCoords.resize(gridSize);
    yCoords.resize(gridSize);
    values.resize(gridSize * gridSize);

    for (int i = 0; i < gridSize; i++) {
        xCoords[i] = minVal + i * step;
        yCoords[i] = minVal + i * step;
    }

    for (int row = 0; row < gridSize; row++) {
        for (int col = 0; col < gridSize; col++) {
            float x = xCoords[col];
            float y = yCoords[row];
            values[row * gridSize + col] = x * x + y * y;
        }
    }
}

/**
 * 打印等值线计算结果
 */
static void PrintContourLines(
    const std::vector<CContourLineResult>& results)
{
    std::cout << "========================================" << std::endl;
    std::cout << "  等值线计算结果 (Contour Lines)" << std::endl;
    std::cout << "========================================" << std::endl;

    for (size_t i = 0; i < results.size(); i++) {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "\n等值级别 Level = " << results[i].Level
                  << " ，共 " << results[i].Curves.size() << " 条曲线"
                  << std::endl;

        for (size_t c = 0; c < results[i].Curves.size(); c++) {
            const MCurve& curve = results[i].Curves[c];
            std::cout << "  曲线 #" << c
                      << "（" << curve.Points.size() << " 个点）:";

            int showCount = std::min((int)curve.Points.size(), 6);
            for (int p = 0; p < showCount; p++) {
                std::cout << " (" << std::setprecision(3)
                          << curve.Points[p].X << ","
                          << curve.Points[p].Y << ")";
            }
            if ((int)curve.Points.size() > showCount)
                std::cout << " ...";
            std::cout << std::endl;
        }
    }
}

/**
 * 打印填色多边形计算结果
 */
static void PrintFilledContours(
    const std::vector<CFilledContourResult>& results)
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "  填色多边形计算结果 (Filled Contours)" << std::endl;
    std::cout << "========================================" << std::endl;

    for (size_t i = 0; i < results.size(); i++) {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "\n区间 [" << results[i].LevelLow
                  << ", " << results[i].LevelHigh << ")"
                  << " ，共 " << results[i].Polygons.size() << " 个多边形"
                  << std::endl;

        for (size_t p = 0; p < results[i].Polygons.size(); p++) {
            const MPolygon& poly = results[i].Polygons[p];
            std::cout << "  多边形 #" << p << "："
                      << poly.OuterRings.size() << " 个外环，"
                      << poly.InnerRings.size() << " 个内环"
                      << std::endl;

            for (size_t r = 0; r < poly.OuterRings.size(); r++) {
                const MRing& ring = poly.OuterRings[r];
                std::cout << "    外环 #" << r
                          << "（" << ring.Points.size() << " 个点，面积="
                          << std::setprecision(4)
                          << ring.SignedArea() << "）:";

                int showCount = std::min((int)ring.Points.size(), 5);
                for (int k = 0; k < showCount; k++) {
                    std::cout << " (" << std::setprecision(3)
                              << ring.Points[k].X << ","
                              << ring.Points[k].Y << ")";
                }
                if ((int)ring.Points.size() > showCount)
                    std::cout << " ...";
                std::cout << std::endl;
            }

            for (size_t r = 0; r < poly.InnerRings.size(); r++) {
                const MRing& ring = poly.InnerRings[r];
                std::cout << "    内环 #" << r
                          << "（" << ring.Points.size() << " 个点，面积="
                          << std::setprecision(4)
                          << ring.SignedArea() << "）:";

                int showCount = std::min((int)ring.Points.size(), 5);
                for (int k = 0; k < showCount; k++) {
                    std::cout << " (" << std::setprecision(3)
                              << ring.Points[k].X << ","
                              << ring.Points[k].Y << ")";
                }
                if ((int)ring.Points.size() > showCount)
                    std::cout << " ...";
                std::cout << std::endl;
            }
        }
    }
}

int main()
{
    std::cout << "============================================" << std::endl;
    std::cout << "  C++ 等值线与填色算法测试" << std::endl;
    std::cout << "  测试函数：f(x,y) = x^2 + y^2" << std::endl;
    std::cout << "  网格范围：[-2, 2] x [-2, 2]" << std::endl;
    std::cout << "============================================" << std::endl;

    // 构造 21×21 网格
    const int gridSize = 21;
    std::vector<float> xCoords, yCoords, values;
    BuildTestGrid(gridSize, xCoords, yCoords, values);

    std::cout << "\n网格大小：" << gridSize << " x " << gridSize << std::endl;
    std::cout << "值域范围：[0, 8]" << std::endl;

    // 初始化等值线生成器
    CContourGenerator generator;
    generator.Initialize(xCoords.data(), gridSize,
                         yCoords.data(), gridSize,
                         values.data());

    // ---- 测试 1：等值线计算 ----
    float contourLevels[] = {0.5f, 1.0f, 2.0f, 4.0f, 6.0f};
    int numContourLevels = sizeof(contourLevels) / sizeof(float);

    std::cout << "\n计算等值线，级别：";
    for (int i = 0; i < numContourLevels; i++)
        std::cout << contourLevels[i] << " ";
    std::cout << std::endl;

    std::vector<CContourLineResult> contourResults =
        generator.ComputeContourLines(contourLevels, numContourLevels);
    PrintContourLines(contourResults);

    // ---- 测试 2：填色多边形计算 ----
    float fillLevels[] = {0.0f, 1.0f, 2.0f, 4.0f, 8.0f};
    int numFillLevels = sizeof(fillLevels) / sizeof(float);

    std::cout << "\n计算填色多边形，区间级别：";
    for (int i = 0; i < numFillLevels; i++)
        std::cout << fillLevels[i] << " ";
    std::cout << std::endl;

    std::vector<CFilledContourResult> fillResults =
        generator.ComputeFilledContours(fillLevels, numFillLevels);
    PrintFilledContours(fillResults);

    std::cout << "\n============================================" << std::endl;
    std::cout << "  测试完成" << std::endl;
    std::cout << "============================================" << std::endl;

    return 0;
}
