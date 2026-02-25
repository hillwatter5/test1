# Makefile
# 功能说明：编译等值线与填色算法测试程序
# 作者：ContourDev
# 创建时间：2026-02-18

CXX      = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2
TARGET   = contour_test
SRCS     = main.cpp ContourGenerator.cpp
OBJS     = $(SRCS:.cpp=.o)

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
