#pragma once
#include "complex.h"

// 支路
struct Branch
{
    int node1{};  // 起始节点
    int node2{};  // 终止节点
    complex Y{0, 0};  // 支路导纳
};


// 节点
struct Node
{
    int type;       // 节点类型：1-PQ, 2-PV, 3-平衡
    float U;        // 电压
    float a;        // 相角
    float e;        // 电压实部
    float f;        // 电压虚部
    float P;        // 有功功率
    float Q;        // 无功功率
};
