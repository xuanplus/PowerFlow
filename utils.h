#pragma once
#include "complex.h"

// ֧·
struct Branch
{
    int node1{};  // ��ʼ�ڵ�
    int node2{};  // ��ֹ�ڵ�
    complex Y{0, 0};  // ֧·����
};


// �ڵ�
struct Node
{
    int type;       // �ڵ����ͣ�1-PQ, 2-PV, 3-ƽ��
    float U;        // ��ѹ
    float a;        // ���
    float e;        // ��ѹʵ��
    float f;        // ��ѹ�鲿
    float P;        // �й�����
    float Q;        // �޹�����
};
