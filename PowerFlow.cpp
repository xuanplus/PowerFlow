#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>

#include "complex.h"
#include "file.h"
#include "utils.h"

using namespace std;


/*
 * 输入文件格式：
 * 节点数
 * 支路数
 * -支路相连的节点1-节点2 |
 * -支路阻抗             |   <- 此两项循环支路数次
 * *节点类型    |
 * *电压幅值    |
 * *电压相角    |    <- 此五项循环节点数次，未知的项设为0
 * *有功功率    |
 * *无功功率    |    
 * 计算精度
*/


int main()
{
	// ----读取原始数据----

	// 读取文件
	ifstream inputFile("input.txt");
	// 输出文件
	ofstream outfile("output.txt");

	// 节点数、支路数
	const int n = read_data<int>(inputFile);
	const int m = read_data<int>(inputFile);

	outfile << "***原始数据***\n" << endl;
	outfile << "节点数：" << n << endl;
	outfile << "支路数：" << m << endl;

	// 支路向量
	vector<Branch> branches(m);

	// 节点向量
	vector<Node> nodes(n);

	// 导纳矩阵
	vector<vector<complex>> Y(n, vector<complex>(n, complex(0, 0)));

	// 填充支路向量
	outfile << "\n支路数据：" << endl;
	for (int i = 0; i < m; i++)
	{
		const pair<int, int> node = read_node(inputFile);
		const float real = read_data<float>(inputFile);
		const float imag = read_data<float>(inputFile);

		branches[i].node1 = node.first;
		branches[i].node2 = node.second;
		branches[i].Y = complex(1, 0) / complex(real, imag);
		outfile << "支路" << i + 1 << "：" << "两端节点为" << node.first << "和" << node.second << "，支路导纳 Y=" << branches[i].Y << endl;
	}

	// 填充导纳矩阵
	for (const Branch branch : branches)
	{
		const int i = branch.node1 - 1;
		const int j = branch.node2 - 1;
		const complex y = branch.Y;

		// 对角元素
		Y[i][i] += y;
		Y[j][j] += y;

		// 非对角元素
		Y[i][j] -= y;
		Y[j][i] -= y;
	}

	// 填充节点向量
	outfile << "\n节点数据：" << endl;
	for (int i = 0; i < n; i++)
	{
		// 由文件读取
		nodes[i].type = read_data<int>(inputFile);	// 节点类型
		nodes[i].U = read_data<float>(inputFile);	// 节点电压
		nodes[i].a = read_data<float>(inputFile);	// 电压相角
		nodes[i].P = read_data<float>(inputFile);	// 有功功率
		nodes[i].Q = read_data<float>(inputFile);	// 无功功率

		string type;
		switch (nodes[i].type)
		{
		case 1:
			type = "PQ节点";
			break;
		case 2:
			type = "PV节点";
			break;
		case 3:
			type = "平衡节点";
			break;
		default:
			break;
		}
		outfile << "节点" << i + 1 << "：" << type << "，电压：" << nodes[i].U << "，相角：" << nodes[i].a << "，注入功率 S=" << complex(nodes[i].P, nodes[i].Q) << endl;
	}

	// 输出导纳矩阵
	outfile << "\n导纳矩阵：" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			outfile << right << setw(19) << Y[i][j];
		}
		outfile << endl;
	}

	outfile << "\n***迭代过程***" << endl;

	// 收敛误差
	const float epsilon = read_data<float>(inputFile);

	inputFile.close();

	// ----设置电压初值----

	// 除平衡节点，设电压实部为1.0，虚部为0.0
	for (int i = 0; i < n; i++)
	{
		if (nodes[i].type != 3)
		{
			nodes[i].e = 1.0;
			nodes[i].f = 0.0;
		}
		else
		{
			nodes[i].e = nodes[i].U * cosf(nodes[i].a);
			nodes[i].f = nodes[i].U * sinf(nodes[i].a);
		}
	}

	// ----迭代求解潮流方程----

	// 最大迭代次数
	const int maxIter = 100;

	// 迭代次数
	int iter = 0;

	// 收敛标志
	bool converged = false;

	// 迭代求解
	while (iter < maxIter)
	{
		outfile << "\n第" << iter + 1 << "次迭代：" << endl;

		// 计算节点注入功率
		vector<float> P(n, 0);
		vector<float> Q(n, 0);
		vector<float> U(n, 0);

		for (int i = 0; i < n; i++)
		{
			// 平衡节点不计算
			if (nodes[i].type == 3) continue;

			// PQ节点和PV节点共有的P
			for (int j = 0; j < n; j++)
			{
				P[i] += nodes[i].e * (Y[i][j].real * nodes[j].e - Y[i][j].imag * nodes[j].f);
				P[i] += nodes[i].f * (Y[i][j].real * nodes[j].f + Y[i][j].imag * nodes[j].e);
			}

			// 分别计算PQ节点的Q和PV节点的U
			if (nodes[i].type == 1)
			{
				for (int j = 0; j < n; j++)
				{
					Q[i] += nodes[i].f * (Y[i][j].real * nodes[j].e - Y[i][j].imag * nodes[j].f);
					Q[i] -= nodes[i].e * (Y[i][j].real * nodes[j].f + Y[i][j].imag * nodes[j].e);
				}
			}
			else if (nodes[i].type == 2)
			{
				U[i] = sqrtf(powf(nodes[i].e, 2) + powf(nodes[i].f, 2));
			}
		}

		// 计算不平衡量
		vector<float> deltaP(n, 0);
		vector<float> deltaQ(n, 0);
		vector<float> deltaU(n, 0);

		for (int i = 0; i < n; i++)
		{
			// 平衡节点不计算
			if (nodes[i].type == 3) continue;

			deltaP[i] = nodes[i].P - P[i];
			deltaQ[i] = nodes[i].Q - Q[i];
			deltaU[i] = sqrtf(powf(nodes[i].U, 2) - powf(U[i], 2));
		}

		// 形成雅可比矩阵
		vector<vector<float>> H(n, vector<float>(n, 0));
		vector<vector<float>> N(n, vector<float>(n, 0));
		vector<vector<float>> J(n, vector<float>(n, 0));
		vector<vector<float>> L(n, vector<float>(n, 0));
		vector<vector<float>> R(n, vector<float>(n, 0));
		vector<vector<float>> S(n, vector<float>(n, 0));

		int i = 0;
		int j = 0;

		// 计算雅可比矩阵
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (nodes[i].type == 3 || nodes[j].type == 3) continue;
				if (i != j)
				{
					H[i][j] = -1.0 * Y[i][j].imag * nodes[i].e + Y[i][j].real * nodes[i].f;
					N[i][j] = Y[i][j].real * nodes[i].e + Y[i][j].imag * nodes[i].f;
					J[i][j] = -1 * N[i][j];
					L[i][j] = H[i][j];
				}
				else
				{
					// 计算节点注入电流
					complex I = Y[i][i] * complex(nodes[i].e, nodes[i].f);

					for (int k = 0; k < n; k++)
					{
						if (k != i)
						{
							I += Y[i][k] * complex(nodes[k].e, nodes[k].f);
						}
					}

					const float a = I.real;
					const float b = I.imag;

					H[i][j] = -1.0 * Y[i][j].imag * nodes[i].e + Y[i][j].real * nodes[i].f + b;
					N[i][j] = Y[i][j].real * nodes[i].e + Y[i][j].imag * nodes[i].f + a;
					J[i][j] = N[i][j] * -1 + a * 2;
					L[i][j] = H[i][j] - b * 2;
					R[i][j] = 2 * nodes[i].f;
					S[i][j] = 2 * nodes[i].e;
				}
			}
		}

		// 填充雅可比矩阵
		vector<vector<float>> JAC(2 * n - 2, vector<float>(2 * n - 2, 0));

		for (int i = 0; i < n; i++)
		{
			// 是否经过平衡节点
			bool sss = false;

			for (int j = 0; j < n; j++)
			{
				if (nodes[i].type == 3 || nodes[j].type == 3)
				{
					sss = true;
					continue;
				}
				int a = sss ? (i - 1) * 2 : i * 2;
				int b = sss ? (j - 1) * 2 : j * 2;
				JAC[a][b] = H[i][j];
				JAC[a][b + 1] = N[i][j];
				if (nodes[i].type == 1)
				{
					JAC[a + 1][b] = J[i][j];
					JAC[a + 1][b + 1] = L[i][j];
				}
				else
				{
					JAC[a + 1][b] = R[i][j];
					JAC[a + 1][b + 1] = S[i][j];
				}
			}
		}

		// 输出雅可比矩阵
		outfile << "\n雅可比矩阵：" << endl;
		for (int i = 0; i < 2 * n - 2; i++)
		{
			for (int j = 0; j < 2 * n - 2; j++)
			{
				outfile << right << setw(12) << JAC[i][j];
			}
			outfile << endl;
		}

		// 计算雅可比矩阵的逆矩阵

		int a = 2 * n - 2;
		vector<vector<float>> augmentedMatrix(a, vector<float>(2 * a, 0));

		// 创建增广矩阵
		for (int i = 0; i < a; ++i) {
			for (int j = 0; j < a; ++j) {
				augmentedMatrix[i][j] = JAC[i][j];
			}
			augmentedMatrix[i][i + a] = 1.0f;
		}

		// 高斯-约旦消元法
		for (int i = 0; i < a; ++i) {
			// 选主元
			float pivot = augmentedMatrix[i][i];
			for (int j = 0; j < 2 * a; ++j) {
				augmentedMatrix[i][j] /= pivot;
			}

			// 消元
			for (int k = 0; k < a; ++k) {
				if (k != i) {
					float factor = augmentedMatrix[k][i];
					for (int j = 0; j < 2 * a; ++j) {
						augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
					}
				}
			}
		}

		// 提取逆矩阵
		vector<vector<float>> inverseJAC(a, vector<float>(a, 0));
		for (int i = 0; i < a; ++i) {
			for (int j = 0; j < a; ++j) {
				inverseJAC[i][j] = augmentedMatrix[i][j + a];
			}
		}

		// 输出逆雅可比矩阵
		outfile << "\n逆雅可比矩阵：" << endl;
		for (int i = 0; i < a; i++)
		{
			for (int j = 0; j < a; j++)
			{
				outfile << right << setw(12) << inverseJAC[i][j];
			}
			outfile << endl;
		}

		// 计算节点不平衡量
		vector<float> delta_ef(a, 0);
		vector<float> delta_PQU(a, 0);

		// 填充delta_PQ
		bool flag = false;
		for (int i = 0; i < n; i++)
		{
			// 平衡节点不计算
			if (nodes[i].type == 3)
			{
				flag = true;
				continue;
			}
			// PQ节点带入deltaP和deltaQ
			if (nodes[i].type == 1)
			{
				delta_PQU[flag ? 2 * (i - 1) : 2 * i] = deltaP[i];
				delta_PQU[flag ? 2 * (i - 1) + 1 : 2 * i + 1] = deltaQ[i];
			} 
			// PV节点带入deltaP和deltaU
			else if(nodes[i].type == 2)
			{
				delta_PQU[flag ? 2 * (i - 1) : 2 * i] = deltaP[i];
				delta_PQU[flag ? 2 * (i - 1) + 1 : 2 * i + 1] = deltaU[i];
			}
		}

		for (int i = 0; i < a; i++)
		{
			for (int j = 0; j < a; j++)
			{
				delta_ef[i] += inverseJAC[i][j] * delta_PQU[j];
			}
		}

		// 输出Δe、Δf
		outfile << "\n[ ";
		for (int i = 0; i < n; i++)
		{
			if (nodes[i].type == 3) continue;
			outfile << "Δe" << i + 1 << " " << "Δf" << i + 1 << " ";
		}
		outfile << "]:\n";
		for (int i = 0; i < a; i++)
		{	
			outfile << right << setw(12) << delta_ef[i] << " ";
		}
		outfile << endl;

		// 修正节点电压值
		flag = false;
		for (int i = 0; i < n; i++)
		{
			// 平衡节点不计算
			if (nodes[i].type == 3)
			{
				flag = true;
				continue;
			}

			nodes[i].f += delta_ef[flag ? 2 * (i - 1) : 2 * i];
			nodes[i].e += delta_ef[flag ? 2 * (i - 1) + 1 : 2 * i + 1];
		}

		// 达到目标精度后停止迭代
		if (abs(delta_ef[0]) < epsilon && abs(delta_ef[1]) < epsilon)
		{
			converged = true;
			break;
		}

		iter++;
	}

	outfile << "\n***潮流计算结果***" << endl;

	if (!converged)
	{
		outfile << "\n" << maxIter << "次迭代后未收敛" << endl;
		return 0;
	}
	else
	{
		outfile << "\n迭代次数：" << iter + 1 << endl;
	}

	// 输出收敛后的节点电压
	outfile << "\n收敛后的节点电压：" << endl;
	for (int i = 0; i < n; i++)
	{
		outfile << "节点" << i + 1 << "：";
		outfile << "U=" << complex(nodes[i].e, nodes[i].f) << endl;
	}

	// 输出各节点注入功率
	outfile << "\n各节点注入功率：" << endl;
	for (int i = 0; i < n; i++)
	{
		outfile << "节点" << i + 1 << "：";
		outfile << "S=" << complex(nodes[i].P, nodes[i].Q) << endl;
	}

	// ----计算平衡节点功率和线路功率----

	// 计算平衡节点功率
	complex S_b = complex(0, 0);

	outfile << endl;
	for (int i = 0; i < n; i++)
	{
		if (nodes[i].type != 3) continue;

		for (int j = 0; j < n; j++)
		{
			S_b += complex(nodes[j].e, nodes[j].f).conj() * Y[i][j].conj();
		}

		S_b = S_b * complex(nodes[i].e, nodes[i].f);
	}

	// 计算线路功率
	vector<vector<complex>> S_l(n, vector<complex>(n, complex(0, 0)));

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) continue;
			for (const Branch branch : branches)
			{
				if (branch.node1 == i + 1 && branch.node2 == j + 1)
				{
					S_l[i][j] = complex(nodes[i].e, nodes[i].f) * (complex(0, 0) + (complex(nodes[i].e, nodes[i].f).conj() - complex(nodes[j].e, nodes[j].f).conj()) * branch.Y.conj());
				}
				else if (branch.node1 == j + 1 && branch.node2 == i + 1)
				{
					S_l[i][j] = complex(nodes[j].e, nodes[j].f) * (complex(0, 0) + (complex(nodes[j].e, nodes[j].f).conj() - complex(nodes[i].e, nodes[i].f).conj()) * branch.Y.conj());
				}
			}
		}
	}

	// 输出线路功率
	outfile << "各支路两端流入功率功率：" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (S_l[i][j] == complex(0, 0)) continue;
			outfile << left << "S_" << i + 1 << j + 1 << "=" << S_l[i][j] << endl;
		}
	}
	outfile << endl;

	// 网络总损耗
	complex delta_S = S_b;

	for (int i = 0; i < n; i++)
	{
		delta_S += complex(nodes[i].P, nodes[i].Q);
	}

	outfile << "网络总损耗：\n" <<  "ΔS=" << delta_S << endl;

	outfile.close();

	cout << "潮流计算完成，打开输出文件查看。" << endl;

	return 0;
}
