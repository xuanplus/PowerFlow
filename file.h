#pragma once
#include <string>
#include <fstream>
#include <sstream>


// ���ļ����ж�ȡ����
template<typename T>
T read_data(std::ifstream& file)
{
	T data;
	std::string line;
	getline(file, line);
	std::istringstream iss(line);
	iss >> data;
	return data;
}


// ���ļ���ȡ�ڵ�
std::pair<int, int> read_node(std::ifstream& file)
{
	std::string line;
	std::getline(file, line);
	const size_t pos = line.find('-');
	int node1 = stoi(line.substr(0, pos));
	int node2 = stoi(line.substr(pos + 1));
	return std::make_pair(node1, node2);
}
