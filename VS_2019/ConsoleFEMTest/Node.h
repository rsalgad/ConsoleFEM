#pragma once
//#include "pch.h"
#include <string>
#include <vector>


class Node
{
public:
	Node();
	Node(int ID, double x, double y, double z);
	Node(int ID, double x, double y, double z, double rx, double ry, double rz);
	~Node();
	const int* GetID() const;
	const double* GetX() const;
	const double* GetY() const;
	const double* GetZ() const;
	const double* GetRx() const;
	const double* GetRy() const;
	const double* GetRz() const;
	std::string ToString() const;

	static Node FindNodeByID(int &ID, std::vector<Node> &vecNode);
	static Node FindNodeByCoordinates(double x, double y, double z, std::vector<Node>& vecNode);

private:
	int _ID;
	double _x, _y, _z, _rx, _ry, _rz;
};

