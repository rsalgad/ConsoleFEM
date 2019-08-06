#pragma once
#include <vector>
#include <string>

class Node
{
public:
	Node();
	Node(int ID, double x, double y, double z);
	Node(int ID, double x, double y, double z, double rx, double ry, double rz);
	~Node();
	int GetID();
	double GetX();
	double GetY();
	double GetZ();
	double GetRx();
	double GetRy();
	double GetRz();
	std::string ToString();

	static Node FindNodeByID(int &ID, std::vector<Node> &vecNode);
	static Node FindNodeByCoordinates(double x, double y, double z, std::vector<Node>& vecNode);

private:
	int _ID;
	double _x, _y, _z, _rx, _ry, _rz;
};

