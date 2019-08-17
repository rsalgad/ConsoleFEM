#include "pch.h"


Node::Node()
{
}

Node::Node(int ID, double x, double y, double z)
{
	_ID = ID;
	_x = x;
	_y = y;
	_z = z;
	_rx = 0;
	_ry = 0;
	_rz = 0;
}

Node::Node(int ID, double x, double y, double z, double rx, double ry, double rz)
{
	_ID = ID;
	_x = x;
	_y = y;
	_z = z;
	_rx = rx;
	_ry = ry;
	_rz = rz;
}

Node Node::FindNodeByID(int & ID, std::vector<Node>& vecNode) {
	//this should be optmized
	for (int i = 0; i < vecNode.size(); i++) {
		if (*vecNode[i].GetID() == ID) {
			return vecNode[i];
		}
	}
}

Node Node::FindNodeByCoordinates(double x, double y, double z, std::vector<Node>& vecNode) {
	//this should be optmized
	for (int i = 0; i < vecNode.size(); i++) {
		if (*vecNode[i].GetX() == x && *vecNode[i].GetY() == y && *vecNode[i].GetZ() == z) {
			return vecNode[i];
		}
	}
}

std::string Node::ToString() const {
	std::string row = "";
	row += "(";
	row += std::to_string(_ID);
	row += ")";
	row += "(";
	row += std::to_string(_x);
	row += ",";
	row += std::to_string(_y);
	row += ",";
	row += std::to_string(_z);
	row += ")";
	return row;
}

const int* Node::GetID()  const
{
	return &_ID;
}

const double* Node::GetX() const{
	return &_x;
}

const double* Node::GetY() const {
	return &_y;
}

const double* Node::GetZ() const {
	return &_z;
}

const double* Node::GetRx() const {
	return &_rx;
}

const double* Node::GetRy() const {
	return &_ry;
}

const double* Node::GetRz() const {
	return &_rz;
}

Node::~Node()
{
}
