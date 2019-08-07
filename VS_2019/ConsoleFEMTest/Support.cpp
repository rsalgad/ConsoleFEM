#include "pch.h"
#include "Support.h"
#include "Matrix.h"
#include <vector>


Support::Support()
{
}

Support::Support(int ID, int nodeID)
{
	_ID = ID;
	_nodeID = nodeID;
}

int Support::GetID()
{
	return _ID;
}

int Support::GetNode()
{
	return _nodeID;
}

std::string Support::ToString()
{
	std::string sup = "";
	sup += "(";
	sup += std::to_string(_ID);
	sup += ")";
	sup += "(";
	sup += "Node: ";
	sup += std::to_string(_nodeID);
	sup += ", ";
	for (int i = 0; i < _support[0].size(); i++) {

		if (i != 0)
		{
			sup += ", ";
		}

		sup += "Component ";
		sup += (i + 1);
		sup += ": ";
		int dir = _support[i][0];

		switch (dir)
		{
		case 1:
			sup += "Translation X = ";
			break;
		case 2:
			sup += "Translation Y = ";
			break;
		case 3:
			sup += "Translation Z = ";
			break;
		case 4:
			sup += "Rotation X = ";
			break;
		case 5:
			sup += "Rotation Y = ";
			break;
		case 6:
			sup += "Rotation Z = ";
			break;
		default:
			break;
		}
		sup += std::to_string(_support[i][1]);
	}

	sup += ")";

	return sup;
}

void Support::Set_tX(double val) {
	std::vector<double> vec = { 1, val };
	_support.push_back(vec);
}

void Support::Set_tY(double val) {
	std::vector<double> vec = { 2, val };
	_support.push_back(vec);
}

void Support::Set_tZ(double val) {
	std::vector<double> vec = { 3, val };
	_support.push_back(vec);
}

void Support::Set_rX(double val) {
	std::vector<double> vec = { 4, val };
	_support.push_back(vec);
}

void Support::Set_rY(double val) {
	std::vector<double> vec = { 5, val };
	_support.push_back(vec);
}

void Support::Set_rZ(double val) {
	std::vector<double> vec = { 6, val };
	_support.push_back(vec);
}

bool Support::operator ==(Support const &m2) {
	bool stat = false;
	if (_nodeID == m2._nodeID) {
		stat = true;
	}
	return stat;
}

int Support::NumberOfDOFBeforeNode(int nodeID, std::vector<Support> &sup) {
	int counter = 0;
	for (int i = 0; i < sup.size(); i++) {
		if (sup[i].GetNode() < nodeID) {
			counter += sup[i].GetSupportVector().size();
		}
	}
	return counter;
}

int Support::NumberOfDOFBeforeDOF(int DOF, std::vector<Support*> sup) {
	int counter = 0;
	int size = 6;
	for (int i = 0; i < sup.size(); i++) {
		for (int j = 0; j < sup[i]->GetSupportVector().size(); j++) {
			if (sup[i]->GetSupportVector()[j][1] == 0) { //I only want to do this if the support is a boundary condition
				int supDOF = (sup[i]->GetNode() - 1) * size + (sup[i]->GetSupportVector()[j][0] - 1);
				if (supDOF < DOF) {
					counter++;
				}
			}
		}
	}
	return counter;
}

bool Support::IsDOFConstrained(int DOF, std::vector<Support*> sup) {
	int size = 6;
	for (int i = 0; i < sup.size(); i++) {
		for (int j = 0; j < sup[i]->GetSupportVector().size(); j++) {
			if (sup[i]->GetSupportVector()[j][1] == 0) { //only constrained if displacement is 0
				int supDOF = (sup[i]->GetNode() - 1) * size + (sup[i]->GetSupportVector()[j][0] - 1);
				if (supDOF == DOF) {
					return true;
				}
			}
		}
	}
	return false;
}

int Support::TotalDOFsRestrained(std::vector<Support> &sup) {
	int counter = 0;
	for (int i = 0; i < sup.size(); i++) {
		for (int j = 0; j < sup[i].GetSupportVector().size(); j++) {
			if (sup[i].GetSupportVector()[j][1] == 0) { //if is an actual support
				counter++;
			}
		}
	}
	return counter;
}

void Support::SortByNodeID(std::vector<Support> &sup) {
	std::vector<Support> sup1;
	int index = 0, count = 0;
		
	for (int j = 0; j < sup.size(); j++) { //for each support
		int min = INT32_MAX;
		for (int i = 0; i < sup.size(); i++) { //for each support
			int nodeID = sup[i].GetNode();
			if (count == 0) {
				if (nodeID < min) {
					min = nodeID;
					index = i;
				}
			}
			else {
				if (!(std::find(sup1.begin(), sup1.end(), sup[i]) != sup1.end())) {
					//if is not inside the sup1 vector
					if (nodeID < min) {
						min = nodeID;
						index = i;
					}
				}
			}
		}
		sup1.push_back(sup[index]);
		count++;
	}

	/*
	std::vector<Support> sup1;
	int maxNodeID = 1;
	int count = 0;
	int iter = 0;
	while (true) {
		if (sup[iter].GetNode() == maxNodeID) {
			sup1.push_back(sup[iter]);
			maxNodeID++;
			count++;
		}
		if (count == sup.size()) {
			break;
		}
		iter++;
		if (iter >= sup.size()) {
			iter = 0;
		}
	}
	*/
	sup = sup1;
}

bool Support::IsNodeConstrained(std::vector<Support*> sup, int nodeID) {
	for (int i = 0; i < sup.size(); i++) { //for each support
		if (sup[i]->GetNode() == nodeID) {
			return true;
		}
	}
	return false;
}

std::vector<std::vector<double>> Support::GetSupportVector()
{
	return _support;
}

void Support::SetSupportVector(std::vector<std::vector<double>> vec)
{
	_support = vec;
}

std::vector<int> Support::GetDisplacementLoadIndexes(std::vector<Support*> vecSup) {
	std::vector<int> vec;
	int DOF = 6;
	for (int i = 0; i < vecSup.size(); i++) { //for each support
		for (int j = 0; j < vecSup[i]->GetSupportVector().size(); j++) {
			if (vecSup[i]->GetSupportVector()[j][1] != 0) { //if zero, then it is a support, not a support load!
				int nodeID = vecSup[i]->GetNode();
				int index = (nodeID - 1)*DOF + (vecSup[i]->GetSupportVector()[j][0] - 1);
				int count = Support::NumberOfDOFBeforeDOF(index, vecSup);
				vec.push_back(index - count);
			}
		}
	}
	return vec;
}



Support::~Support()
{
}
