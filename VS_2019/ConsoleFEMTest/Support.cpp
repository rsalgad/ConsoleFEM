#include "pch.h"


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

int Support::NumberOfDOFBeforeDOF(const int* DOF, const std::map<int, Support*>* sup, const int* nDOF) {
	int counter = 0;
	std::map<int, Support*>::const_iterator it = sup->begin();

	while (it != sup->end()) {
		for (int j = 0; j < it->second->GetSupportVector().size(); j++) {
			if (it->second->GetSupportVector()[j][1] == 0) { //I only want to do this if the support is a boundary condition
				int supDOF = (it->second->GetNode() - 1) * (*nDOF) + (it->second->GetSupportVector()[j][0] - 1);
				if (supDOF < (*DOF)) {
					counter++;
				}
			}
		}
		it++;
	}
	return counter;
}

bool Support::IsDOFConstrained(const int* DOF, const std::map<int, Support*>* sup, const int* nDOF) {
	std::map<int, Support*>::const_iterator it = sup->begin();

	while (it != sup->end()) {

		for (int j = 0; j < it->second->GetSupportVector().size(); j++) {
			if (it->second->GetSupportVector()[j][1] == 0) { //only constrained if displacement is 0
				int supDOF = (it->second->GetNode() - 1) * (*nDOF) + (it->second->GetSupportVector()[j][0] - 1);
				if (supDOF == *DOF) {
					return true;
				}
			}
		}
		it++;
	}
	return false;
}

int Support::TotalDOFsRestrained(const std::map<int, Support*>* sup) {
	int counter = 0;
	std::map<int, Support*>::const_iterator it = sup->begin();
	while (it != sup->end()) {
		for (int j = 0; j < it->second->GetSupportVector().size(); j++) {
			if (it->second->GetSupportVector()[j][1] == 0) { //if is an actual support
				counter++;
			}
		}
		it++;
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

bool Support::IsNodeConstrained(const std::map<int, Support*>* sup,const int* nodeID) {
	
	std::map<int, Support*>::const_iterator it = sup->begin();

	while (it != sup->end()) {//for each support
		if (it->second->GetNode() == *nodeID) {
			return true;
		}
		it++;
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

std::vector<int> Support::GetDisplacementLoadIndexes(const int* DOF, const std::map<int, Support*>* vecSup) {
	std::vector<int> vec;

	std::map<int, Support*>::const_iterator it = vecSup->begin();

	while (it != vecSup->end()) {//for each support

		for (int j = 0; j < it->second->GetSupportVector().size(); j++) {
			if (it->second->GetSupportVector()[j][1] != 0) { //if zero, then it is a support, not a support load!
				int nodeID = it->second->GetNode();
				int index = (nodeID - 1) * (*DOF) + (it->second->GetSupportVector()[j][0] - 1);
				int count = Support::NumberOfDOFBeforeDOF(&index, vecSup, DOF);
				vec.push_back(index - count);
			}
		}
		it++;
	}
	return vec;
}



Support::~Support()
{
}
