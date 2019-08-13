#include "pch.h"
#include "Mass.h"
#include "MatrixOperation.h"
#include "StructureManager.h"

Mass::Mass()
{
}

Mass::Mass(int ID, int nodeID)
{
	_ID = ID;
	_nodeID = nodeID;
}

std::string Mass::ToString()
{
	std::string mass = "";
	mass += "(";
	mass += std::to_string(_ID);
	mass += ")";
	mass += "(";
	mass += "Node: ";
	mass += std::to_string(_nodeID);
	mass += ", ";

	for (int i = 0; i < _mass[0].size(); i++) {

		if (i != 0)
		{
			mass += ", ";
		}

		mass += "Component ";
		mass += (i + 1);
		mass += ": ";
		int dir = _mass[i][0];

		switch (dir)
		{
		case 1:
			mass += "Translation X = ";
			break;
		case 2:
			mass += "Translation Y = ";
			break;
		case 3:
			mass += "Translation Z = ";
			break;
		case 4:
			mass += "Rotation X = ";
			break;
		case 5:
			mass += "Rotation Y = ";
			break;
		case 6:
			mass += "Rotation Z = ";
			break;
		default:
			break;
		}

		mass += std::to_string(_mass[i][1]);

	}

	mass += ")";

	return mass;
}

int Mass::GetID()
{
	return _ID;
}

int Mass::GetNode()
{
	return _nodeID;
}

void Mass::SetMx(double val) {
	std::vector<double> vec = { 1, val };
	_mass.push_back(vec);
}

void Mass::SetMy(double val) {
	std::vector<double> vec = { 2, val };
	_mass.push_back(vec);
}

void Mass::SetMz(double val) {
	std::vector<double> vec = { 3, val };
	_mass.push_back(vec);
}

std::vector<std::vector<double>> Mass::GetMassVector() {
	return _mass;
}

void Mass::SetMassVector(std::vector<std::vector<double>> vec) {
	_mass = vec;
}

bool Mass::operator ==(Mass const &m2) {
	bool stat = false;
	if (_nodeID == m2._nodeID) {
		stat = true;
	}
	return stat;
}

void Mass::SortByNodeID(std::vector<Mass> &mass) {
	std::vector<Mass> mass1;
	int index = 0, count = 0;

	for (int j = 0; j < mass.size(); j++) { //for each support
		int min = INT32_MAX;
		for (int i = 0; i < mass.size(); i++) { //for each support
			int nodeID = mass[i].GetNode();
			if (count == 0) {
				if (nodeID < min) {
					min = nodeID;
					index = i;
				}
			}
			else {
				if (!(std::find(mass1.begin(), mass1.end(), mass[i]) != mass1.end())) {
					//if is not inside the sup1 vector
					if (nodeID < min) {
						min = nodeID;
						index = i;
					}
				}
			}
		}
		mass1.push_back(mass[index]);
		count++;
	}
	mass = mass1;
}

bool Mass::HasDOFAppliedMass(std::vector<Mass*> mass, int DOF) {
	int nDOF = 6;
	for (int i = 0; i < mass.size(); i++) { //for each mass
		int nodeID = mass[i]->GetNode();
		for (int j = 0; j < mass[i]->GetMassVector().size(); j++) {
			int calcDof = (nodeID - 1)*nDOF + (mass[i]->GetMassVector()[j][0] - 1);
			if (calcDof == DOF) {
				return true;
			}
		}
	}
	return false;
}

void Mass::AddExplicitMassesOnExistingMatrix(Matrix* shellMass, const StructureManager* structManager, const int* DOF)
{
	std::map<int, Mass*>::const_iterator it = structManager->Masses()->begin();

	while (it != structManager->Masses()->end()) {
		Mass* m = it->second;
		int nodeID = m->GetNode();
		for (int j = 0; j < m->GetMassVector().size(); j++) {
			int dof = (nodeID - 1) * (*DOF) + (m->GetMassVector()[j][0] - 1);
			int val = Support::NumberOfDOFBeforeDOF(&dof, structManager->Supports());
			shellMass->GetMatrixDouble()[dof - val][dof - val] += m->GetMassVector()[j][1];
		}
		it++;
	}
}

Mass::~Mass()
{
}
