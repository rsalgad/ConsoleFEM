#pragma once
#include "Matrix.h"
#include "Support.h"
#include "StructureManager.h"
#include <vector>

class Mass
{
public:
	Mass();
	Mass(int ID, int nodeID);
	std::string ToString();
	int GetID();
	int GetNode();
	void SetMx(double val);
	void SetMy(double val);
	void SetMz(double val);
	bool operator ==(Mass const &m2);
	std::vector<std::vector<double>> GetMassVector();
	void SetMassVector(std::vector<std::vector<double>> vec);
	~Mass();

	static void SortByNodeID(std::vector<Mass> &mass);
	static void AddExplicitMassesOnExistingMatrix(Matrix* shellMass, const StructureManager* structManager, const int* DOF);
	static bool HasDOFAppliedMass(std::vector<Mass*> mass, int DOF);
private:
	int _ID, _nodeID;
	std::vector<std::vector<double>> _mass; //{{dir,load} = {1,20}, {2,10}} where 1 is Fx, 2 = Fy, 3 = Fz, 4 = Mx, 5 = My; 
};

