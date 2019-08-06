#pragma once
#include "Matrix.h"
#include "Support.h"
#include <vector>

class Mass
{
public:
	Mass();
	Mass(int ID, int nodeID);
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
	static Matrix AddExplicitMassesOnExistingMatrix(Matrix &shellMass, std::vector<Mass> &listOfMasses, std::vector<Support>& listOfSups);
	static bool HasDOFAppliedMass(std::vector<Mass> &mass, int DOF);
private:
	int _ID, _nodeID;
	std::vector<std::vector<double>> _mass; //{{dir,load} = {1,20}, {2,10}} where 1 is Fx, 2 = Fy, 3 = Fz, 4 = Mx, 5 = My; 
};

