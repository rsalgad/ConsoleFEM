#pragma once
//#include "pch.h"
#include <string>
#include <map>
#include <vector>

class Support
{
public:
	Support();
	Support(int ID, int nodeID);
	int GetID();
	int GetNode();
	std::string ToString();
	void Set_tX(double val);
	void Set_tY(double val);
	void Set_tZ(double val);
	void Set_rX(double val);
	void Set_rY(double val);
	void Set_rZ(double val);
	bool operator ==(Support const &m2);
	std::vector<std::vector<double>> GetSupportVector();
	void SetSupportVector(std::vector<std::vector<double>> vec);
	~Support();

	static bool IsNodeConstrained(const std::map<int, Support*>* sup, const int* nodeID);
	static void SortByNodeID(std::vector<Support> &sup);
	static int NumberOfDOFBeforeNode(int nodeID, std::vector<Support> &sup);
	static int NumberOfDOFBeforeDOF(const int* DOF, const std::map<int, Support*>* sup, const int* nDOF);
	static bool IsDOFConstrained(const int* DOF, const std::map<int, Support*>* sup, const int* nDOF);
	static int TotalDOFsRestrained(const std::map<int, Support*>* sup);
	static std::vector<int> GetDisplacementLoadIndexes(const int* DOF, const std::map<int, Support*>* vecSup);
private:
	int _ID, _nodeID;
	std::vector<std::vector<double>> _support; //{{dir,sup} = {1,0}, {2,0}} where 1 is X, 2 is Y, 3 is Z, 4 is rotation X, 5 is rotation Y;
};

