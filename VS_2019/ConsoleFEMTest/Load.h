#pragma once
#include <vector>
#include "Matrix.h"
#include "Node.h"
#include "Support.h"
#include "StructureManager.h"
#include "PreAnalysisSetUp.h"

class Load
{
public:
	Load();
	Load(int ID, int nodeID);
	Load(int ID, int nodeID, std::string status);
	std::string ToString();
	int GetID();
	int GetNode();
	std::string GetStatus();
	void SetFx(double val);
	void SetFy(double val);
	void SetFz(double val);
	void SetMx(double val);
	void SetMy(double val);
	void SetMz(double val);
	bool operator ==(Load const &l2);
	std::vector<std::vector<double>> GetLoadVector();
	void SetLoadVector(std::vector<std::vector<double>> vec);
	~Load();

	static void SortByNodeID(std::vector<Load> &sup);
	static Matrix GetTotalForceMatrix(const Matrix* force, const Matrix* react, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const double* biggest, const double* loadFraction);
	static Matrix AssembleLoadMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp );
	static Matrix AssembleLoadMatrixWithFlag(const StructureManager* structManager, const PreAnalysisSetUp* setUp, std::string flag);
	static Matrix AssembleDispLoadMatrix(std::vector<Node> &vecNode, std::vector<Support> &vecSup);
	static Matrix GetReducedLoadMatrix(const Matrix* loadMatrix, const std::map<int, Support*>* mapSup, const int* DOF);
	static std::vector<int> IdentifyIncrementalLoads(std::vector<Load> &vecLoad);
	static Matrix MultiplyIncrementalTerms(Matrix &redLoadMatrix, std::vector<int> &incIndex, std::vector<Support> &vecSups, double mult);
	static Matrix GetTotalForceNotOrganized(const Matrix* force, const Matrix* react, const StructureManager* structManager, const int* stiffSize);
	static double DefineLoadFractionAtLoadStep(std::string type, int& step, int& totalSteps, int stepsPerPeak, double peakInc, int cyclesPerPeak, double iniPeak);
	static double DefineTotalLoadSteps(std::string type, int &step, int cyclicRepeat, double totalTime, double deltaT);
	static double SampleDynamicForceFunction(double amplitude, double period, double phase, double t);
	static std::map<int, Load*> GetNewLoads(const std::map<int, Load*>* oriMap, const Matrix* totalforce, const int* DOF);
	static Matrix GetForceByNodeID(int const* ID, const Matrix* totalForceMatrix, const int* DOF);

private:
	int _ID, _nodeID;
	std::vector<std::vector<double>> _load; //{{dir,load} = {1,20}, {2,10}} where 1 is Fx, 2 = Fy, 3 = Fz, 4 = Mx, 5 = My; 
	std::string _status = "constant";
};


