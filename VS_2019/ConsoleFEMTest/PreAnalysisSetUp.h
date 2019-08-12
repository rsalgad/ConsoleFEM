#pragma once
#include "StructureManager.h"
#include "AnalysisMethod.h"
#include <map>

class PreAnalysisSetUp
{
public:
	PreAnalysisSetUp();
	~PreAnalysisSetUp();
	PreAnalysisSetUp(const StructureManager* structManager, const AnalysisMethod* analysisMethod, int nLoadSteps, int nIterations );
	
	std::vector<int> IdentifyIncrementalLoads();
	void CalculateShellsGlobalDOFVector();
	void CalculateShellsGlobalMassDOFVector();
	const int* ReducedStiffMatrixSize() const;
	const int* StiffMatrixSize() const;
	const int* DOF() const;
	const std::map<int, std::vector<ShellElement*>>* GetShellThreads() const;
	const int* AvailableThreads() const;
	const std::vector<double>* LoadFactors() const;
	const Matrix* ConstForces() const;
	const Matrix* IncForces() const;
	const int* LoadSteps() const;
	const int* Iterations() const;
	const std::vector<int>* DispLoadDOFs() const;
	
private:
	void SetUpThreadsForShells(const std::map<int, ShellElement*>* listOfShells);
	void CalculateReducedStiffMatrixSize();
	void CalculateForceMatrices();
	void CalculateLoadFactors();
	void CalculateDispLoadDOFs();
	const StructureManager* _structDetails = nullptr;
	std::map<int, std::vector<ShellElement*>> _shellThreads;
	const int _nThreads = std::thread::hardware_concurrency();
	int _redStiffMatrixSize, _stiffMatrixSize;
	const int _DOF = 6;
	Matrix _constForces;
	Matrix _incrForces;
	std::vector<double> _loadFactors;
	const AnalysisMethod* _analysisMethod;
	int _nLoadSteps;
	int _nIterations;
	bool _breakAnalysis = false; //indicates if the analysis should be stopped
	std::vector<int> _dispLoadDOFs;
};

