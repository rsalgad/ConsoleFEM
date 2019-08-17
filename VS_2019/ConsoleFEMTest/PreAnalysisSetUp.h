#pragma once
//#include "pch.h"
#include <vector>
#include <map>
#include <thread>
#include "StructureManager.h"
#include "AnalysisMethod.h"
#include "ShellElement.h"
#include "Matrix.h"

class PreAnalysisSetUp
{
public:
	PreAnalysisSetUp();
	~PreAnalysisSetUp();
	PreAnalysisSetUp(StructureManager *structManager, AnalysisMethod *analysisMethod);
	
	std::vector<int> IdentifyIncrementalLoads();

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
	AnalysisMethod* Analysis() const;
	
private:
	void CalculateShellsGlobalDOFVector();
	void CalculateShellsGlobalMassDOFVector();
	void CalculateSpringsGlobalDOFVector();
	void SetUpThreadsForShells();
	void CalculateReducedStiffMatrixSize();
	void CalculateStiffMatrixSize();
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
	AnalysisMethod* _analysisMethod;
	std::vector<int> _dispLoadDOFs;
};

