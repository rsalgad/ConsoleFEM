#pragma once
#include "StructureManager.h"
#include <map>

class PreAnalysisSetUp
{
public:
	PreAnalysisSetUp();
	~PreAnalysisSetUp();
	PreAnalysisSetUp(const StructureManager* structManager);
	
	std::vector<int> IdentifyIncrementalLoads();
	void CalculateShellsGlobalDOFVector();
	void CalculateShellsGlobalMassDOFVector();
	const int* ReducedStiffMatrixSize() const;
	const int* StiffMatrixSize() const;
	const int* DOF() const;
	const std::map<int, std::vector<ShellElement*>>* GetShellThreads() const;
	const int* AvailableThreads() const;
	
private:
	void SetUpThreadsForShells(const std::map<int, ShellElement*>* listOfShells);
	void CalculateReducedStiffMatrixSize();
	const StructureManager* _structDetails = nullptr;
	std::map<int, std::vector<ShellElement*>> _shellThreads;
	const int _nThreads = std::thread::hardware_concurrency();
	int _redStiffMatrixSize, _stiffMatrixSize;
	const int _DOF = 6;
};

