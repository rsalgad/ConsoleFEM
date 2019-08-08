#pragma once
#include "StructureManager.h"

class PreAnalysisSetUp
{
public:
	PreAnalysisSetUp();
	~PreAnalysisSetUp();
	PreAnalysisSetUp(StructureManager* structManager);
	
	std::vector<int> IdentifyIncrementalLoads();
	void CalculateShellsGlobalDOFVector();
	void CalculateShellsGlobalMassDOFVector();



private:
	StructureManager* structDetails = nullptr;
};

