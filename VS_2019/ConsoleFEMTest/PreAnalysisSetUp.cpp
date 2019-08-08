#include "pch.h"
#include "PreAnalysisSetUp.h"

PreAnalysisSetUp::PreAnalysisSetUp()
{
}

PreAnalysisSetUp::~PreAnalysisSetUp()
{
}

PreAnalysisSetUp::PreAnalysisSetUp(StructureManager* structManager)
{
	structDetails = structManager;
}

std::vector<int> PreAnalysisSetUp::IdentifyIncrementalLoads()
{
	std::vector<int> indexes;
	int DOF = 6;

	std::map<int, Load*>::const_iterator it;

	for (it = structDetails->Loads().begin(); it != structDetails->Loads().end(); it++) { //for each load
		if (it->second->GetStatus() == "increment") {
			for (int j = 0; j < it->second->GetLoadVector().size(); j++) {
				indexes.emplace_back((it->second->GetNode() - 1) * DOF + (it->second->GetLoadVector()[j][0] - 1));
			}
		}
	}

	return indexes;
}

void PreAnalysisSetUp::CalculateShellsGlobalDOFVector()
{
	std::map<int, ShellElement*>::const_iterator it = structDetails->ShellElements;
	while (it != structDetails->ShellElements().end()) {
		it->second->CalculateGlobalDOFVector(structDetails->Supports);
		it++;
	}
}

void PreAnalysisSetUp::CalculateShellsGlobalMassDOFVector()
{
	std::map<int, ShellElement*>::const_iterator it = structDetails->ShellElements;
	while (it != structDetails->ShellElements().end()) {
		it->second->CalculateGlobalMassDOFVector(structDetails->Masses, structDetails->Supports);
		it++;
	}
}
