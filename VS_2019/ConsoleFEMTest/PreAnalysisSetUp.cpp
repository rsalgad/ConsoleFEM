#include "pch.h"
#include "AnalysisMethod.h"
#include "PreAnalysisSetUp.h"

PreAnalysisSetUp::PreAnalysisSetUp()
{
}

PreAnalysisSetUp::~PreAnalysisSetUp()
{
}

PreAnalysisSetUp::PreAnalysisSetUp(const StructureManager* structManager, const AnalysisMethod* analysisMethod, int nLoadSteps, int nIterations)
{
	_structDetails = structManager;
	_analysisMethod = analysisMethod;
	_nLoadSteps = nLoadSteps;
	_nIterations = nIterations;
}

std::vector<int> PreAnalysisSetUp::IdentifyIncrementalLoads()
{
	std::vector<int> indexes;
	int DOF = 6;

	std::map<int, Load*>::const_iterator it;

	for (it = _structDetails->Loads()->begin(); it != _structDetails->Loads()->end(); it++) { //for each load
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
	std::map<int, ShellElement*>::const_iterator it = _structDetails->ShellElements;
	while (it != _structDetails->ShellElements()->end()) {
		it->second->CalculateGlobalDOFVector(_structDetails->Supports);
		it++;
	}
}

void PreAnalysisSetUp::CalculateShellsGlobalMassDOFVector()
{
	std::map<int, ShellElement*>::const_iterator it = _structDetails->ShellElements;
	while (it != _structDetails->ShellElements()->end()) {
		it->second->CalculateGlobalMassDOFVector(_structDetails->Masses, _structDetails->Supports);
		it++;
	}
}

const int* PreAnalysisSetUp::ReducedStiffMatrixSize() const
{
	return &_redStiffMatrixSize;
}

const int* PreAnalysisSetUp::StiffMatrixSize() const
{
	return &_stiffMatrixSize;
}

const int* PreAnalysisSetUp::DOF() const
{
	return &_DOF;
}

const std::map<int, std::vector<ShellElement*>>* PreAnalysisSetUp::GetShellThreads() const
{
	return &_shellThreads;
}

const int* PreAnalysisSetUp::AvailableThreads() const
{
	return &_nThreads;
}

const std::vector<double>* PreAnalysisSetUp::LoadFactors() const
{
	return &_loadFactors;
}

const Matrix* PreAnalysisSetUp::ConstForces() const
{
	return &_constForces;
}

const Matrix* PreAnalysisSetUp::IncForces() const
{
	return &_incrForces;
}

const int* PreAnalysisSetUp::LoadSteps() const
{
	switch (_analysisMethod->Type())
	{
	case AnalysisTypes::Elastic:
		return &_nLoadSteps;
		break;
	case AnalysisTypes::Monotonic:
		return &_nLoadSteps;
		break;
	case AnalysisTypes::Cyclic:
		return &_nLoadSteps;
		break;
	case AnalysisTypes::Reverse_Cyclic:
		break;
	case AnalysisTypes::Seismic:
		break;
	case AnalysisTypes::Impulse:
		break;
	default:
		return nullptr;
		break;
	}
}

//<summary>Sets up the list of elements that each thread will be responsible for</summary>
void PreAnalysisSetUp::SetUpThreadsForShells(const std::map<int, ShellElement*>* listOfShells) {
	int nThreads = std::thread::hardware_concurrency();
	double amount = listOfShells->size() / nThreads;
	int size1, size2;
	
	if (fmod(listOfShells->size(), nThreads) != 0) {
		size1 = floor(amount);
		size2 = fmod(amount, nThreads) * nThreads;
	}
	else {
		size1 = amount;
		size2 = 0;
	}

	for (int i = 0; i < nThreads; i++) {
		if (i < (nThreads - 1))
		{
			std::vector<ShellElement*> list;
			list.assign(listOfShells->begin()->second + i * size1, listOfShells->begin()->second + (i + 1) * size1);
			_shellThreads.insert(std::pair<int, std::vector<ShellElement*>>(i + 1, list));
		}
		else {
			std::vector<ShellElement*> list;
			list.assign(listOfShells->begin()->second + i * size1, listOfShells->end()->second);
			_shellThreads.insert(std::pair<int, std::vector<ShellElement*>>(i+1, list));
		}
	}
}

void PreAnalysisSetUp::CalculateReducedStiffMatrixSize()
{
	_redStiffMatrixSize = _structDetails->Nodes()->size() * _DOF - Support::TotalDOFsRestrained(_structDetails->Supports());
}

void PreAnalysisSetUp::CalculateForceMatrices()
{
	_constForces = Load::AssembleLoadMatrixWithFlag(_structDetails, this, "constant");
	_incrForces = Load::AssembleLoadMatrixWithFlag(_structDetails, this, "increment");
}

void PreAnalysisSetUp::CalculateLoadFactors()
{
	switch (_analysisMethod)
	{
	case AnalysisMethod::Elastic:
	{
		_loadFactors.reserve(_nLoadSteps);
		for (int i = 0; i < _nLoadSteps; i++) {
			_loadFactors.push_back((i + 1) * (1.0 / _nLoadSteps));
		}
		break;
	}
	default:
		break;
	}
}
