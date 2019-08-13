#include "pch.h"
#include "AnalysisMethod.h"
#include "PreAnalysisSetUp.h"
#include "SeismicLoad.h"
#include "ImpulseLoad.h"
#include "DynamicAnalysis.h"

PreAnalysisSetUp::PreAnalysisSetUp()
{
}

PreAnalysisSetUp::~PreAnalysisSetUp()
{
}

PreAnalysisSetUp::PreAnalysisSetUp(const StructureManager* structManager, const AnalysisMethod* analysisMethod)
{
	_structDetails = structManager;
	_analysisMethod = analysisMethod;
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
	return _analysisMethod->LoadSteps();
}

const int* PreAnalysisSetUp::Iterations() const
{
	return &_nIterations;
}

const std::vector<int>* PreAnalysisSetUp::DispLoadDOFs() const
{
	return &_dispLoadDOFs;
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
	switch (_analysisMethod->Type())
	{
	case AnalysisTypes::Seismic:
	{
		_incrForces = Load::AssembleLoadMatrixWithFlag(_structDetails, this, "seismic");
		break;
	}
	case AnalysisTypes::Impulse:
	{
		_incrForces = Load::AssembleLoadMatrixWithFlag(_structDetails, this, "impulse");
		break;
	}
	default:
		_incrForces = Load::AssembleLoadMatrixWithFlag(_structDetails, this, "increment");
		break;
	}

	_constForces = Load::AssembleLoadMatrixWithFlag(_structDetails, this, "constant");

}

void PreAnalysisSetUp::CalculateLoadFactors()
{
	switch (_analysisMethod->Type())
	{
	case AnalysisTypes::Cyclic: 
	{
		_loadFactors.reserve(_nLoadSteps);
		for (int step = 0; step < _nLoadSteps; step++) {
			int stepsPerCycle = stepsPerPeak * 2;//return the total number of steps inside each full cycle
			int cycle = step / stepsPerCycle; //return the number of cycles already covered
			int aux = step + 1 - cycle * stepsPerCycle; //returns the loadstep inside the current cycle
			int aux2 = aux / stepsPerPeak; // returns the number of peaks of the current cycle it already covered

			double currentPeak = iniPeak + (cycle / cyclesPerPeak) * peakInc;
			double increment = currentPeak / stepsPerPeak;

			if (aux2 == 0) {
				_loadFactors.push_back((aux)* increment);
			}
			else {
				_loadFactors.push_back(currentPeak - (aux - stepsPerPeak) * increment);
			}
		}
	}
	case AnalysisTypes::Reverse_Cyclic:
	{
		_loadFactors.reserve(_nLoadSteps);
		for (int step = 0; step < _nLoadSteps; step++) {
			int stepsPerCycle = stepsPerPeak * 4;//return the total number of steps inside each full cycle
			int cycle = step / stepsPerCycle; //return the number of cycles already covered
			int aux = step + 1 - cycle * stepsPerCycle; //returns the loadstep inside the current cycle
			int aux2 = aux / stepsPerPeak; // returns the number of peaks of the current cycle it already covered

			double currentPeak = iniPeak + (cycle / cyclesPerPeak) * peakInc;
			double increment = currentPeak / stepsPerPeak;

			if (aux2 == 0) {
				_loadFactors.push_back(aux * increment);
			}
			else if (aux2 == 1) {
				_loadFactors.push_back(currentPeak - (aux - stepsPerPeak) * increment);
			}
			else if (aux2 == 2) {
				_loadFactors.push_back(-(aux - aux2 * stepsPerPeak) * increment);
			}
			else if (aux2 == 3) {
				_loadFactors.push_back(-currentPeak + (aux - aux2 * stepsPerPeak) * increment);
			}
			else {
				_loadFactors.push_back(0);
			}
		}
	}
	case AnalysisTypes::Impulse:
	{
		const double* deltaT = static_cast<DynamicAnalysis*>(_analysisMethod)->DeltaT();
		_loadFactors.reserve(_nLoadSteps);
		for (int step = 0; step < _nLoadSteps; step++) {
			double time = (step + 1) * (*deltaT);

			if (time <= impLoad.GetPoints()[0][0]) {
				double val = impLoad.GetPoints()[0][1] * time / impLoad.GetPoints()[0][0];
				_loadFactors.push_back(val);
			}
			else if (time < impLoad.GetPoints()[impLoad.GetPoints().size() - 1][0]) {
				for (int i = 0; i < impLoad.GetPoints().size() - 1; i++) {
					if (time > impLoad._points[i][0] && time <= impLoad._points[i + 1][0]) {
						double val = (impLoad.GetPoints()[i][1] * (impLoad.GetPoints()[i + 1][0] - time) + impLoad.GetPoints()[i + 1][1] * (time - impLoad.GetPoints()[i][0])) / (impLoad.GetPoints()[i + 1][0] - impLoad.GetPoints()[i][0]);
						_loadFactors.push_back(val);
					}
				}
			}
			else {
				double val = impLoad.GetPoints()[impLoad.GetPoints().size() - 1][1]; //whatever the last point is is remained constant
				_loadFactors.push_back(val);
			}
		}
		break;
	}
	default:
		_loadFactors.reserve(_nLoadSteps);
		for (int step = 0; step < _nLoadSteps; step++) {
			_loadFactors.push_back((step + 1) * (1.0 / _nLoadSteps));
		}
		break;
	}
}

void PreAnalysisSetUp::CalculateDispLoadDOFs()
{
	_dispLoadDOFs = Support::GetDisplacementLoadIndexes(&_DOF, _structDetails->Supports());
}
