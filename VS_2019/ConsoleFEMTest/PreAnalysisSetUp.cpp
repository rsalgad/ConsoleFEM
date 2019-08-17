#include "pch.h"
#include "CyclicAnalysis.h"
#include "ReverseCyclic.h"

PreAnalysisSetUp::PreAnalysisSetUp()
{
}

PreAnalysisSetUp::~PreAnalysisSetUp()
{
}

PreAnalysisSetUp::PreAnalysisSetUp(StructureManager* structManager, AnalysisMethod* analysisMethod)
{
	_structDetails = structManager;
	_analysisMethod = analysisMethod;
	structManager->SortSupportsByNodeID();
	CalculateStiffMatrixSize();
	CalculateReducedStiffMatrixSize();
	SetUpThreadsForShells();
	CalculateShellsGlobalDOFVector();
	CalculateShellsGlobalMassDOFVector();
	CalculateSpringsGlobalDOFVector();
	CalculateForceMatrices();
	CalculateLoadFactors();
	CalculateDispLoadDOFs();
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
	std::map<int, ShellElement*>::const_iterator it = _structDetails->ShellElements()->begin();
	while (it != _structDetails->ShellElements()->end()) {
		it->second->CalculateGlobalDOFVector(_structDetails->Supports(), &_DOF);
		it++;
	}
}

void PreAnalysisSetUp::CalculateShellsGlobalMassDOFVector()
{
	std::map<int, ShellElement*>::const_iterator it = _structDetails->ShellElements()->begin();
	while (it != _structDetails->ShellElements()->end()) {
		it->second->CalculateGlobalMassDOFVector(_structDetails, &_DOF);
		it++;
	}
}

void PreAnalysisSetUp::CalculateSpringsGlobalDOFVector()
{
	std::map<int, Spring3D*>::const_iterator it = _structDetails->SpringElements()->begin();
	while (it != _structDetails->SpringElements()->end()) {
		it->second->CalculateGlobalDOFVector(_structDetails->Supports(), &_DOF);
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
	return _analysisMethod->Iterations();
}

const std::vector<int>* PreAnalysisSetUp::DispLoadDOFs() const
{
	return &_dispLoadDOFs;
}

AnalysisMethod* PreAnalysisSetUp::Analysis() const
{
	return _analysisMethod;
}

//<summary>Sets up the list of elements that each thread will be responsible for</summary>
void PreAnalysisSetUp::SetUpThreadsForShells() {
	double nThreads = std::thread::hardware_concurrency();
	double amount = _structDetails->ShellElements()->size() / nThreads;
	int size1, size2;
	
	if (_structDetails->ShellElements()->size() >= 2) {

		double div = _structDetails->ShellElements()->size() / nThreads;
		if (div < 1) { //don't have enough enough elements to fill all cores
			std::map<int, ShellElement*>::const_iterator it = _structDetails->ShellElements()->begin();
			while (it != _structDetails->ShellElements()->end()) { //for each shell element -> 1 thread
				std::vector<ShellElement*> list;
				list.push_back(it->second);
				_shellThreads.insert(std::pair<int, std::vector<ShellElement*>>(it->first, list));
				it++;
			}
		}
		else {
			int remainingEle = fmod(_structDetails->ShellElements()->size(), nThreads);
			int elePerThread = (_structDetails->ShellElements()->size() - remainingEle) / (nThreads);

			int counter = 1;
			int totCounter = 0;

			for (int i = 0; i < nThreads; i++) { //for each thread
				std::vector<ShellElement*> list;
				if (counter <= remainingEle) {
					for (int j = 0; j < elePerThread + 1; j++) {
						totCounter++;
						list.push_back(_structDetails->ShellElements()->find(totCounter)->second);
					}
				}
				else {
					for (int j = 0; j < elePerThread; j++) {
						totCounter++;
						list.push_back(_structDetails->ShellElements()->find(totCounter)->second);
					}
				}
				_shellThreads.insert(std::pair<int, std::vector<ShellElement*>>(i + 1, list));
				counter++;
			}
		}
	}
	else {
		std::vector<ShellElement*> list;
		list.push_back(_structDetails->ShellElements()->find(1)->second);
		_shellThreads.insert(std::pair<int, std::vector<ShellElement*>>(1, list));
	}
}

void PreAnalysisSetUp::CalculateReducedStiffMatrixSize()
{
	_redStiffMatrixSize = _structDetails->Nodes()->size() * _DOF - Support::TotalDOFsRestrained(_structDetails->Supports());
}

void PreAnalysisSetUp::CalculateStiffMatrixSize()
{
	_stiffMatrixSize = _structDetails->Nodes()->size() * _DOF;
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
		_loadFactors.reserve(*_analysisMethod->LoadSteps());
		CyclicAnalysis* analysis = static_cast<CyclicAnalysis*>(_analysisMethod);

		for (int step = 0; step < *_analysisMethod->LoadSteps(); step++) {
			int stepsPerCycle = *analysis->StepsPerPeak() * 2;//return the total number of steps inside each full cycle
			int cycle = step / stepsPerCycle; //return the number of cycles already covered
			int aux = step + 1 - cycle * stepsPerCycle; //returns the loadstep inside the current cycle
			int aux2 = aux / *analysis->StepsPerPeak(); // returns the number of peaks of the current cycle it already covered

			double currentPeak = *analysis->IniPeak() + (cycle / *analysis->CyclesPerPeak()) * *analysis->PeakInc();
			double increment = currentPeak / *analysis->StepsPerPeak();

			if (aux2 == 0) {
				_loadFactors.push_back((aux)* increment);
			}
			else {
				_loadFactors.push_back(currentPeak - (aux - *analysis->StepsPerPeak()) * increment);
			}
		}
	}
	case AnalysisTypes::Reverse_Cyclic:
	{
		_loadFactors.reserve(*_analysisMethod->LoadSteps());
		ReverseCyclic* analysis = static_cast<ReverseCyclic*>(_analysisMethod);

		for (int step = 0; step < *_analysisMethod->LoadSteps(); step++) {
			int stepsPerCycle = *analysis->StepsPerPeak() * 4;//return the total number of steps inside each full cycle
			int cycle = step / stepsPerCycle; //return the number of cycles already covered
			int aux = step + 1 - cycle * stepsPerCycle; //returns the loadstep inside the current cycle
			int aux2 = aux / *analysis->StepsPerPeak(); // returns the number of peaks of the current cycle it already covered

			double currentPeak = *analysis->IniPeak() + (cycle / *analysis->CyclesPerPeak()) * *analysis->PeakInc();
			double increment = currentPeak / *analysis->StepsPerPeak();

			if (aux2 == 0) {
				_loadFactors.push_back(aux * increment);
			}
			else if (aux2 == 1) {
				_loadFactors.push_back(currentPeak - (aux - *analysis->StepsPerPeak()) * increment);
			}
			else if (aux2 == 2) {
				_loadFactors.push_back(-(aux - aux2 * *analysis->StepsPerPeak()) * increment);
			}
			else if (aux2 == 3) {
				_loadFactors.push_back(-currentPeak + (aux - aux2 * *analysis->StepsPerPeak()) * increment);
			}
			else {
				_loadFactors.push_back(0);
			}
		}
	}
	case AnalysisTypes::Impulse:
	{

		DynamicAnalysis* analysis = static_cast<DynamicAnalysis*>(_analysisMethod);
		ImpulseLoad* impLoad = static_cast<ImpulseLoad*>(analysis->Load());

		_loadFactors.reserve(*_analysisMethod->LoadSteps());
		for (int step = 0; step < *_analysisMethod->LoadSteps(); step++) {
			double time = (step + 1) * (*analysis->DeltaT());

			if (time <= impLoad->GetPoints()[0][0]) {
				double val = impLoad->GetPoints()[0][1] * time / impLoad->GetPoints()[0][0];
				_loadFactors.push_back(val);
			}
			else if (time < impLoad->GetPoints()[impLoad->GetPoints().size() - 1][0]) {
				for (int i = 0; i < impLoad->GetPoints().size() - 1; i++) {
					if (time > impLoad->GetPoints()[i][0] && time <= impLoad->GetPoints()[i + 1][0]) {
						double val = (impLoad->GetPoints()[i][1] * (impLoad->GetPoints()[i + 1][0] - time) + impLoad->GetPoints()[i + 1][1] * (time - impLoad->GetPoints()[i][0])) / (impLoad->GetPoints()[i + 1][0] - impLoad->GetPoints()[i][0]);
						_loadFactors.push_back(val);
					}
				}
			}
			else {
				double val = impLoad->GetPoints()[impLoad->GetPoints().size() - 1][1]; //whatever the last point is is remained constant
				_loadFactors.push_back(val);
			}
		}
		break;
	}
	default:
		_loadFactors.reserve(*_analysisMethod->LoadSteps());
		for (int step = 0; step < *_analysisMethod->LoadSteps(); step++) {
			_loadFactors.push_back((step + 1) * (1.0 / *_analysisMethod->LoadSteps()));
		}
		break;
	}
}

void PreAnalysisSetUp::CalculateDispLoadDOFs()
{
	_dispLoadDOFs = Support::GetDisplacementLoadIndexes(&_DOF, _structDetails->Supports());
}


