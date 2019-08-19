#include "pch.h"

AnalysisSpringRecorder::AnalysisSpringRecorder()
{
}

AnalysisSpringRecorder::AnalysisSpringRecorder(std::vector<int> springsID)
{
	_springsID = springsID;
}

AnalysisSpringRecorder::AnalysisSpringRecorder(char c, const std::map<int, Spring3D*>* allSprings)
{
	//Gets the spring IDs that will be stored
	if (c == 'a') {
		std::map<int, Node*>::iterator it;
		_springsID.reserve(allSprings->size());
		for (int i = 0; i < allSprings->size(); i++) {
			_springsID.push_back(i + 1);
		}
	}
	else {
		std::cout << "Error initializing recorder" << std::endl;
	}

	InitializeAllVectors();
}

AnalysisSpringRecorder::~AnalysisSpringRecorder()
{
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetNewDisp()
{
	return &_newDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetOldDisp()
{
	return &_oldDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetNewMinDisp()
{
	return&_newMinDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetOldMinDisp()
{
	return &_oldMinDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetNewMaxDisp()
{
	return &_newMaxDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetOldMaxDisp()
{
	return &_oldMaxDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetNewPlasticDisp()
{
	return &_newPlasticDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetOldPlasticDisp()
{
	return &_oldPlasticDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetNewUnlDisp()
{
	return &_newUnlDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetOldUnlDisp()
{
	return &_oldUnlDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetNewRelDisp()
{
	return &_newRelDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetOldRelDisp()
{
	return &_oldRelDisp;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetMaxDispIter()
{
	return &_maxDispIter;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetMinDispIter()
{
	return &_minDispIter;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetUnlDispIter()
{
	return &_unlDispIter;
}

std::vector<std::vector<double*>>* AnalysisSpringRecorder::GetRelDispIter()
{
	return &_relDispIter;
}

std::vector<std::vector<std::string*>>* AnalysisSpringRecorder::GetNewStages()
{
	return &_newStages;
}

std::vector<std::vector<std::string*>>* AnalysisSpringRecorder::GetOldStages()
{
	return &_oldStages;
}

std::vector<std::vector<std::string*>>* AnalysisSpringRecorder::GetListStages() 
{
	return &_listStages;
}

void AnalysisSpringRecorder::SwapOldNewVectors(int ele)
{
	_oldDisp[ele - 1] = _newDisp[ele - 1];
	_oldMinDisp[ele - 1] = _newMinDisp[ele - 1];
	_oldMaxDisp[ele - 1] = _newMaxDisp[ele - 1];
	_oldPlasticDisp[ele - 1] = _newPlasticDisp[ele - 1];
	_oldUnlDisp[ele - 1] = _newUnlDisp[ele - 1];
	_oldRelDisp[ele - 1] = _newRelDisp[ele - 1];
	_oldStages[ele - 1] = _newStages[ele - 1];
}

void AnalysisSpringRecorder::InitializeAllVectors()
{
	//Initializes all disps the vectors
	double* d = new double(0);
	std::string* s = new std::string("initial");
	std::vector<double*> vec(3, d);
	std::vector<std::string*> vec1(3, s);

	for (int i = 0; i < _springsID.size(); i++) { //for each element
		_oldDisp.emplace_back(vec);
		_newDisp.emplace_back(vec);
		_oldMinDisp.emplace_back(vec);
		_newMinDisp.emplace_back(vec);
		_oldMaxDisp.emplace_back(vec);
		_newMaxDisp.emplace_back(vec);
		_oldPlasticDisp.emplace_back(vec);
		_newPlasticDisp.emplace_back(vec);
		_oldUnlDisp.emplace_back(vec);
		_newUnlDisp.emplace_back(vec);
		_oldRelDisp.emplace_back(vec);
		_newRelDisp.emplace_back(vec);
		_maxDispIter.emplace_back(vec);
		_minDispIter.emplace_back(vec);
		_unlDispIter.emplace_back(vec);
		_relDispIter.emplace_back(vec);
		_oldStages.emplace_back(vec1);
		_newStages.emplace_back(vec1);
		_listStages.emplace_back(vec1);

	}
}

void AnalysisSpringRecorder::UpdatePerIterDisps()
{
	_maxDispIter = _newMaxDisp;
	_minDispIter = _newMinDisp;
	_unlDispIter = _newUnlDisp;
	_relDispIter = _newRelDisp;
}
