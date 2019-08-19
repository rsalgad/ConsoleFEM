#pragma once
//#include "pch.h"
#include <string>
#include <map>
#include <vector>

class Spring3D;

class AnalysisSpringRecorder
{
public:
	AnalysisSpringRecorder();
	AnalysisSpringRecorder(std::vector<int> springsID);
	AnalysisSpringRecorder(char c, const std::map<int, Spring3D*>* allSprings);
	~AnalysisSpringRecorder();

	std::vector<std::vector<double*>>* GetNewDisp();
	//std::vector<std::vector<double*>>* GetNewDisp() const;
	std::vector<std::vector<double*>>* GetOldDisp();
	std::vector<std::vector<double*>>* GetNewMinDisp();
	std::vector<std::vector<double*>>* GetOldMinDisp();
	std::vector<std::vector<double*>>* GetNewMaxDisp();
	std::vector<std::vector<double*>>* GetOldMaxDisp();
	std::vector<std::vector<double*>>* GetNewPlasticDisp();
	std::vector<std::vector<double*>>* GetOldPlasticDisp();
	std::vector<std::vector<double*>>* GetNewUnlDisp();
	std::vector<std::vector<double*>>* GetOldUnlDisp();
	std::vector<std::vector<double*>>* GetNewRelDisp();
	std::vector<std::vector<double*>>* GetOldRelDisp();
	std::vector<std::vector<double*>>* GetMaxDispIter();
	std::vector<std::vector<double*>>* GetMinDispIter();
	std::vector<std::vector<double*>>* GetUnlDispIter();
	std::vector<std::vector<double*>>* GetRelDispIter();
	std::vector<std::vector<std::string*>>* GetNewStages();
	std::vector<std::vector<std::string*>>* GetOldStages();
	std::vector<std::vector<std::string*>>* GetListStages();

	void SwapOldNewVectors(int ele);
	void InitializeAllVectors();
	void UpdatePerIterDisps();

private:
	std::vector<std::vector<double*>> _newDisp;
	std::vector<std::vector<double*>> _oldDisp;
	std::vector<std::vector<double*>> _newMinDisp;
	std::vector<std::vector<double*>> _oldMinDisp;
	std::vector<std::vector<double*>> _newMaxDisp;
	std::vector<std::vector<double*>> _oldMaxDisp;
	std::vector<std::vector<double*>> _newPlasticDisp;
	std::vector<std::vector<double*>> _oldPlasticDisp;
	std::vector<std::vector<double*>> _newUnlDisp;
	std::vector<std::vector<double*>> _oldUnlDisp;
	std::vector<std::vector<double*>> _newRelDisp;
	std::vector<std::vector<double*>> _oldRelDisp;


	//"maxDispIter", "minDispIter", "unlDispIter", "relDispIter"
	std::vector<std::vector<double*>> _maxDispIter;
	std::vector<std::vector<double*>> _minDispIter;
	std::vector<std::vector<double*>> _unlDispIter;
	std::vector<std::vector<double*>> _relDispIter;

	//Hierarchy-> First map keys: "old", "new", "list"; second map keys: "spring_ID"
	std::vector<std::vector<std::string*>> _oldStages;
	std::vector<std::vector<std::string*>> _newStages;
	std::vector<std::vector<std::string*>> _listStages;

	std::vector<int> _springsID;
};

