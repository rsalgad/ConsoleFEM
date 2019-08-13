#include "pch.h"
#include "Spring3D.h"
#include "AnalysisSpringRecorder.h"

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

void AnalysisSpringRecorder::SetDisplacementValues(std::string type, std::string status, int springID, std::vector<double> newDisps)
{
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double>>>>::iterator it = _disps.find(type);
	std::map<std::string, std::map<int, std::vector<double>>>::iterator it2 = it->second.find(status);
	std::map<int, std::vector<double>>::iterator it3 = it2->second.find(springID);
	if (it3 != it2->second.end()) {
		it3->second = newDisps;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}
}

std::map<std::string, std::map<std::string, std::map<int, std::vector<double>>>> AnalysisSpringRecorder::GetDisplacementMap() const
{
	return _disps;
}

std::map<std::string, std::map<int, std::vector<double>>> AnalysisSpringRecorder::GetDisplacementPerIterMap() const
{
	return _perIterDisps;
}

void AnalysisSpringRecorder::SetStagesValues(std::string flag, int springID, std::vector<std::string> newDisps)
{
	std::map<std::string, std::map<int, std::vector<std::string>>>::iterator it = _stages.find(flag);
	std::map<int, std::vector<std::string>>::iterator it2 = it->second.find(springID);
	if (it2 != it->second.end()) {
		it2->second = newDisps;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}
}

std::map<std::string, std::map<int, std::vector<std::string>>> AnalysisSpringRecorder::GetStagesMap() const
{
	return _stages;
}

void AnalysisSpringRecorder::InitializeAllVectors()
{
	//Initializes all disps the vectors
	std::vector<double> vec(3, 0);
	std::map<int, std::vector<double>> map1;

	for (int i = 0; i < _springsID.size(); i++) { //for each spring ID
		map1.insert(std::pair<int, std::vector<double>>(_springsID[i], vec));
	}
	std::map<std::string, std::map<int, std::vector<double>>> map2;
	map2.insert(std::pair<std::string, std::map<int, std::vector<double>>>("old", map1));
	map2.insert(std::pair<std::string, std::map<int, std::vector<double>>>("new", map1));

	std::map<std::string, std::map<std::string, std::map<int, std::vector<double>>>> map3;
	map3.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double>>>>("disp", map2));
	map3.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double>>>>("minDisp", map2));
	map3.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double>>>>("maxDisp", map2));
	map3.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double>>>>("plasticDisp", map2));
	map3.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double>>>>("unlDisp", map2));
	map3.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double>>>>("relDisp", map2));

	std::map<std::string, std::map<int, std::vector<double>>> map4;
	map4.insert(std::pair<std::string, std::map<int, std::vector<double>>>("maxDispIter", map1));
	map4.insert(std::pair<std::string, std::map<int, std::vector<double>>>("minDispIter", map1));
	map4.insert(std::pair<std::string, std::map<int, std::vector<double>>>("unlDispIter", map1));
	map4.insert(std::pair<std::string, std::map<int, std::vector<double>>>("relDispIter", map1));

	//Initializes all stages the vectors
	std::vector<std::string> vec1(3, "initial");
	std::map<int, std::vector<std::string>> map5;

	for (int i = 0; i < _springsID.size(); i++) { //for each spring ID
		map5.insert(std::pair<int, std::vector<std::string>>(_springsID[i], vec1));
	}

	std::map<std::string, std::map<int, std::vector<std::string>>> map6;
	map6.insert(std::pair<std::string, std::map<int, std::vector<std::string>>>("old", map5));
	map6.insert(std::pair<std::string, std::map<int, std::vector<std::string>>>("new", map5));
	map6.insert(std::pair<std::string, std::map<int, std::vector<std::string>>>("list", map5));
}

void AnalysisSpringRecorder::UpdatePerIterDisps()
{
	_perIterDisps.find("maxDisp")->second = _disps.find("maxDisp")->second.find("new")->second;
	_perIterDisps.find("minDisp")->second = _disps.find("minDisp")->second.find("new")->second;
	_perIterDisps.find("unlDisp")->second = _disps.find("unlDisp")->second.find("new")->second;
	_perIterDisps.find("relDisp")->second = _disps.find("relDisp")->second.find("new")->second;
}
