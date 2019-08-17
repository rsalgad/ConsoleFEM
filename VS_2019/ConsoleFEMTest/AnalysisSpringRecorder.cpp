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

void AnalysisSpringRecorder::SetDisplacementValues(std::string type, std::string status, int springID, std::vector<double*> newDisps)
{
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>::iterator it = _disps.find(type);
	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it2 = it->second.find(status);
	std::map<int, std::vector<double*>>::iterator it3 = it2->second.find(springID);
	if (it3 != it2->second.end()) {
		it3->second = newDisps;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}
}

void AnalysisSpringRecorder::SetIndividualDisplacementValues(std::string type, std::string status, int springID, int pos, double* val)
{
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>::iterator it = _disps.find(type);
	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it2 = it->second.find(status);
	std::map<int, std::vector<double*>>::iterator it3 = it2->second.find(springID);
	if (it3 != it2->second.end()) {
		it3->second[pos] = val;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}
}

void AnalysisSpringRecorder::SwapOldNewDisps(std::string type, int springID)
{
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>::iterator it = _disps.find(type);
	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it2 = it->second.find("old");
	std::map<int, std::vector<double*>>::iterator it3 = it2->second.find(springID);
	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it4 = it->second.find("new");
	std::map<int, std::vector<double*>>::iterator it5 = it4->second.find(springID);
	if (it3 != it2->second.end()) {
		it3->second = it5->second;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}
}

void AnalysisSpringRecorder::SwapOldNewStages(int springID)
{
	std::map<std::string, std::map<int, std::vector<std::string*>>>::iterator it = _stages.find("old");
	std::map<int, std::vector<std::string*>>::iterator it2 = it->second.find(springID);
	std::map<std::string, std::map<int, std::vector<std::string*>>>::iterator it3 = _stages.find("new");
	std::map<int, std::vector<std::string*>>::iterator it4 = it3->second.find(springID);
	if (it2 != it->second.end()) {
		it2->second = it4->second;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}
	
}

std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>> AnalysisSpringRecorder::GetDisplacementMap() const
{
	return _disps;
}

double* AnalysisSpringRecorder::GetDisplacement(std::string type, std::string status, int springID, int pos) const
{
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>::const_iterator it = _disps.find(type);
	std::map<std::string, std::map<int, std::vector<double*>>>::const_iterator it2 = it->second.find(status);
	std::map<int, std::vector<double*>>::const_iterator it3 = it2->second.find(springID);
	if (it3 != it2->second.end()) {
		return it3->second[pos];
	}
	else {
		std::cout << "Error setting value" << std::endl;
		return NULL;
	}
}

std::vector<double*> AnalysisSpringRecorder::GetDisplacementVector(std::string type, std::string status, int springID) const
{
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>::const_iterator it = _disps.find(type);
	std::map<std::string, std::map<int, std::vector<double*>>>::const_iterator it2 = it->second.find(status);
	std::map<int, std::vector<double*>>::const_iterator it3 = it2->second.find(springID);
	if (it3 != it2->second.end()) {
		return it3->second;
	}
	else {
		std::cout << "Error setting value" << std::endl;
		return std::vector<double*>();
	}
}

std::map<std::string, std::map<int, std::vector<double*>>> AnalysisSpringRecorder::GetDisplacementPerIterMap() const
{
	return _perIterDisps;
}

double* AnalysisSpringRecorder::GetDisplacementPerIter(std::string status, int springID, int pos) const
{
	std::map<std::string, std::map<int, std::vector<double*>>>::const_iterator it = _perIterDisps.find(status);
	std::map<int, std::vector<double*>>::const_iterator it2 = it->second.find(springID);
	if (it2 != it->second.end()) {
		return it2->second[pos];
	}
	else {
		std::cout << "Error setting value" << std::endl;
		return NULL;
	}
}

void AnalysisSpringRecorder::SetStagesValues(std::string flag, int springID, std::vector<std::string*> newDisps)
{
	std::map<std::string, std::map<int, std::vector<std::string*>>>::iterator it = _stages.find(flag);
	std::map<int, std::vector<std::string*>>::iterator it2 = it->second.find(springID);
	if (it2 != it->second.end()) {
		it2->second = newDisps;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}
}

void AnalysisSpringRecorder::SetIndividualStagesValues(std::string flag, int springID, int pos, std::string* string)
{
	std::map<std::string, std::map<int, std::vector<std::string*>>>::iterator it = _stages.find(flag);
	std::map<int, std::vector<std::string*>>::iterator it2 = it->second.find(springID);
	if (it2 != it->second.end()) {
		it2->second[pos] = string;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}
}

std::map<std::string, std::map<int, std::vector<std::string*>>> AnalysisSpringRecorder::GetStagesMap() const
{
	return _stages;
}

std::string* AnalysisSpringRecorder::GetStages(std::string flag, int springID, int pos) const
{
	std::map<std::string, std::map<int, std::vector<std::string*>>>::const_iterator it = _stages.find(flag);
	std::map<int, std::vector<std::string*>>::const_iterator it2 = it->second.find(springID);
	if (it2 != it->second.end()) {
		return it2->second[pos];
	}
	else {
		std::cout << "Error setting value" << std::endl;
		return NULL;
	}
}

std::vector<std::string*> AnalysisSpringRecorder::GetStagesVector(std::string flag, int springID) const
{
	std::map<std::string, std::map<int, std::vector<std::string*>>>::const_iterator it = _stages.find(flag);
	std::map<int, std::vector<std::string*>>::const_iterator it2 = it->second.find(springID);
	if (it2 != it->second.end()) {
		return it2->second;
	}
	else {
		std::cout << "Error setting value" << std::endl;
		return std::vector<std::string*>();
	}

}

void AnalysisSpringRecorder::InitializeAllVectors()
{
	//Initializes all disps the vectors
	std::vector<double*> vec(3, new double(0));
	std::map<int, std::vector<double*>> map1;

	for (int i = 0; i < _springsID.size(); i++) { //for each spring ID
		map1.insert(std::pair<int, std::vector<double*>>(_springsID[i], vec));
	}
	std::map<std::string, std::map<int, std::vector<double*>>> map2;
	map2.insert(std::pair<std::string, std::map<int, std::vector<double*>>>("old", map1));
	map2.insert(std::pair<std::string, std::map<int, std::vector<double*>>>("new", map1));

	_disps.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>("disp", map2));
	_disps.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>("minDisp", map2));
	_disps.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>("maxDisp", map2));
	_disps.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>("plasticDisp", map2));
	_disps.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>("unlDisp", map2));
	_disps.insert(std::pair<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>("relDisp", map2));

	_perIterDisps.insert(std::pair<std::string, std::map<int, std::vector<double*>>>("maxDispIter", map1));
	_perIterDisps.insert(std::pair<std::string, std::map<int, std::vector<double*>>>("minDispIter", map1));
	_perIterDisps.insert(std::pair<std::string, std::map<int, std::vector<double*>>>("unlDispIter", map1));
	_perIterDisps.insert(std::pair<std::string, std::map<int, std::vector<double*>>>("relDispIter", map1));

	//Initializes all stages the vectors
	std::vector<std::string*> vec1(3, new std::string("initial"));
	std::map<int, std::vector<std::string*>> map5;

	for (int i = 0; i < _springsID.size(); i++) { //for each spring ID
		map5.insert(std::pair<int, std::vector<std::string*>>(_springsID[i], vec1));
	}

	_stages.insert(std::pair<std::string, std::map<int, std::vector<std::string*>>>("old", map5));
	_stages.insert(std::pair<std::string, std::map<int, std::vector<std::string*>>>("new", map5));
	_stages.insert(std::pair<std::string, std::map<int, std::vector<std::string*>>>("list", map5));
}

void AnalysisSpringRecorder::UpdatePerIterDisps()
{
	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it = _perIterDisps.find("maxDispIter");
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>::iterator it2 = _disps.find("maxDisp");
	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it3 = it2->second.find("new");

	if (it != _perIterDisps.end()) {
		it->second = it3->second;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}

	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it4 = _perIterDisps.find("minDispIter");
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>::iterator it5 = _disps.find("minDisp");
	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it6 = it5->second.find("new");

	if (it4 != _perIterDisps.end()) {
		it4->second = it6->second;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}

	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it7 = _perIterDisps.find("unlDispIter");
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>::iterator it8 = _disps.find("unlDisp");
	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it9 = it8->second.find("new");

	if (it7 != _perIterDisps.end()) {
		it7->second = it9->second;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}

	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it10 = _perIterDisps.find("relDispIter");
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>>::iterator it11 = _disps.find("relDisp");
	std::map<std::string, std::map<int, std::vector<double*>>>::iterator it12 = it11->second.find("new");

	if (it10 != _perIterDisps.end()) {
		it10->second = it12->second;
	}
	else {
		std::cout << "Error setting value" << std::endl;
	}
}
