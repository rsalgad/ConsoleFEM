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
	void SetDisplacementValues(std::string type, std::string status, int springID, std::vector<double*> newDisps);
	void SetIndividualDisplacementValues(std::string type, std::string status, int springID, int pos, double* val);
	void SwapOldNewDisps(std::string type, int springID);
	void SwapOldNewStages(int springID);
	//std::map<std::string, std::map<std::string, std::map<int, std::vector<double>>>>* GetDisplacementMap();
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>> GetDisplacementMap() const;
	double* GetDisplacement(std::string type, std::string status, int springID, int pos) const;
	std::vector<double*> GetDisplacementVector(std::string type, std::string status, int springID) const;

	void SetStagesValues(std::string flag, int springID, std::vector<std::string*> newDisps);
	void SetIndividualStagesValues(std::string flag, int springID, int pos, std::string* string);
	std::map<std::string, std::map<int, std::vector<std::string*>>> GetStagesMap() const;
	std::string* GetStages(std::string flag, int springID, int pos) const;
	std::vector<std::string*> GetStagesVector(std::string flag, int springID) const;
	std::map<std::string, std::map<int, std::vector<double*>>> GetDisplacementPerIterMap() const;
	double* GetDisplacementPerIter(std::string status, int springID, int pos) const;

	void InitializeAllVectors();
	void UpdatePerIterDisps();

private:
	//Hierarchy-> First map keys: "disp", "minDisp", "maxDisp", "plasticDisp", "unlDisp", "relDisp", "maxDispIter", "minDispIter", "unlDispIter", "relDispIter".
	//second map keys: "old", "new"; third map: "spring_ID"
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double*>>>> _disps;

	//"maxDispIter", "minDispIter", "unlDispIter", "relDispIter"
	std::map<std::string, std::map<int, std::vector<double*>>> _perIterDisps;

	//Hierarchy-> First map keys: "old", "new", "list"; second map keys: "spring_ID"
	std::map<std::string, std::map<int, std::vector<std::string*>>> _stages;

	std::vector<int> _springsID;
};

