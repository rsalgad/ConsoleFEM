#pragma once
#include "Spring3D.h"
#include <map>

class AnalysisSpringRecorder
{
public:
	AnalysisSpringRecorder();
	AnalysisSpringRecorder(std::vector<int> springsID);
	AnalysisSpringRecorder(char c, const std::map<int, Spring3D*>* allNodes);
	~AnalysisSpringRecorder();
	void SetDisplacementValues(std::string type, std::string status, int springID, std::vector<double> newDisps);
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double>>>> GetDisplacementMap() const;

	void SetStagesValues(std::string flag, int springID, std::vector<std::string> newDisps);
	std::map<std::string, std::map<int, std::vector<std::string>>> GetStagesMap() const;
	std::map<std::string, std::map<int, std::vector<double>>> GetDisplacementPerIterMap() const;

	void InitializeAllVectors();
	void UpdatePerIterDisps();

private:
	//Hierarchy-> First map keys: "disp", "minDisp", "maxDisp", "plasticDisp", "unlDisp", "relDisp", "maxDispIter", "minDispIter", "unlDispIter", "relDispIter".
	//second map keys: "old", "new"; third map: "spring_ID"
	std::map<std::string, std::map<std::string, std::map<int, std::vector<double>>>> _disps;
	
	std::map < std::string, std::map<int, std::vector<double>>> _perIterDisps;

	//Hierarchy-> First map keys: "old", "new", "list"; second map keys: "spring_ID"
	std::map<std::string, std::map<int, std::vector<std::string>>> _stages;

	std::vector<int> _springsID;
};

