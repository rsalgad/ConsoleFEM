#pragma once
//#include "pch.h"
#include <string>
#include <vector>
#include <map>
#include "OrthotropicElasticMaterial.h"
#include "SpringAxialModel.h"
#include "SpringGeneralModel.h"
#include "Load.h"
#include "Support.h"
#include "Node.h"
#include "ShellElement.h"
#include "Spring3D.h"
#include "Mass.h"
#include "SeismicLoad.h"
#include "ImpulseLoad.h"
#include "MaterialModel.h"
#include "StructureManager.h"
#include "Matrix.h"
#include "NodalRecorder.h"

class FileOperation
{
public:
	FileOperation();
	static void ReadInputFromXML(std::string fileName, StructureManager& structManager, AnalysisMethod* *analysis);
	static void SaveIterationsResult(std::string fileName, const NodalRecorder<Node>* disps, const std::map<int, Node*>* nodes);
	//static void SaveIterationsForceResult(std::string fileName, const std::map<int, Load*>* result);
	static void SaveResultsFile(std::string &fileName, std::vector<std::vector<Node>> &nodePerStep, std::vector<Node> &listOfNodes, std::vector<Matrix> &forcePerStep, std::vector<double> &natFreq, Matrix &modeShape);
	~FileOperation();
};

