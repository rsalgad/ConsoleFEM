#pragma once
#include "MaterialModel.h"
#include "ElasticMaterial.h"
#include "OrthotropicElasticMaterial.h"
#include "SpringAxialModel.h"
#include "SpringGeneralModel.h"
#include "Load.h"
#include "SeismicLoad.h"
#include "ImpulseLoad.h"
#include "Mass.h"
#include "Support.h"
#include "Node.h"
#include "ShellElement.h"
#include "Spring3D.h"
#include "StructureManager.h"
#include <string>

class FileOperation
{
public:
	FileOperation();
	static void ReadInputFromXML(std::string fileName, std::vector<MaterialModel*> &listOfMaterials, std::vector<OrthotropicElasticMaterial> &listOfShellMaterials, std::vector<SpringAxialModel> &listOfSpringAxialMat, std::vector<SpringGeneralModel> &listOfSpringGeneralMat, std::vector<Load> &listOfLoads, std::vector<Support> &listOfSupports, std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShellElements, std::vector<Spring3D> &listOfSpringElements, std::vector<Mass> &listOfMasses, SeismicLoad &sLoad, ImpulseLoad &impLoad, StructureManager &structManager);
	static void SaveIterationsResult(std::string fileName, std::vector<std::vector<Node>> result, std::vector<Node> originalList);
	static void SaveIterationsForceResult(std::string fileName, std::vector<Matrix> result);
	static void SaveResultsFile(std::string &fileName, std::vector<std::vector<Node>> &nodePerStep, std::vector<Node> &listOfNodes, std::vector<Matrix> &forcePerStep, std::vector<double> &natFreq, Matrix &modeShape);
	~FileOperation();
};

