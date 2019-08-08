// ConsoleFEMTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <time.h>
#include <string>
#include <vector>
#include <fstream>
#include <array>
#include "Matrix.h"
#include "Node.h"
#include "ElasticMaterial.h"
#include "OrthotropicElasticMaterial.h"
#include "Load.h"
#include "Mass.h"
#include "Support.h"
#include "Displacement.h"
#include "ShellElement.h"
#include "Spring3D.h"
#include "ImpulseLoad.h"
#include <mutex>
#include <thread>
#include <filesystem>
#include "Solver.h"
#include "SeismicLoad.h"
#include "IterationManager.h"
#include "SpringAxialModel.h"
#include "SpringGeneralModel.h"
#include "FileOperation.h"
#include "StructureManager.h"
#include <Eigen/Eigenvalues> // header file

StructureManager myStrucManager;
std::vector<Node> listOfNodes;
std::vector<MaterialModel*> listOfMaterials;
std::vector<ElasticMaterial> listOfElasticShellMaterials;
std::vector<OrthotropicElasticMaterial> listOfShellMaterials;
std::vector<SpringAxialModel> listOfSpringAxialMat;
std::vector<SpringGeneralModel> listOfSpringGeneralMat;
std::vector<ShellElement> listOfShellElements;
std::vector<Spring3D> listOfSpringElements;
std::vector<Load> listOfLoads;
std::vector<Mass> listOfMasses;
std::vector<Support> listOfSupports;
SeismicLoad seismicLoad;
ImpulseLoad impulseLoad;

int main(int argc, char *argv[])
{
	std::string fileName;
	if (argc != 1) {

		fileName = argv[1];
	}
	else {
		//fileName = "OneEleCyclicTest2.xml";
		fileName = "OneEleImpTestPanelSpring";
	}
	
	//FileOperation::SaveResultsFile();

	//IterationManager::TESTPerformDynamicAnalysisWithIterationsMatNonlinearDispBased();

	FileOperation::ReadInputFromXML(fileName, listOfMaterials, listOfShellMaterials, listOfSpringAxialMat, listOfSpringGeneralMat, listOfLoads, listOfSupports, listOfNodes, listOfShellElements, listOfSpringElements, listOfMasses, seismicLoad, impulseLoad, myStrucManager);

	//Support::SortByNodeID(listOfSupports);
	//myStrucManager.SortSupportsByNodeID();
	
	//Load::SortByNodeID(listOfLoads); //I don't think that's needed

	std::vector<int> index = Load::IdentifyIncrementalLoads(listOfLoads); //Why is this here?

	/*
	if (listOfShellElements.size() != 0) {
		for (int i = 0; i < listOfShellElements.size(); i++) {
			listOfShellElements[i].CalculateGlobalDOFVector(listOfSupports);
			listOfShellElements[i].CalculateGlobalMassDOFVector(listOfMasses, listOfSupports);
		}
	}
	if (listOfSpringElements.size() != 0) {
		for (int i = 0; i < listOfSpringElements.size(); i++) {
			listOfSpringElements[i].CalculateGlobalDOFVector(listOfSupports);
		}
	}
	*/

	clock_t t;
	t = clock();
	//IterationManager::PerformAnalysisWithIterations(listOfNodes, listOfShellElements, listOfSpringElements, listOfLoads, listOfSupports, 10, fileName);
	//IterationManager::PerformAnalysisWithIterationsMatNonlinearDispBased(listOfNodes, listOfShellElements, listOfSpringElements, listOfLoads, listOfSupports, 100, 10, "reverse-cyclic", 6, 4, 1, 1, 0.5, fileName);
	IterationManager::PerformDynamicAnalysisWithIterationsMatNonlinearDispBased(listOfNodes, listOfShellElements, listOfSpringElements, listOfLoads, listOfMasses, listOfSupports, seismicLoad, impulseLoad, 100, 10, "impulse", fileName);
	//IterationManager::PerformAnalysisWithIterationsGeomNonlinear(listOfNodes, listOfShellElements, listOfSpringElements, listOfLoads, listOfSupports, 10, 10);
	t = clock() - t;
	std::cout << "Time = " << t << std::endl;

	return 0;
}