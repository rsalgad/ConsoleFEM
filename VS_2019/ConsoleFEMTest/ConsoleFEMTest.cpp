// ConsoleFEMTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include "FileOperation.h"
#include "IterationManager.h"
#include <map>
#include <time.h>
#include <fstream>
#include <array>
#include <mutex>
#include <thread>
#include <filesystem>
#include <Eigen/Eigenvalues> // header file

StructureManager myStrucManager;
AnalysisMethod* myAnalysisMethod = nullptr;
SeismicLoad seismicLoad;
ImpulseLoad impulseLoad;

int main(int argc, char *argv[])
{
	std::string fileName;
	bool hasFile = false;

	/*
	double sum = 0;
	clock_t t1;
	t1= clock();
	
	#pragma omp parallel for reduction(+:sum)
		for (int i = 0; i < 10000000000; i++) {
			sum += 1234 * 1234;
		}

		/*
	for (int i = 0; i < 1000000000; i++) {
		sum += 1234 * 1234;
	}
		
	t1 = clock() - t1;
	printf("The sum is %f\n", sum);
	std::cout << "Time " << t1 << std::endl;
	*/
	if (argc != 1) {
		fileName = argv[1];
		hasFile = true;
	}
	else {
		std::cout << "No Input File Specified." << std::endl;
		fileName = "BigStructure_Impulse3";
	}
	
	if (true) { //input file specified

		FileOperation::ReadInputFromXML(fileName, myStrucManager, &myAnalysisMethod);
		PreAnalysisSetUp setUp = PreAnalysisSetUp(&myStrucManager, myAnalysisMethod);

		//Support::SortByNodeID(listOfSupports);
		//myStrucManager.SortSupportsByNodeID();

		//Load::SortByNodeID(listOfLoads); //I don't think that's needed

		//std::vector<int> index = Load::IdentifyIncrementalLoads(listOfLoads); //Why is this here?

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
		IterationManager::PerformDynamicAnalysisWithIterationsMatNonlinearDispBased(&myStrucManager, &setUp, fileName);
		//IterationManager::PerformAnalysisWithIterationsGeomNonlinear(listOfNodes, listOfShellElements, listOfSpringElements, listOfLoads, listOfSupports, 10, 10);
		t = clock() - t;
		std::cout << "Time = " << t << std::endl;
	}

	return 0;
}