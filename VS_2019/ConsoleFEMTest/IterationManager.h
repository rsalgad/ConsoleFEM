#pragma once
#include <vector>
#include "Node.h"
#include "ShellElement.h"
#include "Spring3D.h"
#include "Load.h"
#include "SeismicLoad.h"
#include "ImpulseLoad.h"
#include "Mass.h"
#include "Support.h"
#include "TimeIntegrationMethod.h"

class IterationManager
{
public:
	IterationManager();

	static void PerformAnalysisWithIterations(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Spring3D> &listOfSprings, std::vector<Load> &listOfLoads, std::vector<Support> &listOfSups, int nLoadSteps, std::string &fileName);
	//static void PerformAnalysisWithIterationsGeomNonlinear(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Spring3D> &listOfSprings, std::vector<Load> &listOfLoads, std::vector<Support> &listOfSups, int nIter, int nLoadSteps);
	static void PerformAnalysisWithIterationsMatNonlinearDispBased(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Spring3D> &listOfSprings, std::vector<Load> &listOfLoads, std::vector<Support> &listOfSups, int nIter, int nLoadSteps, std::string type, double cyclicRepeat, int stepsPerPeak, double peakInc, int cyclesPerPeak, double iniPeak, std::string &fileName);
	static void PerformDynamicAnalysisWithIterationsMatNonlinearDispBased(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Spring3D> &listOfSprings, std::vector<Load> &listOfLoads, std::vector<Mass> &listOfMasses, std::vector<Support> &listOfSups, SeismicLoad &sLoad, ImpulseLoad &impLoad, int nIter, int nLoadSteps, std::string type, std::string& fileName);
	static void TESTPerformDynamicAnalysisWithIterationsMatNonlinearDispBased();


	~IterationManager();
private:
	static void PostConvergenceProcedures(std::vector<Spring3D> &listOfSprings, std::vector<std::vector<std::string>> &listOfSpringLoadingStages, std::vector<std::vector<std::string>> &newSpringStages,
		Matrix &shellRestricStiff, std::vector<Support> &listOfSups, std::vector<std::vector<double>> &newListOfDisps, std::vector<std::vector<double>> &listOfMinDisps, std::vector<std::vector<double>> &listOfMaxDisps,
		std::vector<std::vector<double>> &newListOfPlasticDisps, std::vector<std::vector<double>> &unlDispPerIter, std::vector<std::vector<double>> &relDispPerIter, Matrix &dNew, std::vector<Node> &listOfNodes,
		Matrix& F_iter, double highStiff, std::vector<Matrix> &forcePerStep, std::vector<Matrix> &velPerStep, std::vector<Matrix> &accPerStep, Matrix &curVel, Matrix &acc, std::vector<std::vector<double>> &maxDispPerIter,
		std::vector<std::vector<double>> &minDispPerIter, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp, Matrix &prevDisp, Matrix &prevVel, Matrix &prevAcc);
	static void CalculateAccelerationAndVelocity(Matrix &dNew, std::vector<std::vector<int>> &totalMassDOFVec, Matrix &redPrevDisp, Matrix &redVel,
		Matrix &redAcc, TimeIntegrationMethod& IntMethod, double& deltaT, Matrix &prevAcc, Matrix &acc, Matrix &curVel, Matrix &prevVel, std::vector<Support> &listOfSups, std::vector<Node> &listOfNodes, Matrix &completeD_prime);
	static bool CheckConvergenceCriteria(Matrix &D, double limit);
};
