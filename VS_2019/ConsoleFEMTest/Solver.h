#pragma once
#include "Matrix.h"
#include <vector>
#include "ShellElement.h"
#include "Spring3D.h"
#include "Node.h"
#include "Mass.h"
#include "StructureManager.h"
#include "PreAnalysisSetUp.h"

class Solver
{
public:
	Solver();

	//static Matrix CompleteStiffnessMatrixWithThreads(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Spring3D> &listOfSprings, int nThreads);
	static Matrix ReducedStiffnessMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static Matrix ReducedShellStiffMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static Matrix CompleteShellMassMatrixThreads(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Support> &listOfSupports, int nThreads, std::vector<std::vector<ShellElement>> &shellElemVecs);
	static void CompleteSpringStiffMatrixThreadsDispBasedAfterShells(std::vector<Spring3D> &listOfSprings, Matrix &m, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp);
	static void CompleteSpringRestrictedStiffMatrixThreadsDispBasedAfterShells(std::vector<Spring3D> &listOfSprings, Matrix &m, std::vector<Support> &sup, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp);
	static void DisplacementLoadStiffness(Matrix& stiff, std::vector<Support> &listOfSup, double &biggest);
	static void DisplacementLoadForce(Matrix& force, std::vector<Support> &listOfSup, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<Spring3D> &listOfSpring, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, double &biggest, double loadFraction, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp);
	static void PlasticDisplacementLoadForce(Matrix& force, std::vector<Support> &listOfSup, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<Spring3D> &listOfSpring, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp);
	static Matrix ShellRestrictedStiffMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static void CalculateNaturalFrequenciesAndModeShapes(Matrix &stiffMatrix, Matrix &massMatrix, std::vector<double> &natFreq, Matrix &modeShapes, std::vector<std::vector<int>> &totalMassDOFVec);
	static std::vector<double> RayleighDampingConstants(int mode1, double damp1, int mode2, double damp2, std::vector<double> &natFreq);
	static Matrix RayleighDampingMatrix(Matrix &m, Matrix &k, std::vector<double>& constants);
	static Matrix GetTotalModalMatrix(Matrix &m, std::vector<Support> &vecSup, std::vector<Node> &vecNode, std::vector<ShellElement> &vecShell, std::vector<Mass> &vecMass);
	static bool HasDOFMass(std::vector<ShellElement> &shellVec, std::vector<Support> &supVec, int DOF);
	~Solver();
};

