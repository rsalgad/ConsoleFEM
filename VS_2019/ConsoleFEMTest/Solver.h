#pragma once
#include "Matrix.h"
#include <vector>
#include "ShellElement.h"
#include "Spring3D.h"
#include "Node.h"
#include "Mass.h"
#include "StructureManager.h"
#include "PreAnalysisSetUp.h"
#include "AnalysisSpringRecorder.h"

class Solver
{
public:
	Solver();

	//static Matrix CompleteStiffnessMatrixWithThreads(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Spring3D> &listOfSprings, int nThreads);
	static Matrix ReducedStiffnessMatrix(const Matrix* shellStiff, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecorder, double* highStiff);
	static Matrix ReducedShellStiffMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static Matrix CompleteShellMassMatrixThreads(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static Matrix ReducedSpringStiffMatrix(const int* redSize, const std::map<int, Spring3D*>* listOfSprings, const AnalysisSpringRecorder* springRecorder);
	static Matrix ReducedRestrictStiffnessMatrix(const Matrix* shellStiff, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecorder);
	static Matrix ReducedSpringRestrictedStiffMatrix(const int* redSize, const std::map<int, Spring3D*>* listOfSprings, const AnalysisSpringRecorder* springRecorder);
	static void DisplacementLoadStiffness(Matrix& stiff, const std::vector<int>* dispDOFs, double* highStiff);
	static void DisplacementLoadForce(const Matrix* force, const std::map<int, Support*>* listOfSup, const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecord, const double* loadFraction, const double* highStiff);
	static void PlasticDisplacementLoadForce(Matrix& force, std::vector<Support> &listOfSup, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<Spring3D> &listOfSpring, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp);
	static Matrix ShellRestrictedStiffMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static void CalculateNaturalFrequenciesAndModeShapes(Matrix &stiffMatrix, Matrix &massMatrix, std::vector<double> &natFreq, Matrix &modeShapes, std::vector<std::vector<int>> &totalMassDOFVec);
	static std::vector<double> RayleighDampingConstants(int mode1, double damp1, int mode2, double damp2, std::vector<double> &natFreq);
	static Matrix RayleighDampingMatrix(Matrix &m, Matrix &k, std::vector<double>& constants);
	static Matrix GetTotalModalMatrix(Matrix &m, std::vector<Support> &vecSup, std::vector<Node> &vecNode, std::vector<ShellElement> &vecShell, std::vector<Mass> &vecMass);
	static bool HasDOFMass(std::vector<ShellElement> &shellVec, std::vector<Support> &supVec, int DOF);
	static Matrix ReducedForceMatrix(const Matrix* FConst, const Matrix* FIncr, const std::map<int, Support*>* listOfSups, const PreAnalysisSetUp* setUp, const int* step, const double* highStiff, const AnalysisSpringRecorder* springRecord);
	
	~Solver();
};

