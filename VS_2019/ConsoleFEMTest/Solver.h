#pragma once
//#include "pch.h"
#include "Matrix.h"
#include "StructureManager.h"
#include "PreAnalysisSetUp.h"
#include "AnalysisSpringRecorder.h"
#include "Spring3D.h"
#include <map>
#include <vector>

class Solver
{
public:
	Solver();

	//static Matrix CompleteStiffnessMatrixWithThreads(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Spring3D> &listOfSprings, int nThreads);
	static Matrix ReducedStiffnessMatrix(Matrix& shellStiff, const StructureManager* structManager, const PreAnalysisSetUp* setUp, AnalysisSpringRecorder* springRecorder, double* highStiff);
	static Matrix ReducedShellStiffMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static Matrix CompleteShellMassMatrixThreads(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static Matrix ReducedSpringStiffMatrix(const int* redSize, const std::map<int, Spring3D*>* listOfSprings, AnalysisSpringRecorder* springRecorder);
	static Matrix ReducedRestrictStiffnessMatrix(Matrix& shellStiff, const StructureManager* structManager, const PreAnalysisSetUp* setUp, AnalysisSpringRecorder* springRecorder);
	static Matrix ReducedSpringRestrictedStiffMatrix(const PreAnalysisSetUp* setUp, const std::map<int, Spring3D*>* listOfSprings, AnalysisSpringRecorder* springRecorder);
	static void DisplacementLoadStiffness(Matrix& stiff, const std::vector<int>* dispDOFs, double* highStiff);
	static void DisplacementLoadForce(Matrix& force, const StructureManager* structManager, const PreAnalysisSetUp* setUp, AnalysisSpringRecorder* springRecord, const double* loadFraction, const double* highStiff);
	static void PlasticDisplacementLoadForce(Matrix& force, const StructureManager* structManager, const PreAnalysisSetUp* setUp, AnalysisSpringRecorder* springRecord);
	static Matrix ShellRestrictedStiffMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static Matrix CalculateNaturalFrequenciesAndModeShapes(Matrix& stiffMatrix, Matrix& massMatrix, std::vector<double>* natFreq, const std::vector<std::vector<int>>* totalMassDOFVec, const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static std::vector<double> RayleighDampingConstants(const PreAnalysisSetUp* setUp, const std::vector<double>* natFreq);
	static Matrix RayleighDampingMatrix( Matrix& m, Matrix& k, const std::vector<double>* constants);
	static Matrix GetTotalModalMatrix(Matrix& m, const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static bool HasDOFMass(const int* DOF, const StructureManager* structManager, const int* nDOF);
	static Matrix ReducedForceMatrix(Matrix& FConst, Matrix& FIncr, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const int* step, const double* highStiff, AnalysisSpringRecorder* springRecord);
	static Matrix ReducedDynamicStiffnessMatrix(const PreAnalysisSetUp* setUp, Matrix& mass,Matrix& damp, Matrix& add, Matrix& mult2);
	static Matrix CalculateDynamicForce(Matrix& prevDisp, Matrix& prevVel, Matrix& prevAcc, Matrix& add, Matrix& damp, Matrix& totMass, Matrix& FInc, Matrix& mRed, Matrix& add11, const PreAnalysisSetUp* setUp, const double* time);

	~Solver();
};

