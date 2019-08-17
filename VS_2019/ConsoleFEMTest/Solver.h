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
	static Matrix ReducedStiffnessMatrix(const Matrix* shellStiff, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecorder, double* highStiff);
	static Matrix ReducedShellStiffMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static Matrix CompleteShellMassMatrixThreads(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static Matrix ReducedSpringStiffMatrix(const int* redSize, const std::map<int, Spring3D*>* listOfSprings, const AnalysisSpringRecorder* springRecorder);
	static Matrix ReducedRestrictStiffnessMatrix(const Matrix* shellStiff, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecorder);
	static Matrix ReducedSpringRestrictedStiffMatrix(const int* redSize, const std::map<int, Spring3D*>* listOfSprings, const AnalysisSpringRecorder* springRecorder);
	static void DisplacementLoadStiffness(Matrix& stiff, const std::vector<int>* dispDOFs, double* highStiff);
	static void DisplacementLoadForce(const Matrix* force, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecord, const double* loadFraction, const double* highStiff);
	static void PlasticDisplacementLoadForce(const Matrix* force, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecord);
	static Matrix ShellRestrictedStiffMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static Matrix CalculateNaturalFrequenciesAndModeShapes(const Matrix* stiffMatrix, const Matrix* massMatrix, std::vector<double>* natFreq, const std::vector<std::vector<int>>* totalMassDOFVec, const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static std::vector<double> RayleighDampingConstants(const PreAnalysisSetUp* setUp, const std::vector<double>* natFreq);
	static Matrix RayleighDampingMatrix(const Matrix* m, const Matrix* k, const std::vector<double>* constants);
	static Matrix GetTotalModalMatrix(const Matrix* m, const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static bool HasDOFMass(const int* DOF, const StructureManager* structManager, const int* nDOF);
	static Matrix ReducedForceMatrix(const Matrix* FConst, const Matrix* FIncr, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const int* step, const double* highStiff, const AnalysisSpringRecorder* springRecord);
	static Matrix ReducedDynamicStiffnessMatrix(const PreAnalysisSetUp* setUp, const Matrix* mass, const Matrix* damp, Matrix* add, Matrix* mult2);
	static Matrix CalculateDynamicForce(const Matrix* prevDisp, const Matrix* prevVel, const Matrix* prevAcc, Matrix* add, const Matrix* damp, const Matrix* totMass, const Matrix* FInc, const Matrix* mRed, Matrix* add11, const PreAnalysisSetUp* setUp, const double* time);

	~Solver();
};

