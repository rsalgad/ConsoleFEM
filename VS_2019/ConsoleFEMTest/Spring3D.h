#pragma once
//#include "pch.h"
#include "StructureManager.h"
#include "Node.h"
#include "PreAnalysisSetUp.h"
#include "SpringMaterialModels.h"
#include "AnalysisSpringRecorder.h"
#include "Support.h"
#include "Matrix.h"

#include <vector>
#include <map>
#include <string>
#include <mutex>

class Spring3D
{
public:
	Spring3D();
	Spring3D(int ID, Node* n1, Node* n2, std::vector<SpringMaterialModels*> listOfMaterials, char axialDir, char shearDir);
	std::string ToString() const;
	int GetID() const;
	Node GetNode1();
	Node GetNode2();
	std::vector<SpringMaterialModels*> GetListOfMaterials();
	std::vector<int> GetListOfGlobalMaterialDirections();
	std::vector<std::vector<int>> GetGlobalDOFVector();
	std::vector<std::vector<int>> GetGlobalRestDOFVector();
	Matrix GetLocalStifnessMatrixDispBased(std::vector<double> &listOfPos, std::vector<double> &listOfMinDisp, std::vector<double> &listOfMaxDisp, std::vector<double> &listOfPlasticDisp, std::vector<std::string> &listOfLoadStage, std::vector<std::string> &listOfStage, std::vector<double> &listOfUnlDisp, std::vector<double> &listOfRelDisp);
	Matrix GetTransformationMatrix();
	Matrix GetVectorBetweenNodes();
	Matrix GetGlobalStiffMatrixDispBased(const AnalysisSpringRecorder* springRecorder);
	std::vector<double> GlobalDOFVectorNotUsed();
	void CalculateGlobalDOFVector(const std::map<int, Support*>* supList, const int* DOF);
	Matrix GetElementGlobalDisplacementVector(Matrix completeD);
	Matrix GetUnitZVector();
	Matrix GetAxialUnitVector();
	Matrix GetShearUnitVector();
	char GetOutOfPlaneDirection();
	~Spring3D();

	//static void AssembleCompleteGlobalMatrixThreads(std::vector<Spring3D> &vecEle, Matrix& complete, std::mutex& mu);
	static void AssembleSpringGlobalMatrixOnComplete(std::vector<Spring3D> &vecEle, Matrix& complete, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp);
	static Matrix AssembleSpringGlobalRestrictedMatrixOnComplete(const int* redSize, const std::map<int, Spring3D*>* vecEle, const AnalysisSpringRecorder* springRecorder);
	//static void AssembleSpringGlobalMatrixOnCompleteDispBased(std::vector<Spring3D> &vecEle, Matrix& k11, Matrix &k12, Matrix &k21, Matrix& k22);
	static Matrix AssembleSpringGlobalMatrixOnReducedSizedMatrix(const int* redSize, const std::map<int, Spring3D*>* vecEle, const AnalysisSpringRecorder* springRecorder);
	static bool CheckMaterialNonlinearityConvergenceDispBased(const StructureManager* structManager, const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecord, bool* breakAnalysis);
	static std::vector<int> GetPlasticDispIndexes(const std::map<int, std::vector<double>>* listOfPlasticDisp, const StructureManager* structManager, const int* DOF);
	static void UpdateSpringLoadStages(const std::map<int, Spring3D*>* vecEle, AnalysisSpringRecorder* springRecord);

private:
	int _ID;
	Node* _n1;
	Node* _n2;
	char _axialDir, _shearDir, _outPlaneDir;
	Matrix _axialUnitVec = Matrix(3, 1);
	Matrix _shearUnitVec = Matrix(3, 1);
	std::vector<std::vector<int>> _globalDOFList;
	std::vector<std::vector<int>> _globalRestDOFList;
	std::vector<SpringMaterialModels*> _listOfMaterials;
};

