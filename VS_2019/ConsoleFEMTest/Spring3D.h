#pragma once
#include "Node.h"
#include "Matrix.h"
#include "SpringMaterialModels.h"
#include "Support.h"
#include <vector>
#include <mutex>

class Spring3D
{
public:
	Spring3D();
	Spring3D(int ID, Node* n1, Node* n2, std::vector<SpringMaterialModels*> listOfMaterials, char axialDir, char shearDir);
	std::string ToString();
	int GetID();
	Node GetNode1();
	Node GetNode2();
	std::vector<SpringMaterialModels*> GetListOfMaterials();
	std::vector<int> GetListOfGlobalMaterialDirections();
	std::vector<std::vector<int>> GetGlobalDOFVector();
	std::vector<std::vector<int>> GetGlobalRestDOFVector();
	Matrix GetLocalStifnessMatrixDispBased(std::vector<double> &listOfPos, std::vector<double> &listOfMinDisp, std::vector<double> &listOfMaxDisp, std::vector<double> &listOfPlasticDisp, std::vector<std::string> &listOfLoadStage, std::vector<std::string> &listOfStage, std::vector<double> &listOfUnlDisp, std::vector<double> &listOfRelDisp);
	Matrix GetTransformationMatrix();
	Matrix GetVectorBetweenNodes();
	Matrix GetGlobalStiffMatrixDispBased(std::vector<double> &listOfPos, std::vector<double> &listOfMinDisp, std::vector<double> &listOfMaxDisp, std::vector<double> &listOfPlasticDisp, std::vector<std::string> &listOfLoadStage, std::vector<std::string> &listOfStage, std::vector<double> &listOfUnlDisp, std::vector<double> &listOfRelDisp);
	std::vector<double> GlobalDOFVectorNotUsed();
	void CalculateGlobalDOFVector(std::vector<Support*> supList);
	Matrix GetElementGlobalDisplacementVector(Matrix completeD);
	Matrix GetUnitZVector();
	Matrix GetAxialUnitVector();
	Matrix GetShearUnitVector();
	char GetOutOfPlaneDirection();
	~Spring3D();

	//static void AssembleCompleteGlobalMatrixThreads(std::vector<Spring3D> &vecEle, Matrix& complete, std::mutex& mu);
	static void AssembleSpringGlobalMatrixOnComplete(std::vector<Spring3D> &vecEle, Matrix& complete, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp);
	static void AssembleSpringGlobalRestrictedMatrixOnCompleteDispBased(std::vector<Spring3D> &vecEle, Matrix& complete, std::vector<Support> &sup, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp);
	//static void AssembleSpringGlobalMatrixOnCompleteDispBased(std::vector<Spring3D> &vecEle, Matrix& k11, Matrix &k12, Matrix &k21, Matrix& k22);
	static void AssembleSpringGlobalMatrixOnCompleteDispBased(std::vector<Spring3D> &vecEle, Matrix& complete, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp);
	static bool CheckMaterialNonlinearityConvergenceDispBased(std::vector<Spring3D> &vecEle, std::vector<std::vector<double>> &oldListOfDisp, std::vector<std::vector<double>> &newListOfDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &oldListOfMaxDisp, std::vector<std::vector<double>> &oldListOfMinDisp, std::vector<std::vector<double>> &oldListOfPlasticDisp, std::vector<std::vector<double>> &newListOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStage, double convLimit, bool &breakAnalysis, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp, std::vector<std::vector<double>> &oldListOfUnlDisp, std::vector<std::vector<double>> &oldListOfRelDisp, std::vector<std::vector<std::string>> &oldSpringStages, std::vector<std::vector<std::string>> &newSpringStages);
	static std::vector<int> GetPlasticDispIndexes(std::vector<std::vector<double>>& listOfPlasticDisp, std::vector<Spring3D> &listOfSpring, std::vector<Support*> listOfSup);
	static void UpdateSpringLoadStages(std::vector<Spring3D> &vecSup, std::vector<std::vector<std::string>> &listOfStages, std::vector<std::vector<std::string>> &stages);

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

