#pragma once
#include "Node.h"
#include "Matrix.h"
#include "ElasticMaterial.h"
#include "OrthotropicElasticMaterial.h"
#include "Support.h"
#include "Mass.h"
#include <vector>
#include <thread>
#include <mutex>

class ShellElement
{
public:
	ShellElement();
	ShellElement(int ID, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, Node* n6, Node* n7, Node* n8, Node* n9, double thick, int layers, std::vector<OrthotropicElasticMaterial> matList);
	double GetWidth();
	double GetHeight();
	double GetTotalThickness();
	double GetLayers();
	double GetLayerThickness();
	int GetID();
	double GetMassRho();
	void SetNodeList(std::vector<Node *> nodeList);
	bool IsDOFRotational(int &DOF);
	std::vector<int> GlobalDOFVector();
	void CalculateGlobalDOFVector(std::vector<Support> &sup);
	void CalculateGlobalMassDOFVector(std::vector<Mass> &massList, std::vector<Support> &supList);
	std::vector<std::vector<int>> GetGlobalDOFVector();
	std::vector<std::vector<int>> GetGlobalMassDOFVector();
	std::vector<std::vector<int>> GetGlobalRestDOFVector();
	Matrix GetJacobian(double &e, double &n, double &c);
	Matrix GetBMatrix(char type, double &e, double &n, double &c, Matrix &Jacobian);
	Matrix GetTotalBMatrix(double &e, double &n, double &c, Matrix &Jacobian);
	Matrix GetDMatrix(OrthotropicElasticMaterial mat);
	Matrix GetGlobalStiffMatrixSeparate();
	Matrix GetGlobalStiffMatrix();
	Matrix GetGlobalStiffMatrixFollowingBook();
	std::vector<Node*> GetNodeList();

	Matrix Get_v1(Matrix &v3);
	Matrix Get_v2(Matrix &v3, Matrix &v1);
	Matrix Get_v3(bool first, double e, double n);
	Matrix GetTransformationMatrix(Matrix &jacob);
	Matrix GetTransformationMatrix3x3Book(Matrix &jacob);
	Matrix GetTransformationMatrix3x3Hyrnyk(Matrix &jacob);
	Matrix GetMassMatrixTheory();
	Matrix GetMassMatrix();

	~ShellElement();

	static double GetRotZStiffTerm(Matrix &m);
	static Matrix GetShearDMatrix(OrthotropicElasticMaterial mat);
	static Matrix GetMembFlexDMatrix(OrthotropicElasticMaterial mat);
	static double ShapeFunctionN(int i, double &e, double &n);
	static double ShapeFunctionNr(int i, double &e, double &n);
	static double ShapeFunctionNDerivative(int i, char var, double &e, double &n);
	//static Matrix AssembleCompleteGlobalMatrix(std::vector<Node> &vecNode, std::vector<ShellElement> &vecEle);
	//static void AssembleCompleteGlobalMatrixThreads(std::vector<ShellElement> &vecEle, Matrix& complete, std::mutex& mu);
	static void AssembleCompleteGlobalMatrixThreads(std::vector<ShellElement> &vecEle, Matrix& complete, std::vector<Support> &sup, std::mutex& mu);
	static void AssembleCompleteRestrictedGlobalMatrixThreads(std::vector<ShellElement> &vecEle, Matrix& complete, std::vector<Support> &sup, std::mutex& mu);
	static void GetNinthNodeDisplacement(Matrix &totalDisplacementMatrix, std::vector<ShellElement> &listElements);
	static Matrix CondensedReducedStiffMatrixForModal(Matrix &m, std::vector<std::vector<int>> &totalMassDOFVec);
	static Matrix GetMassMatrixNonZeroMassOnly(Matrix &m, std::vector<std::vector<int>> &totalMassDOFVec);
	static void AssembleCompleteGlobalMassMatrixThreads(std::vector<ShellElement> &vecEle, Matrix& complete, std::vector<Support> &sup, std::mutex& mu);
	static std::vector<std::vector<int>> CrescentOrderDOFVector(std::vector<std::vector<int>>& oriVec);
	static Matrix ConvertAccFromReducedToTotal(Matrix &acc, std::vector<std::vector<int>> &totalMassDOFVec, int size);
	static Matrix ReducedAccelerationForceMatrix(Matrix &m, std::vector<std::vector<int>>& totalMassDOFVec);
	static bool IsDOFInTheNinthNode(int DOF, std::vector<ShellElement> &vecEle);
	static std::vector<std::vector<int>> GetTotalGlobalMassDOFVector(std::vector<ShellElement> &vecEle);

private:
	int _ID, _layers;
	double _thick;
	double _massRho = 8 * pow(10, -10);
	//double _massRho = 0;
	Matrix _v1 = Matrix(3, 1);
	Matrix _v2 = Matrix(3, 1);
	Matrix _v3 = Matrix(3, 1);
	std::vector<Node *> _nodeList;
	std::vector<OrthotropicElasticMaterial> _matList;
	std::vector<std::vector<int>> _globalMassDOFList;
	std::vector<std::vector<int>> _globalDOFList;
	std::vector<std::vector<int>> _globalRestDOFList;
};
