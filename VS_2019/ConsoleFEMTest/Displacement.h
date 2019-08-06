#pragma once
#include "Matrix.h"
#include "Support.h"
#include "Node.h"
#include "Spring3D.h"
#include "ShellElement.h"
#include <vector>

class Displacement
{
public:
	Displacement();
	static Matrix GetTotalDisplacementMatrix(Matrix &m, std::vector<Support> &vecSup, std::vector<Node> &vecNode);
	static Matrix GetDisplacementByNodeID(int const &ID, Matrix &totalDisplacementMatrix);
	static Matrix GetRestrainedDispMatrix(std::vector<Support> &vecSup);
	static std::vector<Node> GetNewNodalCoordinates(std::vector<Node> &oriVec, Matrix &totalDisp);
	static Matrix GetTotalDisplacementNotOrganized(Matrix &m, std::vector<Support> &vecSup, std::vector<Node> &vecNode, std::vector<Spring3D> &listOfSpring, std::vector<std::vector<double>> &listOfPlasticDisp);
	static void UpdatePositionVectorsOfSprings(std::vector<std::vector<double>> &oldPos, std::vector<std::vector<double>> &newPos, Matrix &newDisp, std::vector<Spring3D> &vecEle, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &oldListOfMinDisp, std::vector<std::vector<double>> &oldListOfMaxDisp, std::vector<std::vector<double>> &oldListOfPlasticDisp, std::vector<std::vector<double>> &newListOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStages, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp, std::vector<std::vector<double>> &oldListOfUnlDisp, std::vector<std::vector<double>> &oldListOfRelDisp, std::vector<std::vector<double>> &maxDispPerIter, std::vector<std::vector<double>> &minDispPerIter, std::vector<std::vector<double>> &unlDispPerIter, std::vector<std::vector<double>> &relDispPerIter, std::vector<std::vector<std::string>> &oldSpringStages, std::vector<std::vector<std::string>> &newSpringStages);
	static void ZeroOutPositionVectorsOfSprings(std::vector<std::vector<double>> &oldPos, std::vector<std::vector<double>> &newPos, std::vector<Spring3D> &vecEle, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &oldListOfMinDisp, std::vector<std::vector<double>> &oldListOfMaxDisp, std::vector<std::vector<double>> &oldListOfPlasticDisp, std::vector<std::vector<double>> &newListOfPlasticDisp, std::vector<std::vector<std::string>> &listOfSpringLoadingStages, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp, std::vector<std::vector<double>> &oldListOfUnlDisp, std::vector<std::vector<double>> &oldListOfRelDisp, std::vector<std::vector<double>> &maxDispPerIter, std::vector<std::vector<double>> &minDispPerIter, std::vector<std::vector<double>> &unlDispPerIter, std::vector<std::vector<double>> &relDispPerIter, std::vector<std::vector<std::string>> &oldSpringStages, std::vector<std::vector<std::string>> &newSpringStages);
	~Displacement();
};
