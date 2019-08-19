#pragma once
//#include "pch.h"
#include "Matrix.h"
#include "StructureManager.h"
#include "PreAnalysisSetUp.h"
#include "Support.h"
#include "Node.h"
#include "Spring3D.h"
#include "AnalysisSpringRecorder.h"

#include <map>
#include <vector>

class Displacement
{
public:
	Displacement();
	static Matrix GetTotalDisplacementMatrix(Matrix &m, const StructureManager* structManager, const PreAnalysisSetUp* setUp);
	static Matrix GetDisplacementByNodeID(int const &ID, Matrix &totalDisplacementMatrix);
	static Matrix GetRestrainedDispMatrix(std::vector<Support> &vecSup);
	static std::map<int, Node> GetNewNodalCoordinates(const std::map<int, Node*>* oriMap, Matrix& totalDisp);
	static Matrix GetTotalDisplacementNotOrganized(Matrix& m, const StructureManager* structManager, const std::vector<std::vector<double*>> listOfPlasticDisp, const int* DOF);
	static void UpdatePositionVectorsOfSprings(Matrix* newDisp, const std::map<int, Spring3D*>* vecEle, AnalysisSpringRecorder* springRecorder, const int* DOF);
	~Displacement();
};

