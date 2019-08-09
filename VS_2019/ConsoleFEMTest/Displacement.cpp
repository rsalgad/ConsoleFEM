#include "pch.h"
#include "Displacement.h"
#include "Spring3D.h"
#include "ShellElement.h"
#include <iostream>
#include <map>


Displacement::Displacement()
{
}

Matrix Displacement::GetRestrainedDispMatrix(std::vector<Support> &vecSup) {
	int count = 0;
	
	for (int i = 0; i < vecSup.size(); i++) {
		for (int j = 0; j < vecSup[i].GetSupportVector().size(); j++) {
			count++;
		}
	}

	Matrix m(count, 1);

	int count1 = 0;
	for (int i = 0; i < vecSup.size(); i++) {
		for (int j = 0; j < vecSup[i].GetSupportVector().size(); j++) {

			m.GetMatrixDouble()[count1][0] = vecSup[i].GetSupportVector()[j][1];
			count1++;
		}
	}

	return m;
}


Matrix Displacement::GetTotalDisplacementMatrix(Matrix &m, const StructureManager* structManager, const PreAnalysisSetUp* setUp) {
	//This function returns the full displacement vector from the reduced version
	int count = 0;
	double** disp = Matrix::CreateMatrixDouble(*setUp->StiffMatrixSize(), 1);

	std::vector<int> vecPos;
	std::vector<double> vecDisp;

	std::map<int, Support*>::const_iterator it = structManager->Supports()->begin();
	
	//for each support
	while (it != structManager->Supports()->end()) {
		for (int j = 0; j < it->second->GetSupportVector().size(); j++) { //this will store all the global positions and displacement values in vectors
			if (it->second->GetSupportVector()[j][1] == 0) { //only do this if the support is a boundary condition
				int nodeID = it->second->GetNode();
				int dir = it->second->GetSupportVector()[j][0];
				double disp = it->second->GetSupportVector()[j][1];
				vecPos.emplace_back((nodeID - 1) * (*setUp->DOF()) + (dir - 1));
				vecDisp.emplace_back(disp);
			}
		}
		it++;
	}

	int index = 0;
	for (int i = 0; i < *setUp->StiffMatrixSize(); i++) { //for each term of the global D matrix
		if (std::find(vecPos.begin(), vecPos.end(), i) != vecPos.end()) {//if i is within the vec of Position of supports
			disp[i][0] = vecDisp[index];
			index++;
		}
		else {
			disp[i][0] = m.GetMatrixDouble()[i-index][0];
		}
	}
	return Matrix(disp, *setUp->StiffMatrixSize(), 1);
}

Matrix Displacement::GetTotalDisplacementNotOrganized(Matrix &m, std::vector<Support*> vecSup, std::vector<Node> &vecNode, std::vector<Spring3D> &listOfSpring, std::vector<std::vector<double>> &listOfPlasticDisp) {
	int DOF = 6;
	Matrix d(vecNode.size()*DOF, 1);

	for (int i = 0; i < m.GetDimX(); i++) { //copying the terms of the m matrix to the d matrix
		d.GetMatrixDouble()[i][0] = m.GetMatrixDouble()[i][0];
	}

	for (int i = 0; i < listOfPlasticDisp.size(); i++) { //for each spring
		for (int j = 0; j < listOfPlasticDisp[i].size(); j++) { //for each DOF of the spring
			int nodeID = listOfSpring[i].GetNode2().GetID();
			int DOFNode = (nodeID - 1) * DOF + j;
			int count = Support::NumberOfDOFBeforeDOF(DOFNode, vecSup);
			d.GetMatrixDouble()[DOFNode - count][0] -= listOfPlasticDisp[i][j]; //remove the plastic displacement from the displacements
		}
	}

	/*
	for (int i = 0; i < d.GetDimX() - m.GetDimX(); i++) {
		d.GetMatrixDouble()[m.GetDimX() + i][0] = 0;
	}
	*/
	return d;
}

Matrix Displacement::GetDisplacementByNodeID(int const &ID, Matrix &totalDisplacementMatrix) {
	int DOF = 6; //for Shell Elements
	double** disp = Matrix::CreateMatrixDouble(DOF);

	for (int i = 0; i < DOF; i++) {
		disp[i][0] = totalDisplacementMatrix.GetMatrixDouble()[(ID - 1)*DOF + i][0];
	}
	return Matrix(disp, DOF, 1);
}

std::vector<Node> Displacement::GetNewNodalCoordinates(std::vector<Node> &oriVec, Matrix &totalDisp) {
	std::vector<Node> newVec;
	newVec.reserve(oriVec.size());
	for (int i = 0; i < oriVec.size(); i++) { //for each node we have in the problem
		int nodeID = i + 1;
		Matrix dispNode = GetDisplacementByNodeID(nodeID, totalDisp);
		newVec.emplace_back(nodeID, oriVec[i].GetX() + dispNode.GetMatrixDouble()[0][0], oriVec[i].GetY() + dispNode.GetMatrixDouble()[1][0], oriVec[i].GetZ() + dispNode.GetMatrixDouble()[2][0], oriVec[i].GetRx() + dispNode.GetMatrixDouble()[3][0], oriVec[i].GetRy() + dispNode.GetMatrixDouble()[4][0], oriVec[i].GetRz() + dispNode.GetMatrixDouble()[5][0]);
	}
	return newVec;
}

std::map<int, Node*> Displacement::GetNewNodalCoordinates(const std::map<int, Node*>* oriMap, Matrix& totalDisp) {
	std::map<int, Node*> newMap;
	std::map<int, Node*>::const_iterator it = oriMap.begin();
	while (it != oriMap.end()) { //for each node we have in the problem
		int nodeID = it->first;
		Matrix dispNode = GetDisplacementByNodeID(nodeID, totalDisp);
		newMap.insert(std::pair<int, Node*> (nodeID, new Node(nodeID, it->second->GetX() + dispNode.GetMatrixDouble()[0][0], it->second->GetY() + dispNode.GetMatrixDouble()[1][0], it->second->GetZ() + dispNode.GetMatrixDouble()[2][0], it->second->GetRx() + dispNode.GetMatrixDouble()[3][0], it->second->GetRy() + dispNode.GetMatrixDouble()[4][0], it->second->GetRz() + dispNode.GetMatrixDouble()[5][0])));
	}
	return newMap;
}

//This function is responsible for updating all the parameters relevant tot he identification of which region the spring is at in its cyclic material model
void Displacement::UpdatePositionVectorsOfSprings(std::vector<std::vector<double>> &oldPos, std::vector<std::vector<double>> &newPos, Matrix &newDisp, std::vector<Spring3D> &vecEle, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &oldListOfMinDisp, std::vector<std::vector<double>> &oldListOfMaxDisp, std::vector<std::vector<double>> &oldListOfPlasticDisp, std::vector<std::vector<double>> &newListOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStages, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp, std::vector<std::vector<double>> &oldListOfUnlDisp, std::vector<std::vector<double>> &oldListOfRelDisp, std::vector<std::vector<double>> &maxDispPerIter, std::vector<std::vector<double>> &minDispPerIter, std::vector<std::vector<double>> &unlDispPerIter, std::vector<std::vector<double>> &relDispPerIter, std::vector<std::vector<std::string>> &oldSpringStages, std::vector<std::vector<std::string>> &newSpringStages) {

	int DOF = 6;

	for (int i = 0; i < vecEle.size(); i++) { //for each spring element
		Spring3D* spring = &vecEle[i];
		std::vector<double> pos;
		pos.reserve(3);
		pos.emplace_back(newDisp.GetMatrixDouble()[(spring->GetNode2().GetID() - 1)*DOF][0]); //x
		pos.emplace_back(newDisp.GetMatrixDouble()[(spring->GetNode2().GetID() - 1)*DOF + 1][0]); //y
		pos.emplace_back(newDisp.GetMatrixDouble()[(spring->GetNode2().GetID() - 1)*DOF + 2][0]); //z

		std::vector<double> newList;
		newList.reserve(3);
		newList.emplace_back(pos[0]);
		newList.emplace_back(pos[1]);
		newList.emplace_back(pos[2]);

		oldPos[i] = newPos[i]; //make the current 'new' spring pos as the 'old' pos.
		newPos[i] = newList; //update the new spring pos

		oldListOfMinDisp[i] = listOfMinDisp[i]; //updates the old list of min disp with the actual values before changing them
		oldListOfMaxDisp[i] = listOfMaxDisp[i];
		oldListOfPlasticDisp[i] = newListOfPlasticDisp[i];

		oldListOfUnlDisp[i] = listOfUnlDisp[i];
		oldListOfRelDisp[i] = listOfRelDisp[i];

		oldSpringStages[i] = newSpringStages[i];

		std::vector<int> matPos = spring->GetListOfGlobalMaterialDirections();

		for (int j = 0; j < newList.size(); j++) {
			if (newList[j] > maxDispPerIter[i][j]) {
				listOfMaxDisp[i][j] = newList[j];
			}
			else {
				listOfMaxDisp[i][j] = maxDispPerIter[i][j];
			}
			
			if (newList[j] < minDispPerIter[i][j]) {
				listOfMinDisp[i][j] = newList[j];
			}
			else {
				listOfMinDisp[i][j] = minDispPerIter[i][j];
			}

			//becauase 'stage' didn't account for the minD and maxD changes, a new stage is calcualted after the minD and maxD values are already accounted for in the previous lines.
			newSpringStages[i][j] = spring->GetListOfMaterials()[matPos[j]]->GetLoadingStage(newList[j], listOfMaxDisp[i][j], listOfMinDisp[i][j], listOfLoadStages[i][j], unlDispPerIter[i][j], relDispPerIter[i][j]);

			spring->GetListOfMaterials()[matPos[j]]->UpdateUnlAndRelDisps(newSpringStages[i][j], listOfLoadStages[i][j], newList[j], listOfMaxDisp[i][j], listOfMinDisp[i][j], listOfUnlDisp[i][j], listOfRelDisp[i][j], unlDispPerIter[i][j], relDispPerIter[i][j]);

			newListOfPlasticDisp[i][j] = spring->GetListOfMaterials()[matPos[j]]->GetPlasticDisplacement(newList[j], listOfMaxDisp[i][j], listOfMinDisp[i][j], listOfLoadStages[i][j], newSpringStages[i][j], unlDispPerIter[i][j], relDispPerIter[i][j]);
		}
	}
}

void Displacement::ZeroOutPositionVectorsOfSprings(std::vector<std::vector<double>> &oldPos, std::vector<std::vector<double>> &newPos, std::vector<Spring3D> &vecEle, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &oldListOfMinDisp, std::vector<std::vector<double>> &oldListOfMaxDisp, std::vector<std::vector<double>> &oldListOfPlasticDisp, std::vector<std::vector<double>> &newListOfPlasticDisp, std::vector<std::vector<std::string>> &listOfSpringLoadingStages, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp, std::vector<std::vector<double>> &oldListOfUnlDisp, std::vector<std::vector<double>> &oldListOfRelDisp, std::vector<std::vector<double>> &maxDispPerIter, std::vector<std::vector<double>> &minDispPerIter, std::vector<std::vector<double>> &unlDispPerIter, std::vector<std::vector<double>> &relDispPerIter, std::vector<std::vector<std::string>> &oldSpringStages, std::vector<std::vector<std::string>> &newSpringStages) {

	std::vector<double> vec;
	vec.reserve(3);
	vec.emplace_back(0);
	vec.emplace_back(0);
	vec.emplace_back(0);

	std::vector<std::string> strVec;
	strVec.reserve(3);
	strVec.emplace_back("initial");
	strVec.emplace_back("initial");
	strVec.emplace_back("initial");

	oldPos.reserve(vecEle.size());
	newPos.reserve(vecEle.size());
	listOfMinDisp.reserve(vecEle.size());
	listOfMaxDisp.reserve(vecEle.size());
	oldListOfMaxDisp.reserve(vecEle.size());
	oldListOfMinDisp.reserve(vecEle.size());
	oldListOfPlasticDisp.reserve(vecEle.size());
	newListOfPlasticDisp.reserve(vecEle.size());
	listOfSpringLoadingStages.reserve(vecEle.size());
	listOfUnlDisp.reserve(vecEle.size());
	listOfRelDisp.reserve(vecEle.size());
	oldListOfUnlDisp.reserve(vecEle.size());
	oldListOfRelDisp.reserve(vecEle.size());
	maxDispPerIter.reserve(vecEle.size());
	minDispPerIter.reserve(vecEle.size());
	unlDispPerIter.reserve(vecEle.size());
	relDispPerIter.reserve(vecEle.size());

	for (int i = 0; i < vecEle.size(); i++) { //for each element
		oldPos.emplace_back(vec);
		newPos.emplace_back(vec);
		listOfMinDisp.emplace_back(vec);
		listOfMaxDisp.emplace_back(vec);
		oldListOfMaxDisp.emplace_back(vec);
		oldListOfMinDisp.emplace_back(vec);
		oldListOfPlasticDisp.emplace_back(vec);
		newListOfPlasticDisp.emplace_back(vec);
		listOfSpringLoadingStages.emplace_back(strVec);
		oldSpringStages.emplace_back(strVec);
		newSpringStages.emplace_back(strVec);
		listOfUnlDisp.emplace_back(vec);
		listOfRelDisp.emplace_back(vec);
		oldListOfUnlDisp.emplace_back(vec);
		oldListOfRelDisp.emplace_back(vec);
		maxDispPerIter.emplace_back(vec);
		minDispPerIter.emplace_back(vec);
		unlDispPerIter.emplace_back(vec);
		relDispPerIter.emplace_back(vec);
	}
}

Displacement::~Displacement()
{
}
