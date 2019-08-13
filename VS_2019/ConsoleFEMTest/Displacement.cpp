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

Matrix Displacement::GetTotalDisplacementNotOrganized(const Matrix* m, const StructureManager* structManager, const std::map<int, std::vector<double>>* listOfPlasticDisp, const int* DOF) {
	Matrix d(structManager->Nodes()->size()*(*DOF), 1);

	for (int i = 0; i < m->GetDimX(); i++) { //copying the terms of the m matrix to the d matrix
		d.GetMatrixDouble()[i][0] = m->GetMatrixDouble()[i][0];
	}

	std::map<int, std::vector<double>>::const_iterator it = listOfPlasticDisp->begin();

	while (it != listOfPlasticDisp->end()) {
		int springID = it->first;
		for (int j = 0; j < it->second.size(); j++) { //for each DOF of the spring
			int nodeID = structManager->SpringElements()->find(springID)->second->GetNode2().GetID();
			int DOFNode = (nodeID - 1) * (*DOF) + j;
			int count = Support::NumberOfDOFBeforeDOF(&DOFNode, structManager->Supports());
			d.GetMatrixDouble()[DOFNode - count][0] -= it->second[j]; //remove the plastic displacement from the displacements
		}
		it++;
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

std::map<int, Node*> Displacement::GetNewNodalCoordinates(const std::map<int, Node*>* oriMap, Matrix& totalDisp) {
	std::map<int, Node*> newMap;
	std::map<int, Node*>::const_iterator it = oriMap->begin();
	while (it != oriMap->end()) { //for each node we have in the problem
		int nodeID = it->first;
		Matrix dispNode = GetDisplacementByNodeID(nodeID, totalDisp);
		newMap.insert(std::pair<int, Node*> (nodeID, new Node(nodeID, it->second->GetX() + dispNode.GetMatrixDouble()[0][0], it->second->GetY() + dispNode.GetMatrixDouble()[1][0], it->second->GetZ() + dispNode.GetMatrixDouble()[2][0], it->second->GetRx() + dispNode.GetMatrixDouble()[3][0], it->second->GetRy() + dispNode.GetMatrixDouble()[4][0], it->second->GetRz() + dispNode.GetMatrixDouble()[5][0])));
	}
	return newMap;
}

//This function is responsible for updating all the parameters relevant tot he identification of which region the spring is at in its cyclic material model
void Displacement::UpdatePositionVectorsOfSprings(Matrix* newDisp, const std::map<int, Spring3D*>* vecEle, const AnalysisSpringRecorder* springRecorder, const int* DOF) {

	std::map<int, Spring3D*>::const_iterator it = vecEle->begin();
	while (it != vecEle->end()) {
		Spring3D* spring = it->second;


		springRecorder->GetDisplacementMap().find("disp")->second.find("new")->second.find(spring->GetID());
		std::vector<double> newList;
		newList.reserve(3);
		newList.emplace_back(newDisp->GetMatrixDouble()[(spring->GetNode2().GetID() - 1) * (*DOF)][0]); //x
		newList.emplace_back(newDisp->GetMatrixDouble()[(spring->GetNode2().GetID() - 1) * (*DOF) + 1][0]); //y
		newList.emplace_back(newDisp->GetMatrixDouble()[(spring->GetNode2().GetID() - 1) * (*DOF) + 2][0]); //z

		//make the current 'new' spring pos as the 'old' pos.
		springRecorder->GetDisplacementMap().find("disp")->second.find("old")->second.find(spring->GetID())->second = springRecorder->GetDisplacementMap().find("disp")->second.find("new")->second.find(spring->GetID())->second;
		
		//update the new spring pos
		springRecorder->GetDisplacementMap().find("disp")->second.find("new")->second.find(spring->GetID())->second = newList;

		//updates the old list of min disp with the actual values before changing them
		springRecorder->GetDisplacementMap().find("minDisp")->second.find("old")->second.find(spring->GetID())->second = springRecorder->GetDisplacementMap().find("minDisp")->second.find("new")->second.find(spring->GetID())->second;
		springRecorder->GetDisplacementMap().find("maxDisp")->second.find("old")->second.find(spring->GetID())->second = springRecorder->GetDisplacementMap().find("maxDisp")->second.find("new")->second.find(spring->GetID())->second;
		springRecorder->GetDisplacementMap().find("plasticDisp")->second.find("old")->second.find(spring->GetID())->second = springRecorder->GetDisplacementMap().find("plasticDisp")->second.find("new")->second.find(spring->GetID())->second;
		springRecorder->GetDisplacementMap().find("unlDisp")->second.find("old")->second.find(spring->GetID())->second = springRecorder->GetDisplacementMap().find("unlDisp")->second.find("new")->second.find(spring->GetID())->second;
		springRecorder->GetDisplacementMap().find("relDisp")->second.find("old")->second.find(spring->GetID())->second = springRecorder->GetDisplacementMap().find("relDisp")->second.find("new")->second.find(spring->GetID())->second;
		springRecorder->GetStagesMap().find("old")->second.find(spring->GetID())->second = springRecorder->GetStagesMap().find("old")->second.find(spring->GetID())->second;

		std::vector<int> matPos = spring->GetListOfGlobalMaterialDirections();

		for (int j = 0; j < newList.size(); j++) { //for each DOF
			if (newList[j] > springRecorder->GetDisplacementPerIterMap().find("maxDispIter")->second.find(spring->GetID())->second[j]) {
				springRecorder->GetDisplacementMap().find("maxDisp")->second.find("new")->second.find(spring->GetID())->second[j] = newList[j];
			}
			else {
				springRecorder->GetDisplacementMap().find("maxDisp")->second.find("new")->second.find(spring->GetID())->second[j] = springRecorder->GetDisplacementPerIterMap().find("maxDispIter")->second.find(spring->GetID())->second[j];
			}

			if (newList[j] < springRecorder->GetDisplacementPerIterMap().find("minDispIter")->second.find(spring->GetID())->second[j]) {
				springRecorder->GetDisplacementMap().find("minDisp")->second.find("new")->second.find(spring->GetID())->second[j] = newList[j];
			}
			else {
				springRecorder->GetDisplacementMap().find("minDisp")->second.find("new")->second.find(spring->GetID())->second[j] = springRecorder->GetDisplacementPerIterMap().find("minDispIter")->second.find(spring->GetID())->second[j];
			}

			//because 'stage' didn't account for the minD and maxD changes, a new stage is calcualted after the minD and maxD values are already accounted for in the previous lines.
			springRecorder->GetStagesMap().find("new")->second.find(spring->GetID())->second[j] 
				= spring->GetListOfMaterials()[matPos[j]]->GetLoadingStage(newList[j], springRecorder->GetDisplacementMap().find("maxDisp")->second.find("new")->second.find(spring->GetID())->second[j],
					springRecorder->GetDisplacementMap().find("minDisp")->second.find("new")->second.find(spring->GetID())->second[j],
					springRecorder->GetStagesMap().find("list")->second.find(spring->GetID())->second[j],
					springRecorder->GetDisplacementPerIterMap.find("unlDispPerIter")->second.find(spring->GetID())->second[j],
					springRecorder->GetDisplacementPerIterMap.find("relDispPerIter")->second.find(spring->GetID())->second[j]);

			spring->GetListOfMaterials()[matPos[j]]->UpdateUnlAndRelDisps(springRecorder->GetStagesMap().find("new")->second.find(spring->GetID())->second[j],
				springRecorder->GetStagesMap().find("list")->second.find(spring->GetID())->second[j],
				newList[j], springRecorder->GetDisplacementMap().find("maxDisp")->second.find("new")->second.find(spring->GetID())->second[j],
				springRecorder->GetDisplacementMap().find("minDisp")->second.find("new")->second.find(spring->GetID())->second[j],
				springRecorder->GetDisplacementMap().find("unlDisp")->second.find("new")->second.find(spring->GetID())->second[j],
				springRecorder->GetDisplacementMap().find("relDisp")->second.find("new")->second.find(spring->GetID())->second[j],
				springRecorder->GetDisplacementPerIterMap.find("unlDispPerIter")->second.find(spring->GetID())->second[j],
				springRecorder->GetDisplacementPerIterMap.find("relDispPerIter")->second.find(spring->GetID())->second[j]);

			springRecorder->GetDisplacementMap().find("plasticDisp")->second.find("new")->second.find(spring->GetID())->second[j] = spring->GetListOfMaterials()[matPos[j]]->GetPlasticDisplacement(newList[j],
				springRecorder->GetDisplacementMap().find("maxDisp")->second.find("new")->second.find(spring->GetID())->second[j],
				springRecorder->GetDisplacementMap().find("minDisp")->second.find("new")->second.find(spring->GetID())->second[j],
				springRecorder->GetStagesMap().find("list")->second.find(spring->GetID())->second[j],
				springRecorder->GetStagesMap().find("new")->second.find(spring->GetID())->second[j],
				springRecorder->GetDisplacementPerIterMap.find("unlDispPerIter")->second.find(spring->GetID())->second[j],
				springRecorder->GetDisplacementPerIterMap.find("relDispPerIter")->second.find(spring->GetID())->second[j]);
		}
		it++;
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
