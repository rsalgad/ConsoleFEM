#include "pch.h"


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

Matrix Displacement::GetTotalDisplacementNotOrganized(Matrix& m, const StructureManager* structManager, const std::vector<std::vector<double*>> listOfPlasticDisp, const int* DOF) {
	Matrix d(structManager->Nodes()->size()*(*DOF), 1);

	for (int i = 0; i < m.GetDimX(); i++) { //copying the terms of the m matrix to the d matrix
		d.GetMatrixDouble()[i][0] = m.GetMatrixDouble()[i][0];
	}

	for (int i = 0; i < listOfPlasticDisp.size(); i++) {

		for (int j = 0; j < listOfPlasticDisp[i].size(); j++) { //for each DOF of the spring
			int nodeID = *structManager->SpringElements()->find(i + 1)->second->GetNode2().GetID();
			int DOFNode = (nodeID - 1) * (*DOF) + j;
			int count = Support::NumberOfDOFBeforeDOF(&DOFNode, structManager->Supports(), DOF);
			d.GetMatrixDouble()[DOFNode - count][0] -= *listOfPlasticDisp[i][j]; //remove the plastic displacement from the displacements
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

std::map<int, Node> Displacement::GetNewNodalCoordinates(const std::map<int, Node*>* oriMap, Matrix& totalDisp) {
	std::map<int, Node> newMap;
	std::map<int, Node*>::const_iterator it = oriMap->begin();
	while (it != oriMap->end()) { //for each node we have in the problem
		int nodeID = it->first;
		Matrix dispNode = GetDisplacementByNodeID(nodeID, totalDisp);
		newMap.insert(std::pair<int, Node> (nodeID, Node(nodeID, *it->second->GetX() + dispNode.GetMatrixDouble()[0][0], *it->second->GetY() + dispNode.GetMatrixDouble()[1][0], *it->second->GetZ() + dispNode.GetMatrixDouble()[2][0], *it->second->GetRx() + dispNode.GetMatrixDouble()[3][0], *it->second->GetRy() + dispNode.GetMatrixDouble()[4][0], *it->second->GetRz() + dispNode.GetMatrixDouble()[5][0])));
		it++;
	}
	return newMap;
}

//This function is responsible for updating all the parameters relevant tot he identification of which region the spring is at in its cyclic material model
void Displacement::UpdatePositionVectorsOfSprings(Matrix* newDisp, const std::map<int, Spring3D*>* vecEle, AnalysisSpringRecorder* springRecorder, const int* DOF) {


	for (int i = 0; i < vecEle->size(); i++){
		std::map<int, Spring3D*>::const_iterator it = vecEle->find(i + 1);
		Spring3D* spring = it->second;
		int springID = spring->GetID();

		std::vector<double*> newList;
		newList.reserve(3);
		double* d1 = new double(newDisp->GetMatrixDouble()[(*spring->GetNode2().GetID() - 1) * (*DOF)][0]);
		double* d2 = new double(newDisp->GetMatrixDouble()[(*spring->GetNode2().GetID() - 1) * (*DOF) + 1][0]);
		double* d3 = new double(newDisp->GetMatrixDouble()[(*spring->GetNode2().GetID() - 1) * (*DOF) + 2][0]);

		newList.emplace_back(d1); //x
		newList.emplace_back(d2);
		newList.emplace_back(d3);

		//make the current 'new' spring pos as the 'old' pos.
		springRecorder->SwapOldNewVectors(springID);

		//update the new spring pos
		(*springRecorder->GetNewDisp())[i][0] = newList[0];
		(*springRecorder->GetNewDisp())[i][1] = newList[1];
		(*springRecorder->GetNewDisp())[i][2] = newList[2];

		std::vector<int> matPos = spring->GetListOfGlobalMaterialDirections();

		for (int j = 0; j < newList.size(); j++) { //for each DOF
			if (*newList[j] > *(*springRecorder->GetMaxDispIter())[i][j]) {
				(*springRecorder->GetNewMaxDisp())[i][j] = newList[j];
			}
			else {
				(*springRecorder->GetNewMaxDisp())[i][j] = (*springRecorder->GetMaxDispIter())[i][j];
			}

			if (*newList[j] < *(*springRecorder->GetMinDispIter())[i][j]) {
				(*springRecorder->GetNewMinDisp())[i][j] = newList[j];
			}
			else {
				(*springRecorder->GetNewMinDisp())[i][j] = (*springRecorder->GetMinDispIter())[i][j];
			}

			//because 'stage' didn't account for the minD and maxD changes, a new stage is calcualted after the minD and maxD values are already accounted for in the previous lines.
			std::string* s1 = new std::string(spring->GetListOfMaterials()[matPos[j]]->GetLoadingStage(*newList[j], *(*springRecorder->GetNewMaxDisp())[i][j], *(*springRecorder->GetNewMinDisp())[i][j],
				*(*springRecorder->GetListStages())[i][j],
				*(*springRecorder->GetUnlDispIter())[i][j],
				*(*springRecorder->GetRelDispIter())[i][j]));

			(*springRecorder->GetNewStages())[i][j] = s1;

			spring->GetListOfMaterials()[matPos[j]]->UpdateUnlAndRelDisps(*(*springRecorder->GetNewStages())[i][j],
				*(*springRecorder->GetListStages())[i][j],
				*newList[j],
				(*springRecorder->GetNewMaxDisp())[i][j],
				(*springRecorder->GetNewMinDisp())[i][j],
				(*springRecorder->GetNewUnlDisp())[i][j],
				(*springRecorder->GetNewRelDisp())[i][j],
				*(*springRecorder->GetUnlDispIter())[i][j],
				*(*springRecorder->GetRelDispIter())[i][j]);

			double* val = new double (spring->GetListOfMaterials()[matPos[j]]->GetPlasticDisplacement(*newList[j],
				*(*springRecorder->GetNewMaxDisp())[i][j],
				*(*springRecorder->GetNewMinDisp())[i][j],
				*(*springRecorder->GetListStages())[i][j],
				*(*springRecorder->GetNewStages())[i][j],
				*(*springRecorder->GetUnlDispIter())[i][j],
				*(*springRecorder->GetRelDispIter())[i][j]));
			(*springRecorder->GetNewPlasticDisp())[i][j] = val;
		}
	}
}

Displacement::~Displacement()
{
}
