#include "pch.h"
#include "Solver.h"
#include "Matrix.h"
#include "ShellElement.h"
#include "Spring3D.h"
#include "Support.h"
#include "Mass.h"
#include "MatrixOperation.h"
#include "StructureManager.h"
#include "PreAnalysisSetUp.h"
#include <iostream>
#include <mutex>
#include <thread>

Solver::Solver()
{
}
/*
Matrix Solver::CompleteStiffnessMatrixWithThreads(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Spring3D> &listOfSprings, int nThreads)
{
	int DOF = 6;
	unsigned int size = listOfNodes.size()*DOF;
	Matrix m(size, size);
	double amount = listOfShells.size() / nThreads;
	int size1, size2;
	if (fmod(listOfShells.size(), nThreads) != 0) {
		size1 = floor(amount);
		size2 = fmod(amount, nThreads) * nThreads;
	}
	else {
		size1 = amount;
		size2 = 0;
	}

	std::vector<std::vector<ShellElement>> shellElemVecs;

	for (int i = 0; i < nThreads; i++) {
		if (i < (nThreads - 1))
		{
			std::vector<ShellElement> list;
			list.assign(listOfShells.begin() + i * size1, listOfShells.begin() + (i + 1)*size1);
			shellElemVecs.push_back(list);
		}
		else {
			std::vector<ShellElement> list;
			list.assign(listOfShells.begin() + i * size1, listOfShells.end());
			shellElemVecs.push_back(list);
		}
	}
	
	std::mutex matrixLock;
	std::vector<std::thread> threadList;
	for (int i = 0; i < (nThreads - 1); i++) {
		threadList.emplace_back(ShellElement::AssembleCompleteGlobalMatrixThreads, std::ref(shellElemVecs[i]), std::ref(m), std::ref(matrixLock));
	}
	ShellElement::AssembleCompleteGlobalMatrixThreads(shellElemVecs.back(), m, matrixLock);

	for (int i = 0; i < threadList.size(); i++) {
		threadList[i].join();
	}
	if (listOfSprings.size() != 0) {
		Spring3D::AssembleSpringGlobalMatrixOnComplete(listOfSprings, m);
	}
	return m;
}
*/

//<summary>Calculates the complete stiffness matrix with threads based on the displacement based theory</summary>

Matrix Solver::ReducedStiffnessMatrix(const StructureManager* structManager,const PreAnalysisSetUp* setUp)
{	
	Matrix m(setUp->ReducedStiffMatrixSize(), setUp->ReducedStiffMatrixSize());

	if (structManager->ShellElements()->size() != 0) {
		std::mutex matrixLock;
		std::vector<std::thread> threadList;
		std::map<int, std::vector<ShellElement*>>::const_iterator it;
		
		for (int i = 0; i < (setUp->AvailableThreads - 1); i++) {
			it = setUp->GetShellThreads()->find(i + 1);
			if (it != setUp->GetShellThreads()->end()) {
				threadList.emplace_back(ShellElement::AssembleCompleteGlobalMatrixThreads, &it->second, std::ref(m), std::ref(matrixLock));
			}
		}

		it = setUp->GetShellThreads()->end();
		ShellElement::AssembleCompleteGlobalMatrixThreads(&it->second, m, matrixLock);

		for (int i = 0; i < threadList.size(); i++) {
			threadList[i].join();
		}
	}


	//Why is this commented out? The springs are not included in the Stiffness Matrix?
	/*
	if (listOfSprings.size() != 0) {
		Spring3D::AssembleSpringGlobalMatrixOnCompleteDispBased(listOfSprings, m);
	}
	*/
	return m;
}

//<summary>Calculates the complete stiffness matrix of the shell elements with threads based on the force based theory</summary>
Matrix Solver::CompleteShellStiffMatrixThreads(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Support> &listOfSupports, int nThreads, std::vector<std::vector<ShellElement>> &shellElemVecs) {
	int DOF = 6;
	int size = listOfNodes.size()*DOF - Support::TotalDOFsRestrained(listOfSupports);
	Matrix m(size, size);

	if (listOfShells.size() != 0) {
		std::mutex matrixLock;
		std::vector<std::thread> threadList;
		for (int i = 0; i < (nThreads - 1); i++) {
			threadList.emplace_back(ShellElement::AssembleCompleteGlobalMatrixThreads, std::ref(shellElemVecs[i]), std::ref(m), std::ref(listOfSupports), std::ref(matrixLock));
		}
		ShellElement::AssembleCompleteGlobalMatrixThreads(shellElemVecs.back(), m, listOfSupports, matrixLock);

		for (int i = 0; i < threadList.size(); i++) {
			threadList[i].join();
		}
	}
	return m;
}

Matrix Solver::CompleteShellMassMatrixThreads(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Support> &listOfSupports, int nThreads, std::vector<std::vector<ShellElement>> &shellElemVecs) {
	int DOF = 6;
	int size = listOfNodes.size()*DOF - Support::TotalDOFsRestrained(listOfSupports);
	Matrix m(size, size);

	if (listOfShells.size() != 0) {
		std::mutex matrixLock;
		std::vector<std::thread> threadList;
		for (int i = 0; i < (nThreads - 1); i++) {
			threadList.emplace_back(ShellElement::AssembleCompleteGlobalMassMatrixThreads, std::ref(shellElemVecs[i]), std::ref(m), std::ref(listOfSupports), std::ref(matrixLock));
		}
		ShellElement::AssembleCompleteGlobalMassMatrixThreads(shellElemVecs.back(), m, listOfSupports, matrixLock);

		for (int i = 0; i < threadList.size(); i++) {
			threadList[i].join();
		}
	}
	return m;
}

//<summary>Calculates the complete stiffness matrix of the shell elements with threads based on the force based theory</summary>
Matrix Solver::CompleteShellRestrictedStiffMatrixThreads(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Support> &listOfSupports, int nThreads, std::vector<std::vector<ShellElement>> &shellElemVecs) {
	int DOF = 6;
	int sizeRow = Support::TotalDOFsRestrained(listOfSupports);
	int sizeCol = listOfNodes.size()*DOF;
	Matrix m(sizeRow, sizeCol);

	if (listOfShells.size() != 0) {
		std::mutex matrixLock;
		std::vector<std::thread> threadList;
		for (int i = 0; i < (nThreads - 1); i++) {
			threadList.emplace_back(ShellElement::AssembleCompleteRestrictedGlobalMatrixThreads, std::ref(shellElemVecs[i]), std::ref(m), std::ref(listOfSupports), std::ref(matrixLock));
		}
		ShellElement::AssembleCompleteRestrictedGlobalMatrixThreads(shellElemVecs.back(), m, listOfSupports, matrixLock);

		for (int i = 0; i < threadList.size(); i++) {
			threadList[i].join();
		}
	}
	return m;
}

void Solver::CompleteSpringStiffMatrixThreadsDispBasedAfterShells(std::vector<Spring3D> &listOfSprings, Matrix &m, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp) {
	if (listOfSprings.size() != 0) {
		Spring3D::AssembleSpringGlobalMatrixOnCompleteDispBased(listOfSprings, m, listOfDisp, listOfMinDisp, listOfMaxDisp, listOfPlasticDisp, listOfLoadStage, listOfStage, listOfUnlDisp, listOfRelDisp);
	}
}

void Solver::CompleteSpringRestrictedStiffMatrixThreadsDispBasedAfterShells(std::vector<Spring3D> &listOfSprings, Matrix &m, std::vector<Support> &sup, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp) {
	if (listOfSprings.size() != 0) {
		Spring3D::AssembleSpringGlobalRestrictedMatrixOnCompleteDispBased(listOfSprings, m, sup, listOfDisp, listOfMinDisp, listOfMaxDisp, listOfPlasticDisp, listOfLoadStage, listOfStage, listOfUnlDisp, listOfRelDisp);
	}
}

void Solver::DisplacementLoadStiffness(Matrix& stiff, std::vector<Support> &listOfSup, double &biggest) {
	double highStiff = biggest * pow(10, 6);
	std::vector<int> indexes = Support::GetDisplacementLoadIndexes(listOfSup); //the DOF indexes of the supports with displacement loads
	
	for (int i = 0; i < indexes.size(); i++) {
		stiff.GetMatrixDouble()[indexes[i]][indexes[i]] += highStiff;
	}
}

void Solver::DisplacementLoadForce(Matrix& force, std::vector<Support> &listOfSup, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<Spring3D> &listOfSpring, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, double &biggest, double loadFraction, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp) {
	double highStiff = biggest * pow(10, 6);
	std::vector<int> indexes = Support::GetDisplacementLoadIndexes(listOfSup);
	
	int count = 0;

	for (int i = 0; i < listOfSup.size(); i++) {
		for (int j = 0; j < listOfSup[i].GetSupportVector().size(); j++) {
			if (listOfSup[i].GetSupportVector()[j][1] != 0){ //no need to do this if the support is just a boundary condition
				double dispLoad = listOfSup[i].GetSupportVector()[j][1];
				force.GetMatrixDouble()[indexes[count]][0] += highStiff * dispLoad * loadFraction;
				count++;
			}
		}
	}

	std::vector<int> plasticIndexes = Spring3D::GetPlasticDispIndexes(listOfPlasticDisp, listOfSpring, listOfSup);
	int count2 = 0;
	for (int i = 0; i < listOfPlasticDisp.size(); i++) {
		std::vector<int> matPos = listOfSpring[i].GetListOfGlobalMaterialDirections();
		for (int j = 0; j < listOfPlasticDisp[i].size(); j++) {
			if (listOfPlasticDisp[i][j] != 0) { //no need to do this if no plastic disp
				double dispLoad = listOfPlasticDisp[i][j];
				double stiff = listOfSpring[i].GetListOfMaterials()[matPos[j]]->GetSecantStiffnessFromDisplacement(listOfDisp[i][j], listOfPlasticDisp[i][j], listOfMaxDisp[i][j], listOfMinDisp[i][j], listOfLoadStage[i][j], listOfStage[i][j], listOfUnlDisp[i][j], listOfRelDisp[i][j]);
				//double force22 = listOfSpring[i].GetListOfMaterials()[matPos[j]]->GetForceFromDisplacement(dispLoad, listOfMaxDisp[i][j], listOfMinDisp[i][j]);
				double force2 = stiff * dispLoad;
				force.GetMatrixDouble()[plasticIndexes[count2]][0] += force2; //no *loadfraction since I want to apply the entire plastic disp at once
				count2++;
			}
		}
	}
}

void Solver::PlasticDisplacementLoadForce(Matrix& force, std::vector<Support> &listOfSup, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<Spring3D> &listOfSpring, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp) {
	std::vector<int> plasticIndexes = Spring3D::GetPlasticDispIndexes(listOfPlasticDisp, listOfSpring, listOfSup);
	int count2 = 0;
	for (int i = 0; i < listOfPlasticDisp.size(); i++) {
		std::vector<int> matPos = listOfSpring[i].GetListOfGlobalMaterialDirections();
		for (int j = 0; j < listOfPlasticDisp[i].size(); j++) {
			if (listOfPlasticDisp[i][j] != 0) { //no need to do this if no plastic disp
				double dispLoad = listOfPlasticDisp[i][j];
				double stiff = listOfSpring[i].GetListOfMaterials()[matPos[j]]->GetSecantStiffnessFromDisplacement(listOfDisp[i][j], listOfPlasticDisp[i][j], listOfMaxDisp[i][j], listOfMinDisp[i][j], listOfLoadStage[i][j], listOfStage[i][j], listOfUnlDisp[i][j], listOfRelDisp[i][j]);
				//double force2 = listOfSpring[i].GetListOfMaterials()[matPos[j]]->GetForceFromDisplacement(dispLoad, listOfMaxDisp[i][j], listOfMinDisp[i][j]);
				double force2 = stiff * dispLoad;
				force.GetMatrixDouble()[plasticIndexes[count2]][0] += force2; //no *loadfraction since I want to apply the entire plastic disp at once
				count2++;
			}
		}
	}
}

void Solver::CalculateNaturalFrequenciesAndModeShapes(Matrix &stiffMatrix, Matrix &massMatrix, std::vector<double> &natFreq, Matrix &modeShapes, std::vector<std::vector<int>> &totalMassDOFVec) {
	
	Matrix reducedStiff = ShellElement::CondensedReducedStiffMatrixForModal(stiffMatrix, totalMassDOFVec);
	
	Matrix reducedMass = ShellElement::GetMassMatrixNonZeroMassOnly(massMatrix, totalMassDOFVec);

	Matrix sqrtMass = MatrixOperation::Sqrt(reducedMass);
	Matrix invSqrtMass = MatrixOperation::GetInverse(sqrtMass);

	Matrix mult2 = invSqrtMass * reducedStiff * invSqrtMass;

	Eigen::MatrixXd eigenMatrix;
	MatrixOperation::ConvertToEigenMatrix(mult2, eigenMatrix);

	Eigen::EigenSolver<Eigen::MatrixXd> s(eigenMatrix); // the instance s(A) includes the eigensystem
	Eigen::VectorXcd eigenMatrix1 = s.eigenvalues();
	Eigen::MatrixXcd eigenMatrix2 = s.eigenvectors();

	natFreq.reserve(s.eigenvalues().size());

	for (int i = 0; i < s.eigenvalues().size(); i++) {
		natFreq.emplace_back(sqrt(eigenMatrix1(i).real()));
	}

	modeShapes.SetDimensions(s.eigenvalues().size(), s.eigenvalues().size());
	for (int i = 0; i < s.eigenvalues().size(); i++) { //for each mode (column)
		for (int j = 0; j < s.eigenvalues().size(); j++) { //for each DOF (row)
			modeShapes.GetMatrixDouble()[j][i] = eigenMatrix2(j, i).real();
		}
	}
}

std::vector<double> Solver::RayleighDampingConstants(int mode1, double damp1, int mode2, double damp2, std::vector<double> &natFreq) {

	double alpha1 = (2 * damp2 * natFreq[mode2 - 1] - 2 * damp1 * natFreq[mode1 - 1]) / (pow(natFreq[mode2 - 1], 2) - pow(natFreq[mode1 - 1], 2));
	double alpha0 = 2 * damp1 * natFreq[mode1 - 1] - alpha1 * pow(natFreq[mode1 - 1], 2);

	std::vector<double> vec;
	vec.reserve(2);
	vec.emplace_back(alpha0);
	vec.emplace_back(alpha1);

	return vec;
}

Matrix Solver::RayleighDampingMatrix(Matrix &m, Matrix &k, std::vector<double>& constants) {

	Matrix c = m * constants[0] + k * constants[1];

	return c;
}

Matrix Solver::GetTotalModalMatrix(Matrix &m, std::vector<Support> &vecSup, std::vector<Node> &vecNode, std::vector<ShellElement> &vecShell, std::vector<Mass> &vecMass) {
	//This function returns the full displacement vector from the reduced version
	int DOF = 6; // for shell elements
	int size = vecNode.size() * DOF;
	int count = 0;
	Matrix disp(size, 1);
	Matrix ans(size, 0); //0 is intentional

	std::vector<int> vecPos;
	std::vector<double> vecDisp;

	for (int i = 0; i < vecSup.size(); i++) //for each support
	{
		for (int j = 0; j < vecSup[i].GetSupportVector().size(); j++) { //this will store all the global positions and displacement values in vectors
			if (vecSup[i].GetSupportVector()[j][1] == 0) { //only do this if the support is a boundary condition
				int nodeID = vecSup[i].GetNode();
				int dir = vecSup[i].GetSupportVector()[j][0];
				vecPos.emplace_back((nodeID - 1)*DOF + (dir - 1));
				vecDisp.emplace_back(0);
			}
		}
	}

	int index = 0;
	int i = 0;
	for (int k = 0; k < m.GetDimY(); k++) {
		while (i < size) {
			for (int j = 0; j < 3; j++) {
				bool notSup = !(std::find(vecPos.begin(), vecPos.end(), i + j) != vecPos.end());
				if (notSup && (HasDOFMass(vecShell, vecSup, i + j) || Mass::HasDOFAppliedMass(vecMass, i + j))) { //if it has mass
					disp.GetMatrixDouble()[i + j][0] = m.GetMatrixDouble()[index][k];
					index++;
				}
			}
			i += 6;
		}
		ShellElement::GetNinthNodeDisplacement(disp, vecShell);
		ans = MatrixOperation::AddMatrixRight(ans, disp);
		i = 0;
		index = 0;
	}
	return ans;
}

bool Solver::HasDOFMass(std::vector<ShellElement> &shellVec, std::vector<Support> &supVec, int DOF) {
	int nDOF = 6;
	int count = Support::NumberOfDOFBeforeDOF(DOF, supVec);
	int reducedDOF = DOF - count;
	if (!Support::IsDOFConstrained(DOF, supVec)) {
		for (int i = 0; i < shellVec.size(); i++) { //for each shell
			std::vector<std::vector<int>> vec = shellVec[i].GetGlobalMassDOFVector();
			for (int j = 0; j < vec.size(); j++) {
				if (vec[j][0] == reducedDOF) {
					return false;
				}
			}
		}
	}
	else {
		return false;
	}
	return true;
}

Solver::~Solver()
{
}
