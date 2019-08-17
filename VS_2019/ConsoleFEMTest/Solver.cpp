#include "pch.h"
#include "Solver.h"
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

Matrix Solver::ReducedStiffnessMatrix(const Matrix* shellStiff, const StructureManager* structManager,const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecorder, double* highStiff)
{	
	Matrix springStiff = ReducedSpringStiffMatrix(setUp->ReducedStiffMatrixSize(), structManager->SpringElements(), springRecorder);
	Matrix redStiff = *shellStiff + springStiff;

	if (setUp->DispLoadDOFs() != 0) { //if there are displacement loads
		DisplacementLoadStiffness(redStiff, setUp->DispLoadDOFs(), highStiff); //adds the big stiffness term to account for displacement loads, if any
	}

	return *shellStiff + springStiff;
}

Matrix Solver::ReducedDynamicStiffnessMatrix(const PreAnalysisSetUp* setUp, const Matrix* mass, const Matrix* damp, Matrix* add, Matrix* mult2)
{
	DynamicAnalysis* analysis = static_cast<DynamicAnalysis*>(setUp->Analysis());
	double const1 = (1 - analysis->IntegrationMethod()->GetAlphaF()) * analysis->IntegrationMethod()->GetNewmarkGama() * analysis->IntegrationMethod()->GetWilsonTheta() * (*analysis->DeltaT());
	double constM = 1 - analysis->IntegrationMethod()->GetAlphaM();
	double const2 = analysis->IntegrationMethod()->GetNewmarkBeta() * pow(analysis->IntegrationMethod()->GetWilsonTheta() * (*analysis->DeltaT()), 2);
	*mult2 = *mass * constM; //kept separate to be used in a subsequent iteration
	*add = (*mult2 + (*damp) * const1); //kept separate so it can be used below

	return (*add) * (1 / const2);
}

Matrix Solver::CalculateDynamicForce(const Matrix* prevDisp, const Matrix* prevVel, const Matrix* prevAcc, Matrix* add, const Matrix* damp, const Matrix* totMass, const Matrix* FInc, const Matrix* mRed, Matrix* add11, const PreAnalysisSetUp* setUp, const double* time)
{
	DynamicAnalysis* analysis = static_cast<DynamicAnalysis*>(setUp->Analysis());
	Matrix prevTerm = (*prevDisp * (1 / pow(analysis->IntegrationMethod()->GetWilsonTheta() * (*analysis->DeltaT()), 2)) + *prevVel * (1 / (analysis->IntegrationMethod()->GetWilsonTheta() * (*analysis->DeltaT()))) + *prevAcc * 0.5) * (1 / analysis->IntegrationMethod()->GetNewmarkBeta());
	Matrix firstMatrix = (*add) * prevTerm;

	*add11 = (*prevVel + *prevAcc * (analysis->IntegrationMethod()->GetWilsonTheta() * (*analysis->DeltaT())) * (1 - analysis->IntegrationMethod()->GetAlphaF())); //kept separate to be used below
	Matrix secondMatrix = *damp * *add11;

	Matrix thirdMatrix = *totMass * *prevAcc;

	Matrix load(1, 1);
	Matrix prevLoad(1, 1);

	if (analysis->Type() == AnalysisTypes::Seismic) {
		SeismicLoad* sLoad = static_cast<SeismicLoad*>(analysis->Load());
		load = (*totMass * SeismicLoad::GetSeismicLoadVector(sLoad, FInc, time)) * (-1);
		double prevTime = *time - *analysis->DeltaT();
		prevLoad = (*totMass * SeismicLoad::GetSeismicLoadVector(sLoad, FInc, &prevTime)) * (-1);
	}
	else if (analysis->Type() == AnalysisTypes::Impulse) {
		ImpulseLoad* impLoad = static_cast<ImpulseLoad*>(analysis->Load());
		load = *FInc * ImpulseLoad::LoadFromTime(impLoad, time);
		double prevTime = *time - *analysis->DeltaT();
		prevLoad = *FInc * ImpulseLoad::LoadFromTime(impLoad, &prevTime);
	}

	Matrix FInc1 = (load * (analysis->IntegrationMethod()->GetWilsonTheta() * (1 - analysis->IntegrationMethod()->GetAlphaF())));
	Matrix FInc2(FInc1.GetDimX(), 1);
	if (analysis->IntegrationMethod()->GetAlphaF() == 0) {
		FInc2 = (prevLoad * (1 - analysis->IntegrationMethod()->GetWilsonTheta()));
	}
	else {
		FInc2 = (prevLoad * (analysis->IntegrationMethod()->GetAlphaF()));
	}
	Matrix deltaF = FInc1 + FInc2;//kept separate to be used below
	Matrix FIncDyn = deltaF - *mRed * analysis->IntegrationMethod()->GetAlphaF() * *prevDisp;
	return FIncDyn + firstMatrix - secondMatrix - thirdMatrix;
}

//<summary>Calculates the complete stiffness matrix of the shell elements with threads based on the force based theory</summary>
Matrix Solver::ReducedShellStiffMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp) {
	Matrix m(*setUp->ReducedStiffMatrixSize(), *setUp->ReducedStiffMatrixSize());

	if (structManager->ShellElements()->size() != 0) {
		std::mutex matrixLock;
		std::vector<std::thread> threadList;
		std::map<int, std::vector<ShellElement*>>::const_iterator it;

		for (int i = 0; i < (*setUp->AvailableThreads() - 1); i++) {
			it = setUp->GetShellThreads()->find(i + 2); //with this 2 here I'm assuring that the first element will be run by the thread of the program. Every other element will be split into the other threads
			if (it != setUp->GetShellThreads()->end()) {
				threadList.emplace_back(ShellElement::AssembleCompleteGlobalMatrixThreads, &it->second, std::ref(m), std::ref(matrixLock));
			}
		}

		it = setUp->GetShellThreads()->find(1); //first element of the shellthreads vector
		ShellElement::AssembleCompleteGlobalMatrixThreads(&it->second, m, matrixLock);

		for (int i = 0; i < threadList.size(); i++) {
			threadList[i].join();
		}
	}
	return m;
}

Matrix Solver::CompleteShellMassMatrixThreads(const StructureManager* structManager, const PreAnalysisSetUp* setUp) {
	Matrix m(*setUp->ReducedStiffMatrixSize(), *setUp->ReducedStiffMatrixSize());

	if (structManager->ShellElements()->size() != 0) {
		std::mutex matrixLock;
		std::vector<std::thread> threadList;
		std::map<int, std::vector<ShellElement*>>::const_iterator it;

		for (int i = 0; i < (*setUp->AvailableThreads() - 1); i++) {
			it = setUp->GetShellThreads()->find(i + 2);
			if (it != setUp->GetShellThreads()->end()) {
				threadList.emplace_back(ShellElement::AssembleCompleteGlobalMassMatrixThreads, &it->second, std::ref(m), std::ref(matrixLock));
			}
		}

		it = setUp->GetShellThreads()->begin();
		ShellElement::AssembleCompleteGlobalMassMatrixThreads(&it->second, m, matrixLock);

		for (int i = 0; i < threadList.size(); i++) {
			threadList[i].join();
		}
	}

	Mass::AddExplicitMassesOnExistingMatrix(&m, structManager, setUp->DOF()); //adds any additional mass to the existing shell masses

	return m;
}

//<summary>Calculates the complete stiffness matrix of the shell elements with threads based on the force based theory</summary>
Matrix Solver::ShellRestrictedStiffMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp) {
	int sizeRow = Support::TotalDOFsRestrained(structManager->Supports());
	int sizeCol = *setUp->StiffMatrixSize();
	Matrix m(sizeRow, sizeCol);

	if (structManager->ShellElements()->size() != 0) {
		std::mutex matrixLock;
		std::vector<std::thread> threadList;
		std::map<int, std::vector<ShellElement*>>::const_iterator it;

		for (int i = 0; i < (*setUp->AvailableThreads() - 1); i++) {
			it = setUp->GetShellThreads()->find(i + 2);
			if (it != setUp->GetShellThreads()->end()) {
				threadList.emplace_back(ShellElement::AssembleCompleteRestrictedGlobalMatrixThreads, &it->second, std::ref(m), setUp->ReducedStiffMatrixSize(), std::ref(matrixLock));
			}
		}

		it = setUp->GetShellThreads()->begin();
		ShellElement::AssembleCompleteRestrictedGlobalMatrixThreads(&it->second, m, setUp->ReducedStiffMatrixSize(), matrixLock);

		for (int i = 0; i < threadList.size(); i++) {
			threadList[i].join();
		}
	}
	return m;
}

Matrix Solver::ReducedSpringStiffMatrix(const int* redSize, const std::map<int, Spring3D*>* listOfSprings, const AnalysisSpringRecorder* springRecorder) {
	if (listOfSprings->size() != 0) {
		Matrix m = Spring3D::AssembleSpringGlobalMatrixOnReducedSizedMatrix(redSize, listOfSprings, springRecorder);
		return m;
	}
	else {
		std::cout << "Error calculating Spring Stiffness Matrices" << std::endl;
		return Matrix(0);
	}
}

Matrix Solver::ReducedRestrictStiffnessMatrix(const Matrix* shellStiff, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecorder)
{
	Matrix springStiff = ReducedSpringRestrictedStiffMatrix(setUp->ReducedStiffMatrixSize(), structManager->SpringElements(), springRecorder);
	return *shellStiff + springStiff;
}

Matrix Solver::ReducedSpringRestrictedStiffMatrix(const int* redSize, const std::map<int, Spring3D*>* listOfSprings, const AnalysisSpringRecorder* springRecorder) {
	
	if (listOfSprings->size() != 0) {
		Matrix m = Spring3D::AssembleSpringGlobalRestrictedMatrixOnComplete(redSize, listOfSprings, springRecorder);
		return m;
	}
	else {
		std::cout << "Error calculating Spring Stiffness Matrices" << std::endl;
		return Matrix(0);
	}
}

void Solver::DisplacementLoadStiffness(Matrix& stiff, const std::vector<int>* dispDOFs, double* highStiff) {
	*highStiff = MatrixOperation::GetBiggestDiagTerm(stiff) * pow(10, 6);
	
	for (int i = 0; i < dispDOFs->size(); i++) {
		stiff.GetMatrixDouble()[(*dispDOFs)[i]][(*dispDOFs)[i]] += *highStiff;
	}
}

Matrix Solver::ReducedForceMatrix(const Matrix* FConst, const Matrix* FIncr, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const int* step, const double* highStiff,  const AnalysisSpringRecorder* springRecord) {
	Matrix F = (*FConst) + (*FIncr) * (*setUp->LoadFactors())[*step]; //adds the constant and incremental terms of the applied laods, considering the current loadstep
	
	if (setUp->DispLoadDOFs()->size() != 0) { //if there are displacement loads
		DisplacementLoadForce(&F, structManager, setUp, springRecord, &(*setUp->LoadFactors())[*step], highStiff); //adds the big stiffness term to account for displacement loads, if any
	}
	return F;
}

void Solver::DisplacementLoadForce(const Matrix* force, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecord, const double* loadFraction, const double* highStiff) {
	int count = 0;

	std::map<int, Support*>::const_iterator it = structManager->Supports()->begin();

	while (it != structManager->Supports()->end()) {
		for (int j = 0; j < it->second->GetSupportVector().size(); j++) {
			if (it->second->GetSupportVector()[j][1] != 0) { //no need to do this if the support is just a boundary condition
				double dispLoad = it->second->GetSupportVector()[j][1];
				force->GetMatrixDouble()[(*setUp->DispLoadDOFs())[count]][0] += (*highStiff) * dispLoad * (*loadFraction);
				count++;
			}
		}
		it++;
	}

	PlasticDisplacementLoadForce(force, structManager, setUp, springRecord);
}

void Solver::PlasticDisplacementLoadForce(const Matrix* force, const StructureManager* structManager, const PreAnalysisSetUp* setUp, const AnalysisSpringRecorder* springRecord) {
	std::map<int, std::vector<double>> vec = springRecord->GetDisplacementMap().find("plasticDisp")->second.find("new")->second;
	std::vector<int> plasticIndexes = Spring3D::GetPlasticDispIndexes(&vec, structManager, setUp->DOF());
	int count2 = 0;
	std::map<int, std::vector<double>>::iterator it = vec.begin();

	while (it != vec.end()) {
		
		std::vector<int> matPos = structManager->SpringElements()->find(it->first)->second->GetListOfGlobalMaterialDirections();
		std::map<int, std::vector<double>>::iterator it2 = vec.begin();
		for (int j = 0; j < it->second.size(); j++) {
			if (it->second[j] != 0) { //no need to do this if no plastic disp
				double dispLoad = it->second[j];
				double stiff = structManager->SpringElements()->find(it->first)->second->GetListOfMaterials()[matPos[j]]->GetSecantStiffnessFromDisplacement(springRecord->GetDisplacementMap().find("disp")->second.find("new")->second.find(it->first)->second[j],
					springRecord->GetDisplacementMap().find("plasticDisp")->second.find("new")->second.find(it->first)->second[j], 
					springRecord->GetDisplacementMap().find("maxDisp")->second.find("new")->second.find(it->first)->second[j],
					springRecord->GetDisplacementMap().find("minDisp")->second.find("new")->second.find(it->first)->second[j], 
					springRecord->GetStagesMap().find("list")->second.find(it->first)->second[j],
					springRecord->GetStagesMap().find("new")->second.find(it->first)->second[j],
					springRecord->GetDisplacementMap().find("unlDisp")->second.find("new")->second.find(it->first)->second[j], 
					springRecord->GetDisplacementMap().find("relDisp")->second.find("new")->second.find(it->first)->second[j]);
				//double force2 = listOfSpring[i].GetListOfMaterials()[matPos[j]]->GetForceFromDisplacement(dispLoad, listOfMaxDisp[i][j], listOfMinDisp[i][j]);
				double force2 = stiff * dispLoad;
				force->GetMatrixDouble()[plasticIndexes[count2]][0] += force2; //no *loadfraction since I want to apply the entire plastic disp at once
				count2++;
			}
		}
		it++;
	}
}

Matrix Solver::CalculateNaturalFrequenciesAndModeShapes(const Matrix* stiffMatrix, const Matrix* massMatrix, std::vector<double>* natFreq, const std::vector<std::vector<int>>* totalMassDOFVec, const StructureManager* structManager, const PreAnalysisSetUp* setUp) {
	
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

	natFreq->reserve(s.eigenvalues().size());

	for (int i = 0; i < s.eigenvalues().size(); i++) {
		natFreq->emplace_back(sqrt(eigenMatrix1(i).real()));
	}
	Matrix modeShapes(s.eigenvalues().size(), s.eigenvalues().size());
	for (int i = 0; i < s.eigenvalues().size(); i++) { //for each mode (column)
		for (int j = 0; j < s.eigenvalues().size(); j++) { //for each DOF (row)
			modeShapes.GetMatrixDouble()[j][i] = eigenMatrix2(j, i).real();
		}
	}

	Matrix m = Solver::GetTotalModalMatrix(&modeShapes, structManager, setUp); //Creates a matrix similar to modeShapes, but will all DOFS, not only the DOFS with masses
	return m;
}

std::vector<double> Solver::RayleighDampingConstants(const PreAnalysisSetUp* setUp, const std::vector<double>* natFreq) {

	DynamicAnalysis* analysis = static_cast<DynamicAnalysis*>(setUp->Analysis());

	double alpha1 = (2 * (*analysis->Damp2()) * (*natFreq)[(*analysis->DampMode2()) - 1] - 2 * (*analysis->Damp1()) * (*natFreq)[(*analysis->DampMode1()) - 1]) / (pow((*natFreq)[(*analysis->DampMode2()) - 1], 2) - pow((*natFreq)[(*analysis->DampMode1()) - 1], 2));
	double alpha0 = 2 * (*analysis->Damp1()) * (*natFreq)[(*analysis->DampMode1()) - 1] - alpha1 * pow((*natFreq)[(*analysis->DampMode1()) - 1], 2);

	std::vector<double> vec;
	vec.reserve(2);
	vec.emplace_back(alpha0);
	vec.emplace_back(alpha1);

	return vec;
}

Matrix Solver::RayleighDampingMatrix(const Matrix* m, const Matrix* k, const std::vector<double>* constants) {

	Matrix c = *m * (*constants)[0] + *k * (*constants)[1];

	return c;
}

Matrix Solver::GetTotalModalMatrix(const Matrix* m, const StructureManager* structManager, const PreAnalysisSetUp* setUp) {
	//This function returns the full displacement vector from the reduced version
	int count = 0;
	Matrix disp(*setUp->StiffMatrixSize(), 1);
	Matrix ans(*setUp->StiffMatrixSize(), 0); //0 is intentional

	std::vector<int> vecPos;
	std::vector<double> vecDisp;

	std::map<int, Support*>::const_iterator it = structManager->Supports()->begin();

	while (it != structManager->Supports()->end()) {//for each support
		for (int j = 0; j < it->second->GetSupportVector().size(); j++) { //this will store all the global positions and displacement values in vectors
			if (it->second->GetSupportVector()[j][1] == 0) { //only do this if the support is a boundary condition
				int nodeID = it->second->GetNode();
				int dir = it->second->GetSupportVector()[j][0];
				vecPos.emplace_back((nodeID - 1) * (*setUp->DOF()) + (dir - 1));
				vecDisp.emplace_back(0);
			}
		}
		it++;
	}

	int index = 0;
	int i = 0;
	for (int k = 0; k < m->GetDimY(); k++) {
		while (i < *setUp->StiffMatrixSize()) {
			for (int j = 0; j < 3; j++) { //only for translation DOFs
				bool notSup = !(std::find(vecPos.begin(), vecPos.end(), i + j) != vecPos.end());
				int val = i + j;
				if (notSup && (HasDOFMass(&val, structManager, setUp->DOF() ) || Mass::HasDOFAppliedMass(&val, structManager->Masses(), setUp->DOF()))) { //if it has mass
					disp.GetMatrixDouble()[i + j][0] = m->GetMatrixDouble()[index][k];
					index++;
				}
			}
			i += 6;
		}
		ShellElement::GetNinthNodeDisplacement(&disp, structManager->ShellElements(), setUp->DOF());
		ans = MatrixOperation::AddMatrixRight(ans, disp);
		i = 0;
		index = 0;
	}
	return ans;
}

bool Solver::HasDOFMass(const int* DOF, const StructureManager* structManager, const int* nDOF) {
	int count = Support::NumberOfDOFBeforeDOF(DOF, structManager->Supports(), nDOF);
	int reducedDOF = *DOF - count;
	
	if (!Support::IsDOFConstrained(DOF, structManager->Supports(), nDOF)) {
		std::map<int, ShellElement*>::const_iterator it = structManager->ShellElements()->begin();
		while (it != structManager->ShellElements()->end()) {//for each shell
			std::vector<std::vector<int>> vec = it->second->GetGlobalMassDOFVector();
			for (int j = 0; j < vec.size(); j++) {
				if (vec[j][0] == reducedDOF) {
					return false;
				}
			}
			it++;
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
