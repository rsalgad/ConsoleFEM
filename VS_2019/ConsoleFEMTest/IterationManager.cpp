#include "pch.h"
#include "Solver.h"
#include "IterationManager.h"
#include "FileOperation.h"
#include <map>
#include <fstream>

double IterationManager::_biggestStiffVal = 0;
bool IterationManager::_breakAnalysis = false;

IterationManager::IterationManager()
{

}

//For elastic analysis only (with loadsteps, if wanted) working
void IterationManager::PerformElasticAnalysis(const StructureManager* structManager, const PreAnalysisSetUp* setUp, int nLoadSteps, std::string &fileName) {

	//NodalRecorder<Node>* dispRecorder = new NodalRecorder<Node>(structManager->Nodes());
	AnalysisSpringRecorder springRecord('a', structManager->SpringElements());

	//Perform the initial iteration (will occur even when no iterations are specified)
	//The analysis is always performed on the original, undeformed, configuration of the structure, but with increasing load.
	Matrix shellStiff = Solver::ReducedShellStiffMatrix(structManager, setUp); //calcualtes the reduced stiffness amtrix only accounting for the shell elements
	Matrix mRed = Solver::ReducedStiffnessMatrix(shellStiff, structManager, setUp, &springRecord, NULL);

	//start to perform the load steps
	for (int step = 0; step < *setUp->LoadSteps(); step++) {
		Matrix Fmult = (*setUp->ConstForces() + *setUp->IncForces()) * (*setUp->LoadFactors())[step];
		Matrix F_iter = Load::GetReducedLoadMatrix(Fmult, structManager->Supports(), setUp->DOF());
		Matrix d_iter = MatrixOperation::FullCholesky(mRed, F_iter);
		Matrix completeD_iter = Displacement::GetTotalDisplacementMatrix(d_iter, structManager, setUp);

		//this is just a recorder. It is not being used for the rest of the analysis
		//dispRecorder->Add(Displacement::GetNewNodalCoordinates(structManager->Nodes(), completeD_iter)); //this will put the new nodal coordinates on the 'nodesPerIter'
	}

	//FileOperation::SaveIterationsResult("LoadStep", dispRecorder, structManager->Nodes());
}

void IterationManager::PerformDynamicAnalysisWithIterationsMatNonlinearDispBased(const StructureManager* structManager, const PreAnalysisSetUp* setUp, std::string& fileName) {
	// <Start of setting up of variables useful during the analysis>
	NodalRecorder<Node> dispRecorder(structManager->Nodes()); //stores the new nodal coordinates after each complete iteration
	NodalRecorder<Load> forceRecorder(structManager->Nodes()); //stores the forces on each node after each complete iteration
	//Only used in dynamic analyses
	std::vector<Matrix> velPerStep; //stores the velocities on each node at each completed loadstep
	std::vector<Matrix> accPerStep; //stores the acceleration on each node at each completed loadstep

	// <Start of several list of spring data used to enable the localization of each srping in its cyclic material model
	AnalysisSpringRecorder springRecord('a', structManager->SpringElements());
	// <End of several list of spring data used to enable the localization of each srping in its cyclic material model

	//Only used in dynamic analyses
	Matrix prevDisp(*setUp->ReducedStiffMatrixSize(), 1);
	Matrix prevVel(*setUp->ReducedStiffMatrixSize(), 1);
	Matrix prevAcc(*setUp->ReducedStiffMatrixSize(), 1);
	std::vector<double> natFreq;
	std::vector<std::vector<int>> totalMassDOFVec;
	std::vector<double> rayleighConstants;
	Matrix totalModes;
	Matrix shellMass;
	// <End of setting up of variables useful during the analysis>

	//<Start of calculations that are not required to be performed every loadstep>
	Matrix FConstRed = Load::GetReducedLoadMatrix(*setUp->ConstForces(), structManager->Supports(), setUp->DOF()); //reduced version of the matrix of constant forces
	Matrix FIncRed = Load::GetReducedLoadMatrix(*setUp->IncForces(), structManager->Supports(), setUp->DOF()); //reduced version of the matrix of incremental forces
	Matrix shellStiff = Solver::ReducedShellStiffMatrix(structManager, setUp); //calcualtes the reduced stiffness amtrix only accounting for the shell elements
	Matrix shellRestricStiff = Solver::ShellRestrictedStiffMatrix(structManager, setUp); //calculates the stiffness matrix of the restricted DOFs accounting only for the shell elements

	bool isDynamic = (setUp->Analysis()->Type() == AnalysisTypes::Seismic || setUp->Analysis()->Type() == AnalysisTypes::Impulse);
	if (isDynamic) {

		shellMass = Solver::CompleteShellMassMatrixThreads(structManager, setUp);
		totalMassDOFVec = ShellElement::GetTotalGlobalMassDOFVector(structManager->ShellElements()); //gets a vector of all the DOFS that DO NOT Have mass. It is used everytime the matrices are required to be reduced only to the DOFS that have mass
		Matrix modalTotalStiffMatrix = Solver::ReducedStiffnessMatrix(shellStiff, structManager, setUp, &springRecord, &_biggestStiffVal);
		totalModes = Solver::CalculateNaturalFrequenciesAndModeShapes(modalTotalStiffMatrix, shellMass, &natFreq, &totalMassDOFVec, structManager, setUp); //Creates a matrix similar to modeShapes, but will all DOFS, not only the DOFS with masses
		rayleighConstants = Solver::RayleighDampingConstants(setUp, &natFreq);; //vector to store the rayleigh constants
	}

	//Perform the load steps
	for (int step = 0; step < *setUp->LoadSteps(); step++) {
		std::cout << "" << std::endl;
		std::cout << "Loadstep " << step + 1 << std::endl;
		std::cout << "Iteration " << 1 << std::endl;

		DynamicAnalysis* analysis = nullptr;
		double time;

		if (isDynamic) {
			analysis = static_cast<DynamicAnalysis*>(setUp->Analysis());
			time = (step + 1) * *analysis->DeltaT();
		}

		//<Start of matrix stiffness solving>
		//Create a function tthat does all these operations and return the Reduced Stiff
		Matrix mRed = Solver::ReducedStiffnessMatrix(shellStiff, structManager, setUp, &springRecord, &_biggestStiffVal);

		//<Set up damping matrix>
		Matrix damp(mRed.GetDimX());
		if (isDynamic) {
			if (*analysis->Damp1() != 0 || *analysis->Damp2() != 0) {
				damp = Solver::RayleighDampingMatrix(shellMass, mRed, &rayleighConstants);
			}
		}
		//<End set up damping matrix>

		//<Set dynamic 'stiff' matrix>
		Matrix add, add11, mult2, FDyn, totStiff;
		if (isDynamic) {
			Matrix kDyn = Solver::ReducedDynamicStiffnessMatrix(setUp, shellMass, damp, add, mult2);
			totStiff = mRed * (1 - analysis->IntegrationMethod()->GetAlphaF()) + kDyn;
			//<End set dynamic stiff matrix>

			//<Set up dynamic force>
			FDyn = Solver::CalculateDynamicForce(prevDisp, prevVel, prevAcc, add, damp, shellMass, FIncRed, mRed, add11, setUp, &time);
			//<End of set up dynamic force>
		}

		Matrix F_iter, d_iter;
		if (isDynamic) {
			F_iter = FConstRed + FDyn; //adds the constant and incremental terms of the applied laods, considering the current loadstep
			Solver::PlasticDisplacementLoadForce(F_iter, structManager, setUp, &springRecord);
			d_iter = MatrixOperation::FullCholesky(totStiff, F_iter); //solving the F = Kd system using Cholesky

		}
		else {
			F_iter = Solver::ReducedForceMatrix(FConstRed, FIncRed, structManager, setUp, &step, &_biggestStiffVal, &springRecord);
			d_iter = MatrixOperation::FullCholesky(mRed, F_iter); //solving the F = Kd system using Cholesky
		}

		Matrix redPrevDisp, redVel, redAcc, acc, curVel, completeD_iter;
		if (isDynamic) {

			redPrevDisp = ShellElement::ReducedAccelerationForceMatrix(prevDisp, &totalMassDOFVec, setUp->DOF());
			redVel = ShellElement::ReducedAccelerationForceMatrix(prevVel, &totalMassDOFVec, setUp->DOF());
			redAcc = ShellElement::ReducedAccelerationForceMatrix(prevAcc, &totalMassDOFVec, setUp->DOF());
			completeD_iter;
			IterationManager::CalculateAccelerationAndVelocity(structManager, setUp, d_iter, &totalMassDOFVec, redPrevDisp, redVel, redAcc, prevAcc, acc, curVel, prevVel, completeD_iter);

		}
		else {
			completeD_iter = Displacement::GetTotalDisplacementMatrix(d_iter, structManager, setUp); //displacements of all nodes
		}
		//<End of matrix stiffness solving>

		Displacement::UpdatePositionVectorsOfSprings(&completeD_iter, structManager->SpringElements(), &springRecord, setUp->DOF()); //This function is responsible for updating all the parameters relevant tot he identification of which region the spring is at in its cyclic material model

		if (Spring3D::CheckMaterialNonlinearityConvergenceDispBased(structManager, setUp, &springRecord, &_breakAnalysis)) { //check if springs have converged
			//this means that all spring elements converged. Ready for next loadstep or iteration
			dispRecorder.Add(Displacement::GetNewNodalCoordinates(structManager->Nodes(), completeD_iter)); //this add the increment completeD_iter to the list of nodes coord

			Spring3D::UpdateSpringLoadStages(structManager->SpringElements(), &springRecord); //Updates the list that stores the load stages of each spring

			//<Start of solving the restricted stiffness equations to obtain reactions>
			Matrix mRest = Solver::ReducedRestrictStiffnessMatrix(shellRestricStiff, structManager, setUp, &springRecord);
			Matrix totalDisplacement = Displacement::GetTotalDisplacementNotOrganized(d_iter, structManager, *springRecord.GetNewPlasticDisp(), setUp->DOF());
			Matrix reactions = mRest * totalDisplacement;
			//<End of solving the restricted stiffness equations to obtain reactions>

			Matrix completeForce = Load::GetTotalForceMatrix(F_iter, reactions, structManager, setUp, &_biggestStiffVal, &(*setUp->LoadFactors())[step]); //Now organized
			forceRecorder.Add(Load::GetNewLoads(forceRecorder.Nodes(), completeForce, setUp->DOF())); //storing the total forces on this converged loadstep in the appropriate array
			if (isDynamic) {
				velPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(curVel), curVel.GetDimX(), curVel.GetDimY());
				accPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(acc), acc.GetDimX(), acc.GetDimY());
			}

			//<Start of updating important lists for cyclic localization>
			springRecord.UpdatePerIterDisps();
			//<End of updating important lists for cyclic localization>

			if (isDynamic) {
				prevDisp.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(d_iter));
				prevVel.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(curVel));
				prevAcc.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(acc));
			}
		}
		else if (_breakAnalysis) { //if convergence test returns a 'breakAnalysis' flag.
			std::cout << "Analysis didn't converge at loadstep " << step + 1 << " and iteration " << 1 << std::endl;
			break; //this will exit this for loop
		}
		else { //if the convergence test failed. This means not all spring element converged. Need to do iterations in order to converge before going to next loadStep.

			for (int i = 0; i < *setUp->Iterations() - 1; i++) { //Perform the iterations
				std::cout << "Iteration " << i + 2 << std::endl;

				// <Start of stiffness matrix solving of the displacements>
				mRed = Solver::ReducedStiffnessMatrix(shellStiff, structManager, setUp, &springRecord, &_biggestStiffVal);

				if (isDynamic) {
					if (*analysis->Damp1() != 0 || *analysis->Damp2() != 0) {
						damp = Solver::RayleighDampingMatrix(shellMass, mRed, &rayleighConstants);
					}
				}

				if (isDynamic) {
					Matrix kDyn = Solver::ReducedDynamicStiffnessMatrix(setUp, shellMass, damp, add, mult2);
					totStiff = mRed * (1 - analysis->IntegrationMethod()->GetAlphaF()) + kDyn;
					//<End set dynamic stiff matrix>

					//<Set up dynamic force>
					FDyn = Solver::CalculateDynamicForce(prevDisp, prevVel, prevAcc, add, damp, shellMass, FIncRed, mRed, add11, setUp, &time);
					//<End of set up dynamic force>
				}

				if (isDynamic) {
					F_iter = FConstRed + FDyn; //adds the constant and incremental terms of the applied laods, considering the current loadstep
					Solver::PlasticDisplacementLoadForce(F_iter, structManager, setUp, &springRecord);
					d_iter = MatrixOperation::FullCholesky(totStiff, F_iter);
				}
				else {
					F_iter = Solver::ReducedForceMatrix(FConstRed, FIncRed, structManager, setUp, &step, &_biggestStiffVal, &springRecord);
					d_iter = MatrixOperation::FullCholesky(mRed, F_iter);
				}

				if (isDynamic) {
					IterationManager::CalculateAccelerationAndVelocity(structManager, setUp, d_iter, &totalMassDOFVec, redPrevDisp, redVel, redAcc, prevAcc, acc, curVel, prevVel, completeD_iter);
				}
				else {
					completeD_iter = Displacement::GetTotalDisplacementMatrix(d_iter, structManager, setUp); //displacements of all nodes
				}
				// <End of stiffness matrix solving of the displacements>

				Displacement::UpdatePositionVectorsOfSprings(&completeD_iter, structManager->SpringElements(), &springRecord, setUp->DOF());

				if (Spring3D::CheckMaterialNonlinearityConvergenceDispBased(structManager, setUp, &springRecord, &_breakAnalysis) || (i == *setUp->Iterations() - 2)) { //the last part of this check will go inside this if even if the convergence is not achieved BUT it is the last iteration

					if (i == *setUp->Iterations() - 2) {
						std::cout << "Loadstep " << step + 1 << " did not converge!" << std::endl;
					}

					dispRecorder.Add(Displacement::GetNewNodalCoordinates(structManager->Nodes(), completeD_iter)); //this add the increment completeD_iter to the list of nodes coord

					Spring3D::UpdateSpringLoadStages(structManager->SpringElements(), &springRecord);

					//<Start of solving the restricted stiffness equations to obtain reactions>
					Matrix mRest = Solver::ReducedRestrictStiffnessMatrix(shellRestricStiff, structManager, setUp, &springRecord);
					Matrix totalDisplacement = Displacement::GetTotalDisplacementNotOrganized(d_iter, structManager, *springRecord.GetNewPlasticDisp(), setUp->DOF());
					Matrix reactions = mRest * totalDisplacement;
					//<End of solving the restricted stiffness equations to obtain reactions>

					Matrix completeForce = Load::GetTotalForceMatrix(F_iter, reactions, structManager, setUp, &_biggestStiffVal, &(*setUp->LoadFactors())[step]); //Now organized
					forceRecorder.Add(Load::GetNewLoads(forceRecorder.Nodes(), completeForce, setUp->DOF())); //storing the total forces on this converged loadstep in the appropriate array

					if (isDynamic) {
						velPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(curVel), curVel.GetDimX(), curVel.GetDimY());
						accPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(acc), acc.GetDimX(), acc.GetDimY());
					}

					//<Start of updating important lists for cyclic localization>
					springRecord.UpdatePerIterDisps();
					//<End of updating important lists for cyclic localization>

					if (isDynamic) {
						prevDisp.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(d_iter));
						prevVel.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(curVel));
						prevAcc.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(acc));
					}

					break; //this will exit the for loop
				}

				if (_breakAnalysis) {
					std::cout << "Analysis didn't converge at loadstep " << step + 1 << " and iteration " << i + 2 << std::endl;
					break; //this will exit this for loop
				}
			}
		}

		if (_breakAnalysis) {
			break;  //this will exit this for loop
		}
	}

	FileOperation::SaveResultsFile(fileName, structManager, &dispRecorder, &forceRecorder, natFreq, totalModes);
}

bool IterationManager::CheckConvergenceCriteria(Matrix &D, double limit)
{
	double CF = 0;
	int size = D.GetDimX();
	double sum = 0;
	for (int i = 0; i < size; i++) {
		sum += pow((D.GetMatrixDouble()[i][0]), 2);
	}
	sum = sum / size;
	CF = sqrt(sum);
	std::cout << CF << std::endl;
	if (CF < limit) {
		return true;
	}
	else {
		return false;
	}
}


//abandoned for now
/*
void IterationManager::PerformAnalysisWithIterationsGeomNonlinear(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Spring3D> &listOfSprings, std::vector<Load> &listOfLoads, std::vector<Support> &listOfSups, int nIter, int nLoadSteps) {
	std::vector<Node> originalList = listOfNodes;

	std::vector<std::vector<Node>> nodesPerIter;
	nodesPerIter.reserve(nLoadSteps);
	double div = (1.0 / nLoadSteps);
	std::vector<std::vector<ShellElement>> shellElemVecs = Solver::SetUpThreadsForShells(listOfShells, std::thread::hardware_concurrency());


	//Perform the load steps
	for (int step = 0; step < nLoadSteps; step++) {
		Matrix* mRed = &Solver::CompleteStiffnessMatrixWithThreadsDispBased(listOfNodes, listOfShells, listOfSprings, listOfSups, std::thread::hardware_concurrency(), shellElemVecs);
		Matrix* F = &Load::AssembleLoadMatrix(listOfNodes, listOfLoads);
		Matrix Fmult = *F * div;
		Matrix* F_iter = &Load::GetReducedLoadMatrix(Fmult, listOfSups);
		Matrix* Cholesky = &Matrix::CholeskyDecomposition(*mRed);
		Matrix* CholeskyTransp = &Matrix::Transpose(*Cholesky);
		Matrix* interMatrix = &Matrix::ForwardSubstitution(*Cholesky, *F_iter);
		Matrix* d_iter = &Matrix::BackSubstitution(*CholeskyTransp, *interMatrix);
		Matrix* completeD_iter = &Displacement::GetTotalDisplacementMatrix(*d_iter, listOfSups, listOfNodes);

		std::vector<Node> newNodes = Displacement::GetNewNodalCoordinates(listOfNodes, *completeD_iter);
		std::copy(newNodes.begin(), newNodes.end(), listOfNodes.begin());

		//Perform the iterations
		for (int i = 0; i < nIter - 1; i++) {
			Matrix* mRedNew = &Solver::CompleteStiffnessMatrixWithThreadsDispBased(listOfNodes, listOfShells, listOfSprings, listOfSups, std::thread::hardware_concurrency(), shellElemVecs);
			Matrix F_prime = (*mRedNew) * (*d_iter);
			Matrix FNew = *F_iter - F_prime;
			Matrix* Cholesky_prime = &Matrix::CholeskyDecomposition(*mRedNew);
			Matrix* CholeskyTransp_prime = &Matrix::Transpose(*Cholesky_prime);
			Matrix* interMatrix_prime = &Matrix::ForwardSubstitution(*Cholesky_prime, FNew);
			Matrix* dNew = &Matrix::BackSubstitution(*CholeskyTransp_prime, *interMatrix_prime);
			*d_iter += *dNew;
			Matrix* completeD_prime = &Displacement::GetTotalDisplacementMatrix(*dNew, listOfSups, listOfNodes);
			*completeD_iter += *completeD_prime;

			std::vector<Node> newNodes2 = Displacement::GetNewNodalCoordinates(listOfNodes, *completeD_prime);
			std::copy(newNodes2.begin(), newNodes2.end(), listOfNodes.begin());

			if (CheckConvergenceCriteria(*completeD_prime, 0.00001) == true) {
				Cholesky_prime->DestroyMatrixDouble();
				CholeskyTransp_prime->DestroyMatrixDouble();
				interMatrix_prime->DestroyMatrixDouble();
				mRedNew->DestroyMatrixDouble();
				F_prime.DestroyMatrixDouble();
				FNew.DestroyMatrixDouble();
				dNew->DestroyMatrixDouble();
				completeD_prime->DestroyMatrixDouble();

				break; //this will exit the for loop
			}

			Cholesky_prime->DestroyMatrixDouble();
			CholeskyTransp_prime->DestroyMatrixDouble();
			interMatrix_prime->DestroyMatrixDouble();
			completeD_prime->DestroyMatrixDouble();
			mRedNew->DestroyMatrixDouble();
			F_prime.DestroyMatrixDouble();
			FNew.DestroyMatrixDouble();
			dNew->DestroyMatrixDouble();
		}

		nodesPerIter.emplace_back(listOfNodes);

		mRed->DestroyMatrixDouble();
		F->DestroyMatrixDouble();
		Fmult.DestroyMatrixDouble();
		Cholesky->DestroyMatrixDouble();
		CholeskyTransp->DestroyMatrixDouble();
		interMatrix->DestroyMatrixDouble();
		F_iter->DestroyMatrixDouble();
		d_iter->DestroyMatrixDouble();
		completeD_iter->DestroyMatrixDouble();
	}
}
*/

/*
void IterationManager::PostConvergenceProcedures(std::vector<Spring3D> &listOfSprings, std::vector<std::vector<std::string>> &listOfSpringLoadingStages, std::vector<std::vector<std::string>> &newSpringStages,
	Matrix &shellRestricStiff, std::vector<Support> &listOfSups, std::vector<std::vector<double>> &newListOfDisps, std::vector<std::vector<double>> &listOfMinDisps, std::vector<std::vector<double>> &listOfMaxDisps,
	std::vector<std::vector<double>> &newListOfPlasticDisps, std::vector<std::vector<double>> &unlDispPerIter, std::vector<std::vector<double>> &relDispPerIter, Matrix &dNew, std::vector<Node> &listOfNodes,
	Matrix& F_iter, double highStiff, std::vector<Matrix> &forcePerStep, std::vector<Matrix> &velPerStep, std::vector<Matrix> &accPerStep, Matrix &curVel, Matrix &acc, std::vector<std::vector<double>> &maxDispPerIter, 
	std::vector<std::vector<double>> &minDispPerIter, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp, Matrix &prevDisp, Matrix &prevVel, Matrix &prevAcc) {
	
	dispRecorder->Add(Displacement::GetNewNodalCoordinates(structManager->Nodes(), completeD_iter)); //this add the increment completeD_iter to the list of nodes coord

	Spring3D::UpdateSpringLoadStages(listOfSprings, listOfSpringLoadingStages, newSpringStages); //Updates the list that stores the load stages of each spring

	//<Start of solving the restricted stiffness equations to obtain reactions>
	Matrix mRest = Solver::ReducedRestrictStiffnessMatrix(&shellRestricStiff, structManager, setUp, &springRecord);
	Matrix totalDisplacement = Displacement::GetTotalDisplacementNotOrganized(&d_iter, structManager, &springRecord.GetDisplacementMap().find("plasticDisp")->second.find("new")->second, setUp->DOF());
	Matrix reactions = mRest * totalDisplacement;
	//<End of solving the restricted stiffness equations to obtain reactions>

	Matrix completeForce = Load::GetTotalForceMatrix(&F_iter, &reactions, structManager, setUp, &_biggestStiffVal, &(*setUp->LoadFactors())[step]); //Now organized
	forceRecorder->Add(Load::GetNewLoads(structManager->Loads(), &completeForce, setUp->DOF())); //storing the total forces on this converged loadstep in the appropriate array
	velPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(curVel), curVel.GetDimX(), curVel.GetDimY());
	accPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(acc), acc.GetDimX(), acc.GetDimY());

	//<Start of updating important lists for cyclic localization>
	springRecord.UpdatePerIterDisps();
	//<End of updating important lists for cyclic localization>

	prevDisp.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(dNew));
	prevVel.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(curVel));
	prevAcc.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(acc));
}
*/

void IterationManager::CalculateAccelerationAndVelocity(const StructureManager* structManager, const PreAnalysisSetUp* setUp, Matrix& dNew, const std::vector<std::vector<int>>* totalMassDOFVec, Matrix& redPrevDisp, Matrix& redVel,
	Matrix& redAcc, Matrix& prevAcc, Matrix& acc, Matrix& curVel, Matrix& prevVel, Matrix& completeD_prime){
	
	DynamicAnalysis* analysis = static_cast<DynamicAnalysis*>(setUp->Analysis());

	Matrix redDisp = ShellElement::ReducedAccelerationForceMatrix(dNew, totalMassDOFVec, setUp->DOF());
	Matrix deltaAccWilson = (redDisp - redPrevDisp) * (1 / (analysis->IntegrationMethod()->GetNewmarkBeta() * pow(analysis->IntegrationMethod()->GetWilsonTheta()  * (*analysis->DeltaT()), 2))) - redVel * (1 / (analysis->IntegrationMethod()->GetNewmarkBeta() * analysis->IntegrationMethod()->GetWilsonTheta()  * (*analysis->DeltaT()))) - redAcc * (1 / (2 * analysis->IntegrationMethod()->GetNewmarkBeta()));
	Matrix accRed = deltaAccWilson * (1 / analysis->IntegrationMethod()->GetWilsonTheta()) + redAcc;
	acc = ShellElement::ConvertAccFromReducedToTotal(accRed, totalMassDOFVec, setUp->ReducedStiffMatrixSize());

	Matrix velRed = redAcc * (*analysis->DeltaT()) + (accRed - redAcc) * analysis->IntegrationMethod()->GetNewmarkGama() * (*analysis->DeltaT()) + redVel;
	//Matrix velRed = (redDisp - redPrevDisp) * (NewMarkGama / (NewMarkBeta * deltaT)) - redVel * (NewMarkGama / NewMarkBeta) - redAcc * ((deltaT * NewMarkGama / (NewMarkBeta * 2)) - deltaT) + redVel;
	curVel = ShellElement::ConvertAccFromReducedToTotal(velRed, totalMassDOFVec, setUp->ReducedStiffMatrixSize());
	Matrix redNewDisp = redVel * (*analysis->DeltaT()) + redAcc * (pow((*analysis->DeltaT()), 2) / 2) + (accRed - redAcc) * analysis->IntegrationMethod()->GetNewmarkBeta() * pow((*analysis->DeltaT()), 2) + redPrevDisp;
	Matrix disp = ShellElement::ConvertAccFromReducedToTotal(redNewDisp, totalMassDOFVec, setUp->ReducedStiffMatrixSize());
	dNew = MatrixOperation::FillMatrixBasedOnOtherMatrix(disp, dNew);
	completeD_prime = Displacement::GetTotalDisplacementMatrix(dNew, structManager, setUp); //displacements of all nodes
}

IterationManager::~IterationManager()
{
}
