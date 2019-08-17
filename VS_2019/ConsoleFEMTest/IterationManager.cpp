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
	Matrix mRed = Solver::ReducedStiffnessMatrix(&shellStiff, structManager, setUp, &springRecord, NULL);

	//start to perform the load steps
	for (int step = 0; step < *setUp->LoadSteps(); step++) {
		Matrix Fmult = (*setUp->ConstForces() + *setUp->IncForces()) * (*setUp->LoadFactors())[step];
		Matrix F_iter = Load::GetReducedLoadMatrix(&Fmult, structManager->Supports(), setUp->DOF());
		Matrix d_iter = MatrixOperation::FullCholesky(mRed, F_iter);
		Matrix completeD_iter = Displacement::GetTotalDisplacementMatrix(d_iter, structManager, setUp);

		//this is just a recorder. It is not being used for the rest of the analysis
		//dispRecorder->Add(Displacement::GetNewNodalCoordinates(structManager->Nodes(), completeD_iter)); //this will put the new nodal coordinates on the 'nodesPerIter'
	}

	//FileOperation::SaveIterationsResult("LoadStep", dispRecorder, structManager->Nodes());
}

//Displacement-load method.
void IterationManager::PerformMatNonlinearAnalysis(const StructureManager* structManager, const PreAnalysisSetUp* setUp, std::string &fileName) {
	
	// <Start of setting up of variables useful during the analysis>
	//NodalRecorder<Node>* dispRecorder = new NodalRecorder<Node>(); //stores the new nodal coordinates after each complete iteration
	//NodalRecorder<Load>* forceRecorder = new NodalRecorder<Load>(); //stores the forces on each node after each complete iteration
	
	Matrix unorganizedDisplacement(*setUp->StiffMatrixSize(), 1); //The matrix composed of the displacements of the free nodes on top of the displacements of the restricted nodes, regardless of their DOFs
	Matrix forceMatrix(*setUp->StiffMatrixSize() - Support::TotalDOFsRestrained(structManager->Supports()), 1); //The matrix to store the resulting forces at the end of each loadstep

	// <Start of several list of spring data used to enable the localization of each srping in its cyclic material model
	AnalysisSpringRecorder springRecord('a', structManager->SpringElements());
	// <End of several list of spring data used to enable the localization of each srping in its cyclic material model

	// <End of setting up of variables useful during the analysis>

	//<Start of calculations that are not required to be performed every loadstep>
	Matrix FConstRed = Load::GetReducedLoadMatrix(setUp->ConstForces(), structManager->Supports(), setUp->DOF()); //reduced version of the matrix of constant forces
	Matrix FIncRed = Load::GetReducedLoadMatrix(setUp->IncForces(), structManager->Supports(), setUp->DOF()); //reduced version of the matrix of incremental forces

	Matrix shellStiff = Solver::ReducedShellStiffMatrix(structManager, setUp); //calcualtes the reduced stiffness amtrix only accounting for the shell elements
	Matrix shellRestricStiff = Solver::ShellRestrictedStiffMatrix(structManager, setUp); //calculates the stiffness matrix of the restricted DOFs accounting only for the shell elements

	//Perform the load steps
	for (int step = 0; step < *setUp->LoadSteps(); step++) {
		std::cout << "" << std::endl;
		std::cout << "Loadstep " << step + 1 << std::endl;
		std::cout << "Iteration " << 1 << std::endl;

		//<Start of matrix stiffness solving>
		//Create a function tthat does all these operations and return the Reduced Stiff
		Matrix mRed = Solver::ReducedStiffnessMatrix(&shellStiff, structManager, setUp, &springRecord, &_biggestStiffVal);

		//Create a function that does all these operations and return the Reduced Force
		Matrix F_iter = Solver::ReducedForceMatrix(&FConstRed, &FIncRed, structManager, setUp, &step, &_biggestStiffVal, &springRecord);
		
		Matrix d_iter = MatrixOperation::FullCholesky(mRed, F_iter); //solving the F = Kd system using Cholesky
		Matrix completeD_iter = Displacement::GetTotalDisplacementMatrix(d_iter, structManager, setUp); //displacements of all nodes
		//<End of matrix stiffness solving>

		Displacement::UpdatePositionVectorsOfSprings(&completeD_iter, structManager->SpringElements(), &springRecord, setUp->DOF()); //This function is responsible for updating all the parameters relevant tot he identification of which region the spring is at in its cyclic material model

		if (Spring3D::CheckMaterialNonlinearityConvergenceDispBased(structManager, setUp, &springRecord, &_breakAnalysis)) { //check if springs have converged
			//this means that all spring elements converged. Ready for next loadstep or iteration
			//dispRecorder->Add(Displacement::GetNewNodalCoordinates(structManager->Nodes(), completeD_iter)); //this add the increment completeD_iter to the list of nodes coord

			Spring3D::UpdateSpringLoadStages(structManager->SpringElements(), &springRecord); //Updates the list that stores the load stages of each spring
			
			//<Start of solving the restricted stiffness equations to obtain reactions>
			Matrix mRest = Solver::ReducedRestrictStiffnessMatrix(&shellRestricStiff, structManager, setUp, &springRecord);
			Matrix totalDisplacement = Displacement::GetTotalDisplacementNotOrganized(&d_iter, structManager, &springRecord.GetDisplacementMap().find("plasticDisp")->second.find("new")->second , setUp->DOF());
			Matrix reactions = mRest * totalDisplacement;
			//<End of solving the restricted stiffness equations to obtain reactions>

			Matrix completeForce = Load::GetTotalForceMatrix(&F_iter, &reactions, structManager, setUp, &_biggestStiffVal, &(*setUp->LoadFactors())[step]); //Now organized
			//forceRecorder->Add(Load::GetNewLoads(structManager->Loads(), &completeForce, setUp->DOF())); //storing the total forces on this converged loadstep in the appropriate array
			
			//<Start of updating important lists for cyclic localization>
			springRecord.UpdatePerIterDisps();
			//<End of updating important lists for cyclic localization>
		}
		else if (_breakAnalysis) { //if convergence test returns a 'breakAnalysis' flag.
			std::cout << "Analysis didn't converge at loadstep " << step + 1 << " and iteration " << 1 << std::endl;
			break; //this will exit this for loop
		}
		else { //if the convergence test failed. This means not all spring element converged. Need to do iterations in order to converge before going to next loadStep.
			bool converged = false; //flag that indicates if the analysis has converged

			for (int i = 0; i < *setUp->Iterations() - 1; i++) { //Perform the iterations
				std::cout << "Iteration " << i + 2 << std::endl;

				// <Start of stiffness matrix solving of the displacements>
				Matrix mRed = Solver::ReducedStiffnessMatrix(&shellStiff, structManager, setUp, &springRecord, &_biggestStiffVal);

				Matrix F_iter = Solver::ReducedForceMatrix(&FConstRed, &FIncRed, structManager, setUp, &step, &_biggestStiffVal, &springRecord);

				Matrix dNew = MatrixOperation::FullCholesky(mRed, F_iter);
				Matrix completeD_prime = Displacement::GetTotalDisplacementMatrix(dNew, structManager, setUp);
				// <End of stiffness matrix solving of the displacements>

				Displacement::UpdatePositionVectorsOfSprings(&completeD_prime, structManager->SpringElements(), &springRecord, setUp->DOF());

				if (Spring3D::CheckMaterialNonlinearityConvergenceDispBased(structManager, setUp, &springRecord, &_breakAnalysis)) {
					
					//dispRecorder->Add(Displacement::GetNewNodalCoordinates(structManager->Nodes(), completeD_prime)); //this add the increment completeD_iter to the list of nodes coord

					Spring3D::UpdateSpringLoadStages(structManager->SpringElements(), &springRecord);
					
					//<Start of solving the restricted stiffness equations to obtain reactions>
					Matrix mRest = Solver::ReducedRestrictStiffnessMatrix(&shellRestricStiff, structManager, setUp, &springRecord);
					Matrix totalDisplacement = Displacement::GetTotalDisplacementNotOrganized(&dNew, structManager, &springRecord.GetDisplacementMap().find("plasticDisp")->second.find("new")->second, setUp->DOF());
					Matrix reactions = mRest * totalDisplacement;
					//<End of solving the restricted stiffness equations to obtain reactions>

					Matrix completeForce = Load::GetTotalForceMatrix(&F_iter, &reactions, structManager, setUp, &_biggestStiffVal, &(*setUp->LoadFactors())[step]); //Now organized
					//forceRecorder->Add(Load::GetNewLoads(structManager->Loads(), &completeForce, setUp->DOF())); //storing the total forces on this converged loadstep in the appropriate array

					//<Start of updating important lists for cyclic localization>
					springRecord.UpdatePerIterDisps();
					//<End of updating important lists for cyclic localization>

					converged = true;
					break; //this will exit the for loop
				}

				if (_breakAnalysis) {
					std::cout << "Analysis didn't converge at loadstep " << step + 1 << " and iteration " << i + 2 << std::endl;
					break; //this will exit this for loop
				}

				//This will run if analysis did not converge in this iteration
				completeD_iter.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(completeD_prime));

				//The next two are storing important info on variables declared outside of this for loop so it can be used by outside functions
				unorganizedDisplacement = Displacement::GetTotalDisplacementNotOrganized(&dNew, structManager, &springRecord.GetDisplacementMap().find("plasticDisp")->second.find("new")->second, setUp->DOF());
				forceMatrix.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(F_iter));
			}
			//this will run if all the iterations are exausted and no convergence was found
			if (!converged) {
				std::cout << "Loadstep " << step + 1 << " did not converge!" << std::endl;
				//if the iteration converged, this line of code does not need to be run
				Spring3D::UpdateSpringLoadStages(structManager->SpringElements(), &springRecord);
				
				//<Start of solving the restricted stiffness equations to obtain reactions>
				Matrix mRest = Solver::ReducedRestrictStiffnessMatrix(&shellRestricStiff, structManager, setUp, &springRecord);
				Matrix reactions = mRest * unorganizedDisplacement;
				//<End of solving the restricted stiffness equations to obtain reactions>

				Matrix completeForce = Load::GetTotalForceMatrix(&F_iter, &reactions, structManager, setUp, &_biggestStiffVal, &(*setUp->LoadFactors())[step]); //Now organized
				//forceRecorder->Add(Load::GetNewLoads(structManager->Loads(), &completeForce, setUp->DOF())); //storing the total forces on this converged loadstep in the appropriate array

				//<Start of updating important lists for cyclic localization>
				springRecord.UpdatePerIterDisps();
				//<End of updating important lists for cyclic localization>
			}
		}

		if (_breakAnalysis) {
			break;  //this will exit this for loop
		}

		//TODO:HOW TO SAVE THE DISPLACEMENT IN THE RECORDER IF IT DOES NOT CONVERGE??

	}

	//SAVE RESULTS FILE
}

//Displacement-load method.
void IterationManager::PerformDynamicAnalysisWithIterationsMatNonlinearDispBased(const StructureManager* structManager, const PreAnalysisSetUp* setUp, std::string& fileName) {

	// <Start of setting up of variables useful during the analysis>
	//NodalRecorder<Node>* dispRecorder = new NodalRecorder<Node>(); //stores the new nodal coordinates after each complete iteration
	//NodalRecorder<Load>* forceRecorder = new NodalRecorder<Load>(); //stores the forces on each node after each complete iteration
	std::vector<Matrix> velPerStep; //stores the velocities on each node at each completed loadstep
	std::vector<Matrix> accPerStep; //stores the acceleration on each node at each completed loadstep

	Matrix redDisplacement(*setUp->ReducedStiffMatrixSize(), 1); //The matrix composed of the displacements of the free nodes on top of the displacements of the resitricted nodes, regardless of their DOFs
	Matrix forceMatrix(*setUp->ReducedStiffMatrixSize(), 1); //The matrix to store the resulting forces at the end of each loadstep

	// <Start of several list of spring data used to enable the localization of each srping in its cyclic material model
	AnalysisSpringRecorder springRecord('a', structManager->SpringElements());
	//Matrices that will store the values for the displacements, velocities and accelerations of the previous loadstep

	Matrix prevDisp(*setUp->ReducedStiffMatrixSize(), 1);
	Matrix prevVel(*setUp->ReducedStiffMatrixSize(), 1);
	Matrix prevAcc(*setUp->ReducedStiffMatrixSize(), 1);
	// <End of several list of spring data used to enable the localization of each srping in its cyclic material model
	// <End of setting up of variables useful during the analysis>

	//<Start of calculations that are not required to be performed every loadstep>
	Matrix FConstRed = Load::GetReducedLoadMatrix(setUp->ConstForces(), structManager->Supports(), setUp->DOF()); //reduced version of the matrix of constant forces
	Matrix FIncRed = Load::GetReducedLoadMatrix(setUp->IncForces(), structManager->Supports(), setUp->DOF()); //reduced version of the matrix of incremental forces

	Matrix shellStiff = Solver::ReducedShellStiffMatrix(structManager, setUp); 
	//This is the same as the first stiffness matrix
	Matrix modalTotalStiffMatrix = Solver::ReducedStiffnessMatrix(&shellStiff, structManager, setUp, &springRecord, &_biggestStiffVal);
	Matrix shellRestricStiff = Solver::ShellRestrictedStiffMatrix(structManager, setUp); //calculates the stiffness matrix of the restricted DOFs accounting only for the shell elements
	
	Matrix shellMass = Solver::CompleteShellMassMatrixThreads(structManager, setUp);

	//Matrix modalTotalStiffMatrix(MatrixOperation::CopyMatrixDouble(totalStiff), shellStiff.GetDimX(), shellStiff.GetDimY()); //copy of the total stiffness matrix that will be modified in the modal analysis calculatiosn

	std::vector<double> natFreq; //vector that will store the values of the natural frequencies
	std::vector<std::vector<int>> totalMassDOFVec = ShellElement::GetTotalGlobalMassDOFVector(structManager->ShellElements()); //gets a vector of all the DOFS that DO NOT Have mass. It is used everytime the matrices are required to be reduced only to the DOFS that have mass
	Matrix totalModes = Solver::CalculateNaturalFrequenciesAndModeShapes(&modalTotalStiffMatrix, &shellMass, &natFreq, &totalMassDOFVec, structManager, setUp); //Creates a matrix similar to modeShapes, but will all DOFS, not only the DOFS with masses
	
	std::vector<double> rayleighConstants = Solver::RayleighDampingConstants(setUp, &natFreq);; //vector to store the rayleigh constants
	
	//Perform the load steps
	for (int step = 0; step < *setUp->LoadSteps(); step++) {
		std::cout << "" << std::endl;
		std::cout << "Loadstep " << step + 1 << std::endl;
		std::cout << "Iteration " << 1 << std::endl;

		DynamicAnalysis* analysis = static_cast<DynamicAnalysis*>(setUp->Analysis());
		double time = (step + 1) * *analysis->DeltaT();

		//<Start of matrix stiffness solving>
		//<Set up stiffness matrix>
		Matrix mRed = Solver::ReducedStiffnessMatrix(&shellStiff, structManager, setUp, &springRecord, &_biggestStiffVal);
		//<End Set up stiffness Matrix>

		//<Set up damping matrix>
		Matrix damp(mRed.GetDimX());
		if (*analysis->Damp1() != 0 || *analysis->Damp2() != 0) {
			Matrix damp = Solver::RayleighDampingMatrix(&shellMass, &mRed, &rayleighConstants);
		}
		//<End set up damping matrix>

		//<Set dynamic 'stiff' matrix>
		Matrix add(*setUp->ReducedStiffMatrixSize());
		Matrix mult2(*setUp->ReducedStiffMatrixSize());
		Matrix kDyn = Solver::ReducedDynamicStiffnessMatrix(setUp, &shellMass, &damp, &add, &mult2);
		Matrix totStiff = mRed * (1 - analysis->IntegrationMethod()->GetAlphaF()) + kDyn;
		//<End set dynamic stiff matrix>

		//<Set up dynamic force>
		Matrix add11(*setUp->ReducedStiffMatrixSize(), 1);
		Matrix FDyn = Solver::CalculateDynamicForce(&prevDisp, &prevVel, &prevAcc, &add, &damp, &shellMass, &FIncRed, &mRed, &add11, setUp, &time);
		//<End of set up dynamic force>


		Matrix F_iter = FConstRed + FDyn; //adds the constant and incremental terms of the applied laods, considering the current loadstep
		Solver::PlasticDisplacementLoadForce(&F_iter, structManager, setUp, &springRecord);
		Matrix d_iter = MatrixOperation::FullCholesky(totStiff, F_iter); //solving the F = Kd system using Cholesky
		//<End of matrix stiffness solving>

		//<Calculate aceleration and velocity>
		Matrix redPrevDisp = ShellElement::ReducedAccelerationForceMatrix(&prevDisp, &totalMassDOFVec, setUp->DOF());
		Matrix redVel = ShellElement::ReducedAccelerationForceMatrix(&prevVel, &totalMassDOFVec, setUp->DOF());
		Matrix redAcc = ShellElement::ReducedAccelerationForceMatrix(&prevAcc, &totalMassDOFVec, setUp->DOF());
		Matrix acc(1, 1);
		Matrix curVel(1, 1);
		Matrix completeD_iter;
		IterationManager::CalculateAccelerationAndVelocity(structManager, setUp, &d_iter, &totalMassDOFVec, &redPrevDisp, &redVel, &redAcc, &prevAcc, &acc, &curVel, &prevVel, &completeD_iter);
		//<End calcualtion of acc and vel>
		
		Displacement::UpdatePositionVectorsOfSprings(&completeD_iter, structManager->SpringElements(), &springRecord, setUp->DOF()); //This function is responsible for updating all the parameters relevant tot he identification of which region the spring is at in its cyclic material model

		if (Spring3D::CheckMaterialNonlinearityConvergenceDispBased(structManager, setUp, &springRecord, &_breakAnalysis)) { //check if springs have converged
			//this means that all spring elements converged. Ready for next loadstep or iteration
			//dispRecorder->Add(Displacement::GetNewNodalCoordinates(structManager->Nodes(), completeD_iter)); //this add the increment completeD_iter to the list of nodes coord

			Spring3D::UpdateSpringLoadStages(structManager->SpringElements(), &springRecord); //Updates the list that stores the load stages of each spring

			//<Start of solving the restricted stiffness equations to obtain reactions>
			Matrix mRest = Solver::ReducedRestrictStiffnessMatrix(&shellRestricStiff, structManager, setUp, &springRecord);
			Matrix totalDisplacement = Displacement::GetTotalDisplacementNotOrganized(&d_iter, structManager, &springRecord.GetDisplacementMap().find("plasticDisp")->second.find("new")->second, setUp->DOF());
			Matrix reactions = mRest * totalDisplacement;
			//<End of solving the restricted stiffness equations to obtain reactions>

			Matrix completeForce = Load::GetTotalForceMatrix(&F_iter, &reactions, structManager, setUp, &_biggestStiffVal, &(*setUp->LoadFactors())[step]); //Now organized
			//forceRecorder->Add(Load::GetNewLoads(structManager->Loads(), &completeForce, setUp->DOF())); //storing the total forces on this converged loadstep in the appropriate array
			velPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(curVel), curVel.GetDimX(), curVel.GetDimY());
			accPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(acc), acc.GetDimX(), acc.GetDimY());

			//<Start of updating important lists for cyclic localization>
			springRecord.UpdatePerIterDisps();
			//<End of updating important lists for cyclic localization>

			prevDisp.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(d_iter));
			prevVel.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(curVel));
			prevAcc.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(acc));
		
		}
		else if (_breakAnalysis) { //if convergence test returns a 'breakAnalysis' flag.
			std::cout << "Analysis didn't converge at loadstep " << step + 1 << " and iteration " << 1 << std::endl;
			break; //this will exit this for loop
		}
		else { //if the convergence test failed. This means not all spring element converged. Need to do iterations in order to converge before going to next loadStep.
			bool converged = false; //flag that indicates if the analysis has converged
			//F_iter.DestroyMatrixDouble(); //this need to be destroyed before the next iteration because the previous matrix have different displacement factors.


			for (int i = 0; i < *setUp->Iterations() - 1; i++) { //Perform the iterations
				std::cout << "Iteration " << i + 2 << std::endl;

				// <Start of stiffness matrix solving of the displacements>
				Matrix mRed = Solver::ReducedStiffnessMatrix(&shellStiff, structManager, setUp, &springRecord, &_biggestStiffVal);
				
				//<Set up damping matrix>
				Matrix damp(mRed.GetDimX());
				if (*analysis->Damp1() != 0 || *analysis->Damp2() != 0) {
					Matrix damp = Solver::RayleighDampingMatrix(&shellMass, &mRed, &rayleighConstants);
				}
				//<End set up damping matrix>
				
				//<Set dynamic 'stiff' matrix>
				Matrix kDyn = Solver::ReducedDynamicStiffnessMatrix(setUp, &shellMass, &damp, &add, &mult2);
				Matrix totStiff = mRed * (1 - analysis->IntegrationMethod()->GetAlphaF()) + kDyn;
				//<End set dynamic stiff matrix>

				//<Set up dynamic force>
				Matrix FDyn = Solver::CalculateDynamicForce(&prevDisp, &prevVel, &prevAcc, &add, &damp, &shellMass, &FIncRed, &mRed, &add11, setUp, &time);
				//<End of set up dynamic force>

				Matrix F_iter = FConstRed + FDyn;
				Solver::PlasticDisplacementLoadForce(&F_iter, structManager, setUp, &springRecord);
				Matrix dNew = MatrixOperation::FullCholesky(totStiff, F_iter);
				// <End of stiffness matrix solving of the displacements>

				//<Calculate aceleration and velocity>
				Matrix completeD_prime;
				IterationManager::CalculateAccelerationAndVelocity(structManager, setUp, &dNew, &totalMassDOFVec, &redPrevDisp, &redVel, &redAcc, &prevAcc, &acc, &curVel, &prevVel, &completeD_prime);

				//<End calcualtion of acc and vel>

				Displacement::UpdatePositionVectorsOfSprings(&completeD_prime, structManager->SpringElements(), &springRecord, setUp->DOF()); //This function is responsible for updating all the parameters relevant tot he identification of which region the spring is at in its cyclic material model

				if (Spring3D::CheckMaterialNonlinearityConvergenceDispBased(structManager, setUp, &springRecord, &_breakAnalysis)) {
					//dispRecorder->Add(Displacement::GetNewNodalCoordinates(structManager->Nodes(), completeD_prime)); //this add the increment completeD_iter to the list of nodes coord

					Spring3D::UpdateSpringLoadStages(structManager->SpringElements(), &springRecord); //Updates the list that stores the load stages of each spring

					//<Start of solving the restricted stiffness equations to obtain reactions>
					Matrix mRest = Solver::ReducedRestrictStiffnessMatrix(&shellRestricStiff, structManager, setUp, &springRecord);
					Matrix totalDisplacement = Displacement::GetTotalDisplacementNotOrganized(&dNew, structManager, &springRecord.GetDisplacementMap().find("plasticDisp")->second.find("new")->second, setUp->DOF());
					Matrix reactions = mRest * totalDisplacement;
					//<End of solving the restricted stiffness equations to obtain reactions>

					Matrix completeForce = Load::GetTotalForceMatrix(&F_iter, &reactions, structManager, setUp, &_biggestStiffVal, &(*setUp->LoadFactors())[step]); //Now organized
					//forceRecorder->Add(Load::GetNewLoads(structManager->Loads(), &completeForce, setUp->DOF())); //storing the total forces on this converged loadstep in the appropriate array
					velPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(curVel), curVel.GetDimX(), curVel.GetDimY());
					accPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(acc), acc.GetDimX(), acc.GetDimY());

					//<Start of updating important lists for cyclic localization>
					springRecord.UpdatePerIterDisps();
					//<End of updating important lists for cyclic localization>

					prevDisp.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(d_iter));
					prevVel.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(curVel));
					prevAcc.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(acc));

					converged = true;
					break; //this will exit the for loop
				}

				if (_breakAnalysis) {
					std::cout << "Analysis didn't converge at loadstep " << step + 1 << " and iteration " << i + 2 << std::endl;
					break; //this will exit this for loop
				}

				//This will run if analysis did not converge in this iteration
				//completeD_iter.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(completeD_prime));

				//The next two are storing important info on variables declared outside of this for loop so it can be used by outside functions
				redDisplacement.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(dNew));
				forceMatrix.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(F_iter));
			}
			if (!converged) {
				std::cout << "Loadstep " << step + 1 << " did not converge!" << std::endl;
				//if the iteration converged, this line of code does not need to be run
				//dispRecorder->Add(Displacement::GetNewNodalCoordinates(structManager->Nodes(), completeD_prime)); //this add the increment completeD_iter to the list of nodes coord

				Spring3D::UpdateSpringLoadStages(structManager->SpringElements(), &springRecord); //Updates the list that stores the load stages of each spring

				//<Start of solving the restricted stiffness equations to obtain reactions>
				Matrix mRest = Solver::ReducedRestrictStiffnessMatrix(&shellRestricStiff, structManager, setUp, &springRecord);
				Matrix totalDisplacement = Displacement::GetTotalDisplacementNotOrganized(&redDisplacement, structManager, &springRecord.GetDisplacementMap().find("plasticDisp")->second.find("new")->second, setUp->DOF());
				Matrix reactions = mRest * totalDisplacement;
				//<End of solving the restricted stiffness equations to obtain reactions>

				Matrix completeForce = Load::GetTotalForceMatrix(&F_iter, &reactions, structManager, setUp, &_biggestStiffVal, &(*setUp->LoadFactors())[step]); //Now organized
				//forceRecorder->Add(Load::GetNewLoads(structManager->Loads(), &completeForce, setUp->DOF())); //storing the total forces on this converged loadstep in the appropriate array
				velPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(curVel), curVel.GetDimX(), curVel.GetDimY());
				accPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(acc), acc.GetDimX(), acc.GetDimY());

				//<Start of updating important lists for cyclic localization>
				springRecord.UpdatePerIterDisps();
				//<End of updating important lists for cyclic localization>

				prevDisp.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(d_iter));
				prevVel.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(curVel));
				prevAcc.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(acc));
			}
		}

		if (_breakAnalysis) {
			break;  //this will exit this for loop
		}

		//TODO:HOW TO SAVE THE DISPLACEMENT IN THE RECORDER IF IT DOES NOT CONVERGE??
	}

	//FileOperation::SaveResultsFile(fileName, nodesPerStep, originalList, forcePerStep, natFreq, totalModes);
}

void IterationManager::TESTPerformDynamicAnalysisWithIterationsMatNonlinearDispBased() {

	// <Start of setting up of variables useful during the analysis>
	double deltaT = 0.01;
	double totalTime = 50;
	int totalSteps = totalTime / deltaT;
	// <Start of several list of spring data used to enable the localization of each srping in its cyclic material model


	Matrix prevDisp(1, 1);
	Matrix prevVel(1, 1);
	Matrix prevAcc(1, 1);

	// <End of several list of spring data used to enable the localization of each srping in its cyclic material model

	// <End of setting up of variables useful during the analysis>

	//<Start of calculations that are not required to be performed every loadstep>

	Matrix k(1);
	Matrix m(1);
	k.GetMatrixDouble()[0][0] = 1000;
	m.GetMatrixDouble()[0][0] = 5000;

	std::vector<double> natFreq;

	Matrix sqrtMass = MatrixOperation::Sqrt(m);
	Matrix invSqrtMass = MatrixOperation::GetInverse(sqrtMass);

	Matrix mult1 = invSqrtMass * k;
	Matrix mult2 = mult1 * invSqrtMass;

	Eigen::MatrixXd eigenMatrix;
	MatrixOperation::ConvertToEigenMatrix(mult2, eigenMatrix);

	Eigen::EigenSolver<Eigen::MatrixXd> s(eigenMatrix); // the instance s(A) includes the eigensystem

	natFreq.reserve(s.eigenvalues().size());

	for (int i = 0; i < s.eigenvalues().size(); i++) {
		natFreq.emplace_back(sqrt(s.eigenvalues()(i).real()));
	}

	Matrix modeShapes;
	modeShapes.SetDimensions(s.eigenvalues().size(), s.eigenvalues().size());
	for (int i = 0; i < s.eigenvalues().size(); i++) { //for each mode (column)
		for (int j = 0; j < s.eigenvalues().size(); j++) { //for each DOF (row)
			modeShapes.GetMatrixDouble()[j][i] = s.eigenvectors()(j, i).real();
		}
	}

	/*
	std::vector<double> rayleighConstants;
	rayleighConstants.reserve(2);
	int mode1 = 1, mode2 = 2; //for now
	double damp1 = 0, damp2 = 0; //for now
	rayleighConstants = Solver::RayleighDampingConstants(mode1, damp1, mode2, damp2, natFreq);
	*/

	double NewMarkGama = 1.0 / 2;
	double NewMarkBeta = 1.0 / 4;
	double WilsonTheta = 1.0;

 	std::vector<double> disp;
	disp.reserve(totalSteps);
	std::vector<double> force;
	force.reserve(totalSteps);
	std::vector<double> velVec;
	velVec.reserve(totalSteps);
	std::vector<double> accVec;
	accVec.reserve(totalSteps);

	//Perform the load steps
	for (int step = 0; step < totalSteps; step++) {
		std::cout << "" << std::endl;
		std::cout << "Loadstep " << step + 1 << std::endl;
		std::cout << "Iteration " << 1 << std::endl;

		//<Start of matrix stiffness solving>
		//<Set up stiffness matrix>

		//<End Set up stiffness Matrix>


		//<Set up damping matrix>
		//Matrix damp = Solver::RayleighDampingMatrix(m, k, rayleighConstants);,
		Matrix damp(1);
		//<End set up damping matrix>

		//<Set dynamic 'stiff' matrix>
		double const1 = NewMarkGama * WilsonTheta * deltaT;
		double const2 = NewMarkBeta * pow(WilsonTheta * deltaT, 2);
		Matrix mult1 = damp * const1;
		Matrix add = m + mult1;
		Matrix kDyn = add * (1 / const2);
		//<End set dynamic stiff matrix>

		//<Set up dynamic force>
		Matrix first = add * (1 / NewMarkBeta);
		Matrix multFirst1 = prevDisp * (1 / pow(WilsonTheta * deltaT, 2));
		Matrix multFirst2 = prevVel * (1 / (WilsonTheta * deltaT));
		Matrix multFirst3 = prevAcc * 0.5;
		Matrix add1 = multFirst1 + multFirst2;
		Matrix add2 = add1 + multFirst3;
		Matrix firstMatrix = first * add2;

		Matrix multSecond1 = prevAcc * (WilsonTheta * deltaT);
		Matrix add11 = prevVel + multSecond1;
		Matrix secondMatrix = damp * add11;

		Matrix thirdMatrix = m * prevAcc;

		Matrix FInc1(1);
		FInc1.GetMatrixDouble()[0][0] = (Load::SampleDynamicForceFunction(10, 0.44721359549995793, 0, (step + 1) * deltaT) * WilsonTheta);
		Matrix FInc2(1);
		FInc2.GetMatrixDouble()[0][0] = (Load::SampleDynamicForceFunction(10, 0.44721359549995793, 0, (step)* deltaT) * (1 - WilsonTheta));
		Matrix FInc = FInc1 - FInc2;

		Matrix FDyn1 = FInc + firstMatrix;
		Matrix FDyn2 = FDyn1 - secondMatrix;
		Matrix FDyn = FDyn2 - thirdMatrix;

		//<End of set up dynamic force>

		//double fraction = Load::DefineLoadFractionAtLoadStep(type, step, totalSteps, stepsPerPeak, peakInc, cyclesPerPeak, iniPeak);
		Matrix F_iter = FDyn; //adds the constant and incremental terms of the applied laods, considering the current loadstep
		Matrix totStiff = k + kDyn;
		Matrix d_iter = MatrixOperation::FullCholesky(totStiff, F_iter); //solving the F = Kd system using Cholesky
		//<End of matrix stiffness solving>

		disp.emplace_back(d_iter.GetMatrixDouble()[0][0]);
		force.emplace_back(FInc1.GetMatrixDouble()[0][0]);

		//<Calculate aceleration and velocity>
		Matrix vel1 = d_iter - prevDisp;
		Matrix velFirst = vel1 * (NewMarkGama / (NewMarkBeta * deltaT));
		Matrix velSecond = prevVel * (NewMarkGama / NewMarkBeta);
		Matrix vecThird = prevAcc * ((NewMarkGama * deltaT / (2 * NewMarkBeta)) - deltaT);
		Matrix curVel1 = velFirst - velSecond;
		Matrix curVel2 = curVel1 - vecThird;
		Matrix curVel = curVel2 + prevVel;

		Matrix acc1 = FInc1;
		Matrix acc2 = damp * curVel;
		Matrix acc3 = k * d_iter;
		Matrix Facc1 = acc1 + acc2;
		Matrix Facc = Facc1 - acc3;
		Matrix acc = MatrixOperation::FullCholesky(m, Facc);

		//<End calcualtion of acc and vel>

		velVec.emplace_back(curVel.GetMatrixDouble()[0][0]);
		accVec.emplace_back(acc.GetMatrixDouble()[0][0]);

		prevDisp.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(d_iter));
		prevVel.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(curVel));
		prevAcc.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(acc));

	}

	std::ofstream myfile;
	myfile.open("ResultTEST.txt");
	for (int i = 0; i < disp.size(); i++) {
		myfile << disp[i] << " " << force[i] << std::endl;	
	}
	myfile.close();

	std::ofstream myfile2;
	myfile2.open("ResultTESTAccVel.txt");
	for (int i = 0; i < accVec.size(); i++) {
		myfile2 << velVec[i] << " " << accVec[i] << std::endl;
	}
	myfile2.close();

	std::cout << "End of Analysis" << std::endl;
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

void IterationManager::CalculateAccelerationAndVelocity(const StructureManager* structManager, const PreAnalysisSetUp* setUp, Matrix* dNew, const std::vector<std::vector<int>>* totalMassDOFVec, const Matrix* redPrevDisp, const Matrix* redVel,
	const Matrix* redAcc, const Matrix* prevAcc, Matrix* acc, Matrix* curVel, const Matrix* prevVel, Matrix* completeD_prime){
	
	DynamicAnalysis* analysis = static_cast<DynamicAnalysis*>(setUp->Analysis());

	Matrix redDisp = ShellElement::ReducedAccelerationForceMatrix(dNew, totalMassDOFVec, setUp->DOF());
	Matrix deltaAccWilson = (redDisp - *redPrevDisp) * (1 / (analysis->IntegrationMethod()->GetNewmarkBeta() * pow(analysis->IntegrationMethod()->GetWilsonTheta()  * (*analysis->DeltaT()), 2))) - *redVel * (1 / (analysis->IntegrationMethod()->GetNewmarkBeta() * analysis->IntegrationMethod()->GetWilsonTheta()  * (*analysis->DeltaT()))) - *redAcc * (1 / (2 * analysis->IntegrationMethod()->GetNewmarkBeta()));
	Matrix accRed = deltaAccWilson * (1 / analysis->IntegrationMethod()->GetWilsonTheta()) + *redAcc;
	*acc = ShellElement::ConvertAccFromReducedToTotal(&accRed, totalMassDOFVec, setUp->ReducedStiffMatrixSize());

	Matrix velRed = *redAcc * (*analysis->DeltaT()) + (accRed - *redAcc) * analysis->IntegrationMethod()->GetNewmarkGama() * (*analysis->DeltaT()) + *redVel;
	//Matrix velRed = (redDisp - redPrevDisp) * (NewMarkGama / (NewMarkBeta * deltaT)) - redVel * (NewMarkGama / NewMarkBeta) - redAcc * ((deltaT * NewMarkGama / (NewMarkBeta * 2)) - deltaT) + redVel;
	*curVel = ShellElement::ConvertAccFromReducedToTotal(&velRed, totalMassDOFVec, setUp->ReducedStiffMatrixSize());
	Matrix redNewDisp = *redVel * (*analysis->DeltaT()) + *redAcc * (pow((*analysis->DeltaT()), 2) / 2) + (accRed - *redAcc) * analysis->IntegrationMethod()->GetNewmarkBeta() * pow((*analysis->DeltaT()), 2) + *redPrevDisp;
	Matrix disp = ShellElement::ConvertAccFromReducedToTotal(&redNewDisp, totalMassDOFVec, setUp->ReducedStiffMatrixSize());
	*dNew = MatrixOperation::FillMatrixBasedOnOtherMatrix(&disp, dNew);
	*completeD_prime = Displacement::GetTotalDisplacementMatrix(*dNew, structManager, setUp); //displacements of all nodes
}

IterationManager::~IterationManager()
{
}
