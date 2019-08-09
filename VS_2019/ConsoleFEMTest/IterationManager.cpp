#include "pch.h"
#include "IterationManager.h"
#include "Node.h"
#include "ShellElement.h"
#include "Spring3D.h"
#include "Load.h"
#include "Support.h"
#include "Solver.h"
#include "Displacement.h"
#include "FileOperation.h"
#include "MatrixOperation.h"
#include "TimeIntegrationMethod.h"
#include "Recorder.h"
#include "PreAnalysisSetUp.h"
#include <iostream>
#include <fstream>


IterationManager::IterationManager()
{
}

//For elastic analysis only (with loadsteps, if wanted) working
void IterationManager::PerformAnalysisWithIterations(const StructureManager* structManager, const PreAnalysisSetUp* setUp, int nLoadSteps, std::string &fileName) {

	Recorder<Node*>* dispRecorder = new Recorder<Node*>();
	
	double div = (1.0 / nLoadSteps);

	//Perform the initial iteration (will occur even when no iterations are specified)
	//The analysis is always performed on the original, undeformed, configuration of the structure, but with increasing load.
	Matrix mRed = Solver::ReducedStiffnessMatrix(structManager, setUp);
	Matrix F = Load::AssembleLoadMatrix(structManager, setUp);

	//start to perform the load steps
	for (int step = 0; step < nLoadSteps; step++) {
		Matrix Fmult = F * div * (step + 1);
		Matrix F_iter = Load::GetReducedLoadMatrix(Fmult, structManager->Supports(), setUp->DOF());
		Matrix Cholesky = MatrixOperation::CholeskyDecomposition(mRed);
		Matrix CholeskyTransp = MatrixOperation::Transpose(Cholesky);
		Matrix interMatrix = MatrixOperation::ForwardSubstitution(Cholesky, F_iter);
		Matrix d_iter = MatrixOperation::BackSubstitution(CholeskyTransp, interMatrix);
		Matrix completeD_iter = Displacement::GetTotalDisplacementMatrix(d_iter, structManager, setUp);

		//this is just a recorder. It is not being used for the rest of the analysis
		dispRecorder->Add(Displacement::GetNewNodalCoordinates(structManager->Nodes(), completeD_iter)); //this will put the new nodal coordinates on the 'nodesPerIter'
	}

	FileOperation::SaveIterationsResult("LoadStep", dispRecorder, structManager->Nodes());
}

//Displacement-load method.
void IterationManager::PerformAnalysisWithIterationsMatNonlinearDispBased(std::vector<Node> &listOfNodes, std::vector<ShellElement> &listOfShells, std::vector<Spring3D> &listOfSprings, 
																		  std::vector<Load> &listOfLoads, std::vector<Support> &listOfSups, int nIter, int nLoadSteps, std::string type,
																		  double cyclicRepeat, int stepsPerPeak, double peakInc, int cyclesPerPeak, double iniPeak, std::string &fileName) {
	
	// <Start of setting up of variables useful during the analysis>
	std::vector<Node> originalList = listOfNodes; //Gets a copy of the list of nodes in their original positions in case they get changed during the analysis
	int DOF = 6; //number of DOFs considered
	std::vector<std::vector<Node>> nodesPerStep; //stores the new nodal coordinates after each complete iteration
	std::vector<Matrix> forcePerStep; //stores the forces on each node at each completed loadstep
	int totalSteps = Load::DefineTotalLoadSteps(type, nLoadSteps, cyclicRepeat, 0, 0);
	nodesPerStep.reserve(totalSteps); //reserve the amount of memory that I know it will consume
	bool breakAnalysis = false; //indicates if the analysis should be stopped
	Matrix unorganizedDisplacement(listOfNodes.size()*DOF, 1); //The matrix composed of the displacements of the free nodes on top of the displacements of the resitricted nodes, regardless of their DOFs
	Matrix forceMatrix(listOfNodes.size()*DOF - Support::TotalDOFsRestrained(listOfSups), 1); //The matrix to store the resulting forces at the end of each loadstep

	// <Start of several list of spring data used to enable the localization of each srping in its cyclic material model
	std::vector<std::vector<double>> oldListOfDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> newListOfDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> listOfMinDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> listOfMaxDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> oldListOfMinDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> oldListOfMaxDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> newListOfPlasticDisps; //Used to update the stiffness of the springs
	std::vector<std::vector<double>> oldListOfPlasticDisps; //Used to update the stiffness of the springs
	std::vector<std::vector<double>> oldListOfUnlDisp; //Use to check when to conenct between unload and reload branchs
	std::vector<std::vector<double>> oldListOfRelDisp;
	std::vector<std::vector<double>> listOfUnlDisp; //Use to check when to conenct between unload and reload branchs
	std::vector<std::vector<double>> listOfRelDisp;
	std::vector<std::vector<double>> maxDispPerIter;
	std::vector<std::vector<double>> minDispPerIter;
	std::vector<std::vector<double>> unlDispPerIter;
	std::vector<std::vector<double>> relDispPerIter;
	std::vector<std::vector<std::string>> oldSpringStages;
	std::vector<std::vector<std::string>> newSpringStages;
	std::vector<std::vector<std::string>> listOfSpringLoadingStages;
	// <End of several list of spring data used to enable the localization of each srping in its cyclic material model

	Displacement::ZeroOutPositionVectorsOfSprings(oldListOfDisps, newListOfDisps, listOfSprings, listOfMinDisps, listOfMaxDisps, oldListOfMinDisps, oldListOfMaxDisps, oldListOfPlasticDisps, newListOfPlasticDisps, listOfSpringLoadingStages, listOfUnlDisp, listOfRelDisp, oldListOfUnlDisp, oldListOfRelDisp, maxDispPerIter, minDispPerIter, unlDispPerIter, relDispPerIter, oldSpringStages, newSpringStages);
	// <End of setting up of variables useful during the analysis>

	//<Start of calculations that are not required to be performed every loadstep>
	Matrix Fconst = Load::AssembleLoadMatrixWithFlag(listOfNodes, listOfLoads, "constant"); //assembles the matrix of constant forces
	Matrix Finc = Load::AssembleLoadMatrixWithFlag(listOfNodes, listOfLoads, "increment"); //assembles the matrix of incremental forces
	std::vector<int> indexes = Load::IdentifyIncrementalLoads(listOfLoads); //Identify which DOFs are incremental loads
	Matrix FConstRed = Load::GetReducedLoadMatrix(Fconst, listOfSups); //reduced version of the matrix of constant forces
	Matrix FIncRed = Load::GetReducedLoadMatrix(Finc, listOfSups); //reduced version of the matrix of incremental forces
	
	std::vector<std::vector<ShellElement>> shellElemVecs = Solver::SetUpThreadsForShells(listOfShells, std::thread::hardware_concurrency()); //separates the total shell elements amongst the available threads
	Matrix shellStiff = Solver::CompleteShellStiffMatrixThreads(listOfNodes, listOfShells, listOfSups, std::thread::hardware_concurrency(), shellElemVecs); //calcualtes the reduced stiffness amtrix only accounting for the shell elements
	Matrix shellRestricStiff = Solver::CompleteShellRestrictedStiffMatrixThreads(listOfNodes, listOfShells, listOfSups, std::thread::hardware_concurrency(), shellElemVecs); //calculates the stiffness matrix of the restricted DOFs accounting only for the shell elements

	//Perform the load steps
	for (int step = 0; step < totalSteps; step++) {
		std::cout << "" << std::endl;
		std::cout << "Loadstep " << step + 1 << std::endl;
		std::cout << "Iteration " << 1 << std::endl;

		//<Start of matrix stiffness solving>
		Matrix m(shellStiff.GetDimX(), shellStiff.GetDimY()); //Initializes the 'total' reduced matrix
		Solver::CompleteSpringStiffMatrixThreadsDispBasedAfterShells(listOfSprings, m, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter); //this will update 'm' with the values of the spring elements
		Matrix mRed = m + shellStiff; //complete reduced total stiffness matrix
		double highStiff = MatrixOperation::GetBiggestDiagTerm(mRed);
		Solver::DisplacementLoadStiffness(mRed, listOfSups, highStiff); //adds the big stiffness term to account for displacement loads, if any
		double fraction = Load::DefineLoadFractionAtLoadStep(type, step, totalSteps, stepsPerPeak, peakInc, cyclesPerPeak, iniPeak);
		Matrix F_iter = FConstRed + FIncRed * fraction; //adds the constant and incremental terms of the applied laods, considering the current loadstep
		Solver::DisplacementLoadForce(F_iter, listOfSups, newListOfPlasticDisps, listOfSprings, listOfMinDisps, listOfMaxDisps, highStiff, fraction, newListOfDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter); //adds the big stiffness term times the displacement amount to the force matrix
		Matrix d_iter = MatrixOperation::FullCholesky(mRed, F_iter); //solving the F = Kd system using Cholesky
		Matrix completeD_iter = Displacement::GetTotalDisplacementMatrix(d_iter, listOfSups, listOfNodes); //displacements of all nodes
		//<End of matrix stiffness solving>

		Displacement::UpdatePositionVectorsOfSprings(oldListOfDisps, newListOfDisps, completeD_iter, listOfSprings, listOfMinDisps, listOfMaxDisps, oldListOfMinDisps, oldListOfMaxDisps, oldListOfPlasticDisps, newListOfPlasticDisps, listOfSpringLoadingStages, listOfUnlDisp, listOfRelDisp, oldListOfUnlDisp, oldListOfRelDisp, maxDispPerIter, minDispPerIter, unlDispPerIter,relDispPerIter, oldSpringStages, newSpringStages); //This function is responsible for updating all the parameters relevant tot he identification of which region the spring is at in its cyclic material model

		std::vector<Node> newNodes = Displacement::GetNewNodalCoordinates(originalList, completeD_iter); //this add the increment completeD_iter to the list of nodes coord
		if (Spring3D::CheckMaterialNonlinearityConvergenceDispBased(listOfSprings, oldListOfDisps, newListOfDisps, listOfMaxDisps, listOfMinDisps, oldListOfMaxDisps, oldListOfMinDisps, oldListOfPlasticDisps, newListOfPlasticDisps, listOfSpringLoadingStages, 0.001, breakAnalysis, unlDispPerIter, relDispPerIter, oldListOfUnlDisp, oldListOfRelDisp, oldSpringStages, newSpringStages)) { //check if springs have converged
			//this means that all spring elements converged. Ready for next loadstep or iteration
			Spring3D::UpdateSpringLoadStages(listOfSprings, listOfSpringLoadingStages, newSpringStages); //Updates the list that stores the load stages of each spring
			
			//<Start of solving the restricted stiffness equations to obtain reactions>
			Matrix m2(shellRestricStiff.GetDimX(), shellRestricStiff.GetDimY());
			Solver::CompleteSpringRestrictedStiffMatrixThreadsDispBasedAfterShells(listOfSprings, m2, listOfSups, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter);
			Matrix mRest = shellRestricStiff + m2;
			Matrix totalDisplacement = Displacement::GetTotalDisplacementNotOrganized(d_iter, listOfSups, listOfNodes, listOfSprings, newListOfPlasticDisps);
			Matrix reactions = mRest * totalDisplacement;
			//<End of solving the restricted stiffness equations to obtain reactions>

			Matrix totalForce = Load::GetTotalForceNotOrganized(F_iter, reactions, listOfSups, listOfNodes); //"Not organized" means that the order on the matrix is not the correct DOF order
			Matrix completeForce = Load::GetTotalForceMatrix(totalForce, listOfSups, listOfNodes, highStiff, fraction); //Now organized
			forcePerStep.emplace_back(MatrixOperation::CopyMatrixDouble(completeForce), completeForce.GetDimX(), completeForce.GetDimY()); //storing the total forces on this converged loadstep in the appropriate array
			
			//<Start of updating important lists for cyclic localization>
			maxDispPerIter = listOfMaxDisps;
			minDispPerIter = listOfMinDisps;
			unlDispPerIter = listOfUnlDisp;
			relDispPerIter = listOfRelDisp;
			//<End of updating important lists for cyclic localization>
		}
		else if (breakAnalysis) { //if convergence test returns a 'breakAnalysis' flag.
			std::cout << "Analysis didn't converge at loadstep " << step + 1 << " and iteration " << 1 << std::endl;
			break; //this will exit this for loop
		}
		else { //if the convergence test failed. This means not all spring element converged. Need to do iterations in order to converge before going to next loadStep.
			bool converged = false; //flag that indicates if the analysis has converged
			//F_iter.DestroyMatrixDouble(); //this need to be destroyed before the next iteration because the previous matrix have different displacement factors.


			for (int i = 0; i < nIter - 1; i++) { //Perform the iterations
				std::cout << "Iteration " << i + 2 << std::endl;

				// <Start of stiffness matrix solving of the displacements>
				Matrix m(shellStiff.GetDimX(), shellStiff.GetDimY());
				Solver::CompleteSpringStiffMatrixThreadsDispBasedAfterShells(listOfSprings, m, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter);
				Matrix mRed = m + shellStiff;
				double highStiff = MatrixOperation::GetBiggestDiagTerm(mRed);
				Solver::DisplacementLoadStiffness(mRed, listOfSups, highStiff); //adds the big stiffness term to account for displacement laods
				Matrix F_iter = FConstRed + FIncRed * fraction;
				Solver::DisplacementLoadForce(F_iter, listOfSups, newListOfPlasticDisps, listOfSprings, listOfMinDisps, listOfMaxDisps, highStiff, fraction, newListOfDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter); //adds the big stiffness term times the displacement amount to the force matrix
				Matrix dNew = MatrixOperation::FullCholesky(mRed, F_iter);
				Matrix completeD_prime = Displacement::GetTotalDisplacementMatrix(dNew, listOfSups, listOfNodes);
				// <End of stiffness matrix solving of the displacements>

				Displacement::UpdatePositionVectorsOfSprings(oldListOfDisps, newListOfDisps, completeD_prime, listOfSprings, listOfMinDisps, listOfMaxDisps, oldListOfMinDisps, oldListOfMaxDisps, oldListOfPlasticDisps, newListOfPlasticDisps, listOfSpringLoadingStages, listOfUnlDisp, listOfRelDisp, oldListOfUnlDisp, oldListOfRelDisp, maxDispPerIter, minDispPerIter, unlDispPerIter, relDispPerIter, oldSpringStages, newSpringStages);

				newNodes = Displacement::GetNewNodalCoordinates(originalList, completeD_prime); //this add the increment completeD_iter to the list of nodes coord
				if (Spring3D::CheckMaterialNonlinearityConvergenceDispBased(listOfSprings, oldListOfDisps, newListOfDisps, listOfMaxDisps, listOfMinDisps, oldListOfMaxDisps, oldListOfMinDisps, oldListOfPlasticDisps, newListOfPlasticDisps, listOfSpringLoadingStages, 0.001, breakAnalysis, unlDispPerIter, relDispPerIter, oldListOfUnlDisp, oldListOfRelDisp, oldSpringStages, newSpringStages)) {
					
					Spring3D::UpdateSpringLoadStages(listOfSprings, listOfSpringLoadingStages, newSpringStages);
					
					//<Start of solving the restricted stiffness equations to obtain reactions>
					Matrix m2(shellRestricStiff.GetDimX(), shellRestricStiff.GetDimY());
					Solver::CompleteSpringRestrictedStiffMatrixThreadsDispBasedAfterShells(listOfSprings, m2, listOfSups, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter);
					Matrix mRest = shellRestricStiff + m2;
					Matrix totalDisplacement = Displacement::GetTotalDisplacementNotOrganized(dNew, listOfSups, listOfNodes, listOfSprings, newListOfPlasticDisps);
					Matrix reactions = mRest * totalDisplacement;
					//<End of solving the restricted stiffness equations to obtain reactions>

					Matrix totalForce = Load::GetTotalForceNotOrganized(F_iter, reactions, listOfSups, listOfNodes);
					Matrix completeForce = Load::GetTotalForceMatrix(totalForce, listOfSups, listOfNodes, highStiff, fraction);
					forcePerStep.emplace_back(MatrixOperation::CopyMatrixDouble(completeForce), completeForce.GetDimX(), completeForce.GetDimY());

					//<Start of updating important lists for cyclic localization>
					maxDispPerIter = listOfMaxDisps;
					minDispPerIter = listOfMinDisps;
					unlDispPerIter = listOfUnlDisp;
					relDispPerIter = listOfRelDisp;
					//<End of updating important lists for cyclic localization>

					converged = true;
					break; //this will exit the for loop
				}

				if (breakAnalysis) {
					std::cout << "Analysis didn't converge at loadstep " << step + 1 << " and iteration " << i + 2 << std::endl;
					break; //this will exit this for loop
				}

				//This will run if analysis did not converge in this iteration
				completeD_iter.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(completeD_prime));

				//The next two are storing important info on variables declared outside of this for loop so it can be used by outside functions
				unorganizedDisplacement = Displacement::GetTotalDisplacementNotOrganized(dNew, listOfSups, listOfNodes, listOfSprings, newListOfPlasticDisps);
				forceMatrix.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(F_iter));
			}
			if (!converged) {
				std::cout << "Loadstep " << step + 1 << " did not converge!" << std::endl;
				//if the iteration converged, this line of code does not need to be run
				Spring3D::UpdateSpringLoadStages(listOfSprings, listOfSpringLoadingStages, newSpringStages);
				
				//<Start of solving the restricted stiffness equations to obtain reactions>
				Matrix m2(shellRestricStiff.GetDimX(), shellRestricStiff.GetDimY());
				Solver::CompleteSpringRestrictedStiffMatrixThreadsDispBasedAfterShells(listOfSprings, m2, listOfSups, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter);
				Matrix mRest = shellRestricStiff + m2;
				Matrix reactions = mRest * unorganizedDisplacement;
				//<End of solving the restricted stiffness equations to obtain reactions>

				Matrix totalForce = Load::GetTotalForceNotOrganized(forceMatrix, reactions, listOfSups, listOfNodes);
				Matrix completeForce = Load::GetTotalForceMatrix(totalForce, listOfSups, listOfNodes, highStiff, fraction);
				forcePerStep.emplace_back(MatrixOperation::CopyMatrixDouble(completeForce), completeForce.GetDimX(), completeForce.GetDimY());
				
				//<Start of updating important lists for cyclic localization>
				maxDispPerIter = listOfMaxDisps;
				minDispPerIter = listOfMinDisps;
				unlDispPerIter = listOfUnlDisp;
				relDispPerIter = listOfRelDisp;
				//<End of updating important lists for cyclic localization>
			}
		}

		if (breakAnalysis) {
			break;  //this will exit this for loop
		}

		nodesPerStep.emplace_back(newNodes);
	}

	FileOperation::SaveIterationsResult("LoadStep", nodesPerStep, originalList);
	FileOperation::SaveIterationsForceResult("ForceStep", forcePerStep);
}

//Displacement-load method.
void IterationManager::PerformDynamicAnalysisWithIterationsMatNonlinearDispBased(StructureManager structManager, SeismicLoad &sLoad, 
																			     ImpulseLoad &impLoad, int nIter, int nLoadSteps, std::string type, std::string& fileName) {

	// <Start of setting up of variables useful during the analysis>
	std::vector<Node*> originalList = listOfNodes; //Gets a copy of the list of nodes in their original positions in case they get changed during the analysis
	int DOF = 6; //number of DOFs considered
	std::vector<std::vector<Node>> nodesPerStep; //stores the new nodal coordinates after each complete iteration
	std::vector<Matrix> forcePerStep; //stores the forces on each node at each completed loadstep
	std::vector<Matrix> velPerStep; //stores the velocities on each node at each completed loadstep
	std::vector<Matrix> accPerStep; //stores the acceleration on each node at each completed loadstep
	double reducedDOF = listOfNodes.size()*DOF - Support::TotalDOFsRestrained(listOfSups);
	double deltaT;
	double totalTime;
	if (type == "seismic") {
		deltaT = 0.005;
		totalTime = sLoad.GetTime()[0] * sLoad.GetRecords()[0].size();
	}
	else if (type == "impulse") {
		deltaT = 0.0005;
		totalTime = impLoad.GetPoints()[impLoad.GetPoints().size() - 1][0] + 1.93;
	}
	else {
		deltaT = 0.01;
		totalTime = 10;
	}

	int totalSteps = Load::DefineTotalLoadSteps(type, nLoadSteps, 0, totalTime, deltaT);
	nodesPerStep.reserve(totalSteps); //reserve the amount of memory that I know it will consume
	bool breakAnalysis = false; //indicates if the analysis should be stopped due to critical non-convergence
	Matrix redDisplacement(reducedDOF, 1); //The matrix composed of the displacements of the free nodes on top of the displacements of the resitricted nodes, regardless of their DOFs
	Matrix forceMatrix(reducedDOF, 1); //The matrix to store the resulting forces at the end of each loadstep

	// <Start of several list of spring data used to enable the localization of each srping in its cyclic material model
	std::vector<std::vector<double>> oldListOfDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> newListOfDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> listOfMinDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> listOfMaxDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> oldListOfMinDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> oldListOfMaxDisps; //Used to update the forces on the springs inside a substep
	std::vector<std::vector<double>> newListOfPlasticDisps; //Used to update the stiffness of the springs
	std::vector<std::vector<double>> oldListOfPlasticDisps; //Used to update the stiffness of the springs
	std::vector<std::vector<double>> oldListOfUnlDisp; //Use to check when to conenct between unload and reload branchs
	std::vector<std::vector<double>> oldListOfRelDisp;
	std::vector<std::vector<double>> listOfUnlDisp; //Use to check when to conenct between unload and reload branchs
	std::vector<std::vector<double>> listOfRelDisp;
	std::vector<std::vector<double>> maxDispPerIter;
	std::vector<std::vector<double>> minDispPerIter;
	std::vector<std::vector<double>> unlDispPerIter;
	std::vector<std::vector<double>> relDispPerIter;
	std::vector<std::vector<std::string>> oldSpringStages;
	std::vector<std::vector<std::string>> newSpringStages;
	std::vector<std::vector<std::string>> listOfSpringLoadingStages;

	//Matrices that will store the values for the displacements, velocities and accelerations of the previous loadstep
	Matrix prevDisp(reducedDOF, 1);
	Matrix prevVel(reducedDOF, 1);
	Matrix prevAcc(reducedDOF, 1);
	// <End of several list of spring data used to enable the localization of each srping in its cyclic material model

	//Initialize the values of the lists to zero or a string value. 
	Displacement::ZeroOutPositionVectorsOfSprings(oldListOfDisps, newListOfDisps, listOfSprings, listOfMinDisps, listOfMaxDisps, oldListOfMinDisps, oldListOfMaxDisps, oldListOfPlasticDisps, newListOfPlasticDisps, listOfSpringLoadingStages, listOfUnlDisp, listOfRelDisp, oldListOfUnlDisp, oldListOfRelDisp, maxDispPerIter, minDispPerIter, unlDispPerIter, relDispPerIter, oldSpringStages, newSpringStages);
	// <End of setting up of variables useful during the analysis>

	//<Start of calculations that are not required to be performed every loadstep>
	Matrix Fconst = Load::AssembleLoadMatrixWithFlag(listOfNodes, listOfLoads, "constant"); //assembles the matrix of constant forces
	Matrix Finc = Matrix(1, 1);
	if (type == "seismic") {
		Finc = Load::AssembleLoadMatrixWithFlag(listOfNodes, listOfLoads, "seismic"); //assembles the matrix of seismic forces
	}
	else if (type == "impulse") {
		Finc = Load::AssembleLoadMatrixWithFlag(listOfNodes, listOfLoads, "impulse"); //assembles the matrix of impulse forces
	}
	else {
		Finc = Load::AssembleLoadMatrixWithFlag(listOfNodes, listOfLoads, "increment"); //assembles the matrix of incremental forces
	}
	Matrix FConstRed = Load::GetReducedLoadMatrix(Fconst, listOfSups); //reduced (only free DOFs) version of the matrix of constant forces
	Matrix FIncRed = Load::GetReducedLoadMatrix(Finc, listOfSups); //reduced (only free DOFs) version of the matrix of incremental force

	std::vector<std::vector<ShellElement>> shellElemVecs = Solver::SetUpThreadsForShells(listOfShells, std::thread::hardware_concurrency()); //separates the total shell elements amongst the available threads
	Matrix shellStiff = Solver::CompleteShellStiffMatrixThreads(listOfNodes, listOfShells, listOfSups, std::thread::hardware_concurrency(), shellElemVecs); //calcualtes the reduced stiffness amtrix only accounting for the shell elements
	Matrix m(shellStiff.GetDimX(), shellStiff.GetDimY()); //Initializes the 'total' reduced matrix
	Solver::CompleteSpringStiffMatrixThreadsDispBasedAfterShells(listOfSprings, m, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter); //this will update 'm' with the values of the spring elements
	Matrix modalTotalStiffMatrix = m + shellStiff; //complete reduced total stiffness matrix. Used for modal analysis purposes
	Matrix shellRestricStiff = Solver::CompleteShellRestrictedStiffMatrixThreads(listOfNodes, listOfShells, listOfSups, std::thread::hardware_concurrency(), shellElemVecs); //calculates the stiffness matrix of the restricted DOFs accounting only for the shell elements
	Matrix shellMass = Solver::CompleteShellMassMatrixThreads(listOfNodes, listOfShells, listOfSups, std::thread::hardware_concurrency(), shellElemVecs);
	Matrix totMass = Mass::AddExplicitMassesOnExistingMatrix(shellMass, listOfMasses, listOfSups); //adds any additional mass to the existing shell masses

	//Matrix modalTotalStiffMatrix(MatrixOperation::CopyMatrixDouble(totalStiff), shellStiff.GetDimX(), shellStiff.GetDimY()); //copy of the total stiffness matrix that will be modified in the modal analysis calculatiosn

	std::vector<double> natFreq; //vector that will store the values of the natural frequencies
	Matrix modeShapes; //matrix that will store the displacement values for the different modes
	std::vector<std::vector<int>> totalMassDOFVec = ShellElement::GetTotalGlobalMassDOFVector(listOfShells); //gets a vector of all the DOFS that DO NOT Have mass. It is used everytime the matrices are required to be reduced only to the DOFS that have mass
	Solver::CalculateNaturalFrequenciesAndModeShapes(modalTotalStiffMatrix, totMass, natFreq, modeShapes, totalMassDOFVec); //Calculates the nat freqs and mode shapes and store them on the natFreq vector and the modeShapes Matrix
	Matrix totalModes = Solver::GetTotalModalMatrix(modeShapes, listOfSups, listOfNodes, listOfShells, listOfMasses); //Creates a matrix similar to modeShapes, but will all DOFS, not only the DOFS with masses
	
	std::vector<double> rayleighConstants; //vector to store the rayleigh constants
	rayleighConstants.reserve(2);
	int mode1 = 1, mode2 = 2; //for now
	double damp1 = 0, damp2 = 0; //for now
	//rayleighConstants = Solver::RayleighDampingConstants(mode1, damp1, mode2, damp2, natFreq);
	
	//Wilson-Theta: Gama = 1/2, Beta = 1/6, Theta = 1.42, alpha = 0; HTT-alpha: Gama = 0.6, Beta 0.3025, theta = 1, alpha = 0.1
	TimeIntegrationMethod IntMethod(IntegrationMethod::AverageNewmark);

	//Perform the load steps
	for (int step = 0; step < totalSteps; step++) {
		std::cout << "" << std::endl;
		std::cout << "Loadstep " << step + 1 << std::endl;
		std::cout << "Iteration " << 1 << std::endl;

		double time = (step + 1) * deltaT;

		//<Start of matrix stiffness solving>
		//<Set up stiffness matrix>
		Matrix m(shellStiff.GetDimX(), shellStiff.GetDimY()); //Initializes the 'total' reduced matrix
		Solver::CompleteSpringStiffMatrixThreadsDispBasedAfterShells(listOfSprings, m, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter); //this will update 'm' with the values of the spring elements
		Matrix mRed = m + shellStiff; //complete reduced total stiffness matrix

		//<End Set up stiffness Matrix>

		//<Set up damping matrix>
		Matrix damp(mRed.GetDimX());
		if (damp1 != 0 || damp2 != 0) {
			Matrix damp = Solver::RayleighDampingMatrix(totMass, mRed, rayleighConstants);
		}
		//<End set up damping matrix>

		//<Set dynamic 'stiff' matrix>
		double const1 = (1 - IntMethod.GetAlphaF()) * IntMethod.GetNewmarkGama() * IntMethod.GetWilsonTheta()  * deltaT;
		double constM = 1 - IntMethod.GetAlphaM();
		double const2 = IntMethod.GetNewmarkBeta() * pow(IntMethod.GetWilsonTheta()  * deltaT, 2);
		Matrix mult2 = totMass * constM; //kept separate to be used in a subsequent iteration
		Matrix add = (mult2 + damp * const1); //kept separate so it can be used below
		Matrix kDyn = add * (1 / const2);

		Matrix totStiff = mRed * (1 - IntMethod.GetAlphaF()) + kDyn;

		double highStiff = 0;
		//<End set dynamic stiff matrix>

		//<Set up dynamic force>
		Matrix prevTerm = (prevDisp * (1 / pow(IntMethod.GetWilsonTheta() * deltaT, 2)) + prevVel * (1 / (IntMethod.GetWilsonTheta() * deltaT)) + prevAcc * 0.5) * (1 / IntMethod.GetNewmarkBeta());
		Matrix firstMatrix = add * prevTerm;

		Matrix add11 = (prevVel + prevAcc * (IntMethod.GetWilsonTheta()  * deltaT) * (1 - IntMethod.GetAlphaF())); //kept separate to be used below
		Matrix secondMatrix = damp * add11;

		Matrix thirdMatrix = totMass * prevAcc;

		Matrix load(1, 1);
		Matrix prevLoad(1, 1);
		if (type == "seismic") {
			load = (totMass * SeismicLoad::GetSeismicLoadVector(sLoad, FIncRed, time)) * (-1);
			prevLoad = (totMass * SeismicLoad::GetSeismicLoadVector(sLoad, FIncRed, time - deltaT)) * (-1);
		}
		else if (type == "impulse") {
			load = FIncRed * ImpulseLoad::LoadFromTime(impLoad, time);
			prevLoad = FIncRed * ImpulseLoad::LoadFromTime(impLoad, time - deltaT);
		}
		else {
			double magnitude = 33000;
			double period = 1;
			load = FIncRed * (Load::SampleDynamicForceFunction(magnitude, (1 / period) * 2 * 3.1415, 0, time));
			prevLoad = FIncRed * Load::SampleDynamicForceFunction(magnitude, (1 / period) * 2 * 3.1415, 0, time - deltaT);
		}

		Matrix FInc1 = (load * (IntMethod.GetWilsonTheta()  * (1 - IntMethod.GetAlphaF())));
		Matrix FInc2(FInc1.GetDimX(), 1);
		if (IntMethod.GetAlphaF() == 0) {
			FInc2 = (prevLoad * (1 - IntMethod.GetWilsonTheta()));
		}
		else {
			FInc2 = (prevLoad * (IntMethod.GetAlphaF()));
		}
		Matrix deltaF = FInc1 + FInc2;//kept separate to be used below
		Matrix FIncDyn = deltaF - mRed * IntMethod.GetAlphaF() * prevDisp;
		Matrix FDyn = FIncDyn + firstMatrix - secondMatrix - thirdMatrix;
		//<End of set up dynamic force>


		Matrix F_iter = FConstRed + FDyn; //adds the constant and incremental terms of the applied laods, considering the current loadstep
		Solver::PlasticDisplacementLoadForce(F_iter, listOfSups, newListOfPlasticDisps, listOfSprings, listOfMinDisps, listOfMaxDisps, newListOfDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter);
		Matrix d_iter = MatrixOperation::FullCholesky(totStiff, F_iter); //solving the F = Kd system using Cholesky
		//<End of matrix stiffness solving>

		//<Calculate aceleration and velocity>
		Matrix redPrevDisp = ShellElement::ReducedAccelerationForceMatrix(prevDisp, totalMassDOFVec);
		Matrix redVel = ShellElement::ReducedAccelerationForceMatrix(prevVel, totalMassDOFVec);
		Matrix redAcc = ShellElement::ReducedAccelerationForceMatrix(prevAcc, totalMassDOFVec);
		Matrix acc(1, 1);
		Matrix curVel(1, 1);
		Matrix completeD_iter;
		IterationManager::CalculateAccelerationAndVelocity(d_iter, totalMassDOFVec, redPrevDisp, redVel, redAcc, IntMethod, deltaT, prevAcc, acc, curVel, prevVel, listOfSups, listOfNodes, completeD_iter);
		//<End calcualtion of acc and vel>
		
		Displacement::UpdatePositionVectorsOfSprings(oldListOfDisps, newListOfDisps, completeD_iter, listOfSprings, listOfMinDisps, listOfMaxDisps, 
													 oldListOfMinDisps, oldListOfMaxDisps, oldListOfPlasticDisps, newListOfPlasticDisps, listOfSpringLoadingStages, 
													 listOfUnlDisp, listOfRelDisp, oldListOfUnlDisp, oldListOfRelDisp, maxDispPerIter, minDispPerIter, unlDispPerIter, relDispPerIter, oldSpringStages, newSpringStages); //This function is responsible for updating all the parameters relevant tot he identification of which region the spring is at in its cyclic material model

		std::vector<Node> newNodes = Displacement::GetNewNodalCoordinates(originalList, completeD_iter); //this add the increment completeD_iter to the list of nodes coord
		if (Spring3D::CheckMaterialNonlinearityConvergenceDispBased(listOfSprings, oldListOfDisps, newListOfDisps, listOfMaxDisps, listOfMinDisps, oldListOfMaxDisps, oldListOfMinDisps, oldListOfPlasticDisps, newListOfPlasticDisps, listOfSpringLoadingStages, 0.001, breakAnalysis, unlDispPerIter, relDispPerIter, oldListOfUnlDisp, oldListOfRelDisp, oldSpringStages, newSpringStages)) { //check if springs have converged
			//this means that all spring elements converged. Ready for next loadstep or iteration
			IterationManager::PostConvergenceProcedures(listOfSprings, listOfSpringLoadingStages, newSpringStages, shellRestricStiff, listOfSups, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps,
				unlDispPerIter, relDispPerIter, d_iter, listOfNodes, F_iter, 0, forcePerStep, velPerStep, accPerStep, curVel, acc, maxDispPerIter, minDispPerIter, listOfUnlDisp, listOfRelDisp, prevDisp, prevVel, prevAcc);
		}
		else if (breakAnalysis) { //if convergence test returns a 'breakAnalysis' flag.
			std::cout << "Analysis didn't converge at loadstep " << step + 1 << " and iteration " << 1 << std::endl;
			break; //this will exit this for loop
		}
		else { //if the convergence test failed. This means not all spring element converged. Need to do iterations in order to converge before going to next loadStep.
			bool converged = false; //flag that indicates if the analysis has converged
			//F_iter.DestroyMatrixDouble(); //this need to be destroyed before the next iteration because the previous matrix have different displacement factors.


			for (int i = 0; i < nIter - 1; i++) { //Perform the iterations
				std::cout << "Iteration " << i + 2 << std::endl;

				// <Start of stiffness matrix solving of the displacements>
				Matrix m(shellStiff.GetDimX(), shellStiff.GetDimY());
				Solver::CompleteSpringStiffMatrixThreadsDispBasedAfterShells(listOfSprings, m, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter);
				Matrix mRed = m + shellStiff;
				
				//<Set up damping matrix>
				Matrix damp(mRed.GetDimX());
				if (damp1 != 0 || damp2 != 0) {
					Matrix damp = Solver::RayleighDampingMatrix(totMass, mRed, rayleighConstants);
				}
				//<End set up damping matrix>
				
				//<Set dynamic 'stiff' matrix>
				Matrix add = mult2 + damp * const1;
				Matrix kDyn = add * (1 / const2);
				Matrix totStiff = mRed * (1 - IntMethod.GetAlphaF()) + kDyn;
				double highStiff = 0;
				//<End set dynamic stiff matrix>

				//<Set up dynamic force>
				Matrix firstMatrix = add * prevTerm;
				Matrix secondMatrix = damp * add11;
				Matrix FIncDyn = deltaF - mRed * IntMethod.GetAlphaF() * prevDisp;
				Matrix FDyn = FIncDyn + firstMatrix - secondMatrix - thirdMatrix;
				//<End of set up dynamic force>

				Matrix F_iter = FConstRed + FDyn;
				Solver::PlasticDisplacementLoadForce(F_iter, listOfSups, newListOfPlasticDisps, listOfSprings, listOfMinDisps, listOfMaxDisps, newListOfDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter);
				Matrix dNew = MatrixOperation::FullCholesky(totStiff, F_iter);
				// <End of stiffness matrix solving of the displacements>

				//<Calculate aceleration and velocity>
				Matrix completeD_prime;
				IterationManager::CalculateAccelerationAndVelocity(dNew, totalMassDOFVec, redPrevDisp, redVel, redAcc, IntMethod, deltaT, prevAcc, acc, curVel, prevVel, listOfSups, listOfNodes, completeD_prime);

				//<End calcualtion of acc and vel>

				Displacement::UpdatePositionVectorsOfSprings(oldListOfDisps, newListOfDisps, completeD_prime, listOfSprings, listOfMinDisps, listOfMaxDisps, oldListOfMinDisps, oldListOfMaxDisps, oldListOfPlasticDisps, newListOfPlasticDisps, listOfSpringLoadingStages, listOfUnlDisp, listOfRelDisp, oldListOfUnlDisp, oldListOfRelDisp, maxDispPerIter, minDispPerIter, unlDispPerIter, relDispPerIter, oldSpringStages, newSpringStages);

				newNodes = Displacement::GetNewNodalCoordinates(originalList, completeD_prime); //this add the increment completeD_iter to the list of nodes coord
				if (Spring3D::CheckMaterialNonlinearityConvergenceDispBased(listOfSprings, oldListOfDisps, newListOfDisps, listOfMaxDisps, listOfMinDisps, oldListOfMaxDisps, oldListOfMinDisps, oldListOfPlasticDisps, newListOfPlasticDisps, listOfSpringLoadingStages, 0.001, breakAnalysis, unlDispPerIter, relDispPerIter, oldListOfUnlDisp, oldListOfRelDisp, oldSpringStages, newSpringStages)) {

					IterationManager::PostConvergenceProcedures(listOfSprings, listOfSpringLoadingStages, newSpringStages, shellRestricStiff, listOfSups, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps,
						unlDispPerIter, relDispPerIter, dNew, listOfNodes, F_iter, 0, forcePerStep, velPerStep, accPerStep, curVel, acc, maxDispPerIter, minDispPerIter, listOfUnlDisp, listOfRelDisp, prevDisp, prevVel, prevAcc);

					converged = true;
					break; //this will exit the for loop
				}

				if (breakAnalysis) {
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
				IterationManager::PostConvergenceProcedures(listOfSprings, listOfSpringLoadingStages, newSpringStages, shellRestricStiff, listOfSups, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps,
					unlDispPerIter, relDispPerIter, redDisplacement, listOfNodes, forceMatrix, 0, forcePerStep, velPerStep, accPerStep, curVel, acc, maxDispPerIter, minDispPerIter, listOfUnlDisp, listOfRelDisp, prevDisp, prevVel, prevAcc);
			}
		}

		if (breakAnalysis) {
			break;  //this will exit this for loop
		}

		nodesPerStep.emplace_back(newNodes);
	}

	FileOperation::SaveResultsFile(fileName, nodesPerStep, originalList, forcePerStep, natFreq, totalModes);
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


void IterationManager::PostConvergenceProcedures(std::vector<Spring3D> &listOfSprings, std::vector<std::vector<std::string>> &listOfSpringLoadingStages, std::vector<std::vector<std::string>> &newSpringStages,
	Matrix &shellRestricStiff, std::vector<Support> &listOfSups, std::vector<std::vector<double>> &newListOfDisps, std::vector<std::vector<double>> &listOfMinDisps, std::vector<std::vector<double>> &listOfMaxDisps,
	std::vector<std::vector<double>> &newListOfPlasticDisps, std::vector<std::vector<double>> &unlDispPerIter, std::vector<std::vector<double>> &relDispPerIter, Matrix &dNew, std::vector<Node> &listOfNodes,
	Matrix& F_iter, double highStiff, std::vector<Matrix> &forcePerStep, std::vector<Matrix> &velPerStep, std::vector<Matrix> &accPerStep, Matrix &curVel, Matrix &acc, std::vector<std::vector<double>> &maxDispPerIter, 
	std::vector<std::vector<double>> &minDispPerIter, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp, Matrix &prevDisp, Matrix &prevVel, Matrix &prevAcc) {
	
	
	Spring3D::UpdateSpringLoadStages(listOfSprings, listOfSpringLoadingStages, newSpringStages);

	//<Start of solving the restricted stiffness equations to obtain reactions>
	Matrix m2(shellRestricStiff.GetDimX(), shellRestricStiff.GetDimY());
	Solver::CompleteSpringRestrictedStiffMatrixThreadsDispBasedAfterShells(listOfSprings, m2, listOfSups, newListOfDisps, listOfMinDisps, listOfMaxDisps, newListOfPlasticDisps, listOfSpringLoadingStages, newSpringStages, unlDispPerIter, relDispPerIter);
	Matrix mRest = shellRestricStiff + m2;
	Matrix totalDisplacement = Displacement::GetTotalDisplacementNotOrganized(dNew, listOfSups, listOfNodes, listOfSprings, newListOfPlasticDisps);
	Matrix reactions = mRest * totalDisplacement;
	//<End of solving the restricted stiffness equations to obtain reactions>

	Matrix totalForce = Load::GetTotalForceNotOrganized(F_iter, reactions, listOfSups, listOfNodes); //"Not organized" means that the order on the matrix is not the correct DOF order
	Matrix completeForce = Load::GetTotalForceMatrix(totalForce, listOfSups, listOfNodes, highStiff, 0); //Now organized
	forcePerStep.emplace_back(MatrixOperation::CopyMatrixDouble(completeForce), completeForce.GetDimX(), completeForce.GetDimY()); //storing the total forces on this converged loadstep in the appropriate array
	velPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(curVel), curVel.GetDimX(), curVel.GetDimY());
	accPerStep.emplace_back(MatrixOperation::CopyMatrixDouble(acc), acc.GetDimX(), acc.GetDimY());

	//<Start of updating important lists for cyclic localization>
	maxDispPerIter = listOfMaxDisps;
	minDispPerIter = listOfMinDisps;
	unlDispPerIter = listOfUnlDisp;
	relDispPerIter = listOfRelDisp;
	//<End of updating important lists for cyclic localization>

	prevDisp.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(dNew));
	prevVel.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(curVel));
	prevAcc.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(acc));
}

void IterationManager::CalculateAccelerationAndVelocity(Matrix &dNew, std::vector<std::vector<int>> &totalMassDOFVec, Matrix &redPrevDisp, Matrix &redVel,
	Matrix &redAcc, TimeIntegrationMethod& IntMethod, double& deltaT, Matrix &prevAcc, Matrix &acc, Matrix &curVel, Matrix &prevVel, std::vector<Support> &listOfSups,
	std::vector<Node> &listOfNodes, Matrix &completeD_prime){
	Matrix redDisp = ShellElement::ReducedAccelerationForceMatrix(dNew, totalMassDOFVec);
	Matrix deltaAccWilson = (redDisp - redPrevDisp) * (1 / (IntMethod.GetNewmarkBeta() * pow(IntMethod.GetWilsonTheta()  * deltaT, 2))) - redVel * (1 / (IntMethod.GetNewmarkBeta() * IntMethod.GetWilsonTheta()  * deltaT)) - redAcc * (1 / (2 * IntMethod.GetNewmarkBeta()));
	Matrix accRed = deltaAccWilson * (1 / IntMethod.GetWilsonTheta()) + redAcc;
	acc = ShellElement::ConvertAccFromReducedToTotal(accRed, totalMassDOFVec, prevAcc.GetDimX());

	Matrix velRed = redAcc * deltaT + (accRed - redAcc) * IntMethod.GetNewmarkGama() * deltaT + redVel;
	//Matrix velRed = (redDisp - redPrevDisp) * (NewMarkGama / (NewMarkBeta * deltaT)) - redVel * (NewMarkGama / NewMarkBeta) - redAcc * ((deltaT * NewMarkGama / (NewMarkBeta * 2)) - deltaT) + redVel;
	curVel = ShellElement::ConvertAccFromReducedToTotal(velRed, totalMassDOFVec, prevVel.GetDimX());
	Matrix redNewDisp = redVel * deltaT + redAcc * (pow(deltaT, 2) / 2) + (accRed - redAcc) * IntMethod.GetNewmarkBeta() * pow(deltaT, 2) + redPrevDisp;
	Matrix disp = ShellElement::ConvertAccFromReducedToTotal(redNewDisp, totalMassDOFVec, dNew.GetDimX());
	dNew = MatrixOperation::FillMatrixBasedOnOtherMatrix(disp, dNew);
	completeD_prime = Displacement::GetTotalDisplacementMatrix(dNew, listOfSups, listOfNodes); //displacements of all nodes
}

IterationManager::~IterationManager()
{
}
