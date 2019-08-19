#include "pch.h"
#include <iostream>
#include <mutex>

Spring3D::Spring3D()
{
}

//<summary>Initializes a new Spring3D element</summary>
//<ID>The ID of the element</ID>
//<listOfMaterials>List of materials used in the element. One for each direction</listOfMaterials>
//<axialDir>Defines the direction of the axial direction of the spring. Possible values='x', 'y', and 'z'</axialDir>
//<shearDir>Defines the direction of the shear direction of the spring. Possible values='x', 'y', and 'z'</shearDir>
Spring3D::Spring3D(int ID, Node* n1, Node* n2, std::vector<SpringMaterialModels*> listOfMaterials, char axialDir, char shearDir) {
	_ID = ID;
	_n1 = n1;
	_n2 = n2;
	_listOfMaterials = listOfMaterials;
	_axialDir = axialDir;
	_shearDir = shearDir;
	_outPlaneDir = GetOutOfPlaneDirection();
	_axialUnitVec = GetAxialUnitVector();
	_shearUnitVec = GetShearUnitVector();
}

std::string Spring3D::ToString() const {
	std::string spring = "";
	spring += "(";
	spring += std::to_string(_ID);
	spring += ")";
	spring += "(";
	spring += "Node 1: ";
	spring += std::to_string(*_n1->GetID());
	spring += ", ";
	spring += "Node 2: ";
	spring += std::to_string(*_n2->GetID());
	spring += ", ";
	spring += "Axial Dir: ";
	spring += _axialDir;
	spring += ", ";
	spring += "Shear Dir: ";
	spring += _shearDir;
	spring += ", ";
	spring += "Out of Plane Dir: ";
	spring += _outPlaneDir;
	spring += ")";
	return spring;
}

int Spring3D::GetID() const
{
	return _ID;
}

//<summary>Returns the first node of the element</summary>
Node Spring3D::GetNode1() {
	return *_n1;
}

//<summary>Returns the second node of the element</summary>
Node Spring3D::GetNode2() {
	return *_n2;
}

//<summary>Returns the list of materials of the element</summary>
std::vector<SpringMaterialModels*> Spring3D::GetListOfMaterials() {
	return _listOfMaterials;
}

//<summary>Returns the list with the global direction of each material direction of the element</summary>
std::vector<int> Spring3D::GetListOfGlobalMaterialDirections() {
	std::vector<char> dirVec;
	dirVec.reserve(3);
	dirVec.emplace_back(_axialDir);
	dirVec.emplace_back(_shearDir);
	dirVec.emplace_back(_outPlaneDir);

	char auxList[3] = { 'x', 'y', 'z' };
	std::vector<int> pos;

	for (int i = 0; i < dirVec.size(); i++) {
		std::vector<char>::iterator vec = std::find(dirVec.begin(), dirVec.end(), auxList[i]);
		pos.push_back(vec - dirVec.begin());
	}
	return pos;
}

//<summary>Returns name of the out of plane direction of the spring based on the axial and shear directions specified</summary>
char Spring3D::GetOutOfPlaneDirection() {
	if ((_axialDir == 'x' && _shearDir == 'y') || (_axialDir == 'y' && _shearDir == 'x')) {
		return 'z';
	}
	else if ((_axialDir == 'x' && _shearDir == 'z') || (_axialDir == 'z' && _shearDir == 'x')) {
		return 'y';
	}
	else {
		return 'x';
	}

}

//<summary>Returns the unit vector of the axial direction</summary>
Matrix Spring3D::GetAxialUnitVector() {
	Matrix unitDirAxial;
	if (_axialDir == 'x') {
		unitDirAxial = VectorOperation::UnitVecX();
	}
	else if (_axialDir == 'y') {
		unitDirAxial = VectorOperation::UnitVecY();
	}
	else {
		unitDirAxial = VectorOperation::UnitVecZ();
	}
	return unitDirAxial;
}

//<summary>Returns the unit vector of the shear direction</summary>
Matrix Spring3D::GetShearUnitVector() {
	Matrix unitDirShear;
	if (_shearDir == 'x') {
		unitDirShear = VectorOperation::UnitVecX();
	}
	else if (_shearDir == 'y') {
		unitDirShear = VectorOperation::UnitVecY();
	}
	else {
		unitDirShear = VectorOperation::UnitVecZ();
	}
	return unitDirShear;
}

//<summary>Returns the unit vector of the out of plane direction</summary>
Matrix Spring3D::GetUnitZVector() {

	Matrix mult = VectorOperation::VectorCrossProduct(_axialUnitVec, _shearUnitVec);
	Matrix unit = VectorOperation::UnitVector(mult);

	return unit;
}

//<summary>Returns the vector between the two defining nodes of the element</summary>
Matrix Spring3D::GetVectorBetweenNodes() {
	Matrix vec = VectorOperation::VectorFromNodes(*_n2, *_n1);
	return vec;
}

//<summary>Returns the vector with the global DOF coordinates of the elements of the stiffness matrix of the element</summary>
//<comment>Not used anymore</comment>
std::vector<double> Spring3D::GlobalDOFVectorNotUsed() {
	std::vector<double> vec;
	vec.reserve(6);
	int DOF = 6; //dictated by the DOF of the shell element
	int init = (*_n1->GetID() - 1)*DOF;//if node 1, its 'coordinate' in the stiffness matrix is 0, if 2 -> 5, if 3 -> 10;
	vec.emplace_back(init);
	vec.emplace_back(init + 1);
	vec.emplace_back(init + 2);
	int init2 = (*_n2->GetID() - 1)*DOF;//if node 1, its 'coordinate' in the stiffness matrix is 0, if 2 -> 5, if 3 -> 10;
	vec.emplace_back(init2);
	vec.emplace_back(init2 + 1);
	vec.emplace_back(init2 + 2);
	return vec; // every 3 terms in the vector is related to a node coordinates in the global stiff matrix
}

//<summary>Returns the vector with the global DOF coordinates of the elements of the stiffness matrix of the element</summary>
//<supList>The list of supports in the system</supList>
//<comment>This list already accounts for the DOFs removed due to the supports</comment>
void Spring3D::CalculateGlobalDOFVector(const std::map<int, Support*>* supList, const int* DOF) {
	std::vector<std::vector<int>> vecFree;
	std::vector<std::vector<int>> vecRestr;
	std::vector<Node> nodeVec;
	nodeVec.push_back(*_n1);
	nodeVec.push_back(*_n2);

	int DOFLocal = 3;
	for (int i = 0; i < 2; i++) // for each node
	{
		Node n = nodeVec[i];
		//int count = Support::NumberOfDOFBeforeNode(n.GetID(), supList);
		int init = (*n.GetID() - 1)* *DOF;//this is the location in the Global DOF. If node 1, its 'coordinate' in the stiffness matrix is 0, if 2 -> 5, if 3 -> 10;
		int initLocal = i * DOFLocal; //this is the location in the element stiffness matrix
		
		for (int j = 0; j < DOFLocal; j++) {
			int val = init + j;
			if (!Support::IsDOFConstrained(&val, supList, DOF)) { //if true = DOF of the node is not constrained
				int count = Support::NumberOfDOFBeforeDOF(&val, supList, DOF);
				std::vector<int> vec1 = { init + j - count, initLocal + j };
				vecFree.push_back(vec1);
			}
			else { // DOF of the node is constrained
				int count = Support::NumberOfDOFBeforeDOF(&val, supList, DOF);
				std::vector<int> vec1 = { count, initLocal + j };
				vecRestr.push_back(vec1);
			}
		}
	}
	_globalDOFList = vecFree;// every 5 terms in the vector is related to a node coordinates in the global stiff matrix
	_globalRestDOFList = vecRestr;
}

//<summary>Returns the vector between the two defining nodes of the element</summary>
std::vector<std::vector<int>> Spring3D::GetGlobalDOFVector() {
	return _globalDOFList;
}

std::vector<std::vector<int>> Spring3D::GetGlobalRestDOFVector() {
	return _globalRestDOFList;
}

//<summary>Calculates the local stiffness matrix based on the displacement based</summary>
Matrix Spring3D::GetLocalStifnessMatrixDispBased(std::vector<double*> listOfPos, std::vector<double*> listOfMinDisp, std::vector<double*> listOfMaxDisp, std::vector<double*> listOfPlasticDisp, std::vector<std::string*> listOfLoadStage, std::vector<std::string*> &listOfStage, std::vector<double*> listOfUnlDisp, std::vector<double*> listOfRelDisp) {
	//disp is 0 for initial stiffness
	double** stiff = Matrix::CreateMatrixDouble(3);
	std::vector<int> pos = GetListOfGlobalMaterialDirections();

	for (int i = 0; i < 3; i++) { //for each direction of the Spring3D (i.e., x, y, and z)
		if (*listOfPos[i] == 0){ //this is used to verify if the initial stiffness should be used or if the actual displacement should be used to calcualte the actual stiffness
			stiff[pos[i]][pos[i]] = (*_listOfMaterials[pos[i]]).GetInitialStiffness();
		}
		else {
			stiff[pos[i]][pos[i]] = (*_listOfMaterials[pos[i]]).GetSecantStiffnessFromDisplacement(*listOfPos[i], *listOfPlasticDisp[i], *listOfMaxDisp[i], *listOfMinDisp[i], *listOfLoadStage[i], *listOfStage[i], *listOfUnlDisp[i], *listOfRelDisp[i]);
		}
	}
	return Matrix(stiff, 3);
}

//<summary>Calculates the transformation matrix</summary>
Matrix Spring3D::GetTransformationMatrix() {
	
	Matrix vecZ = GetUnitZVector();
	Matrix step1 = MatrixOperation::AddMatrixRight(_axialUnitVec, _shearUnitVec); //these vectors do not change as the spring deforms... Correct?
	Matrix step2 = MatrixOperation::AddMatrixRight(step1, vecZ);

	Matrix fill(3, 3);
	Matrix step3 = fill - step2;

	double** totalTransform = Matrix::CreateMatrixDouble(3, 6);
	Matrix totalTransformMatrix(totalTransform, 3, 6);

	/*
	for (int i = 0; i < 2; i++) {
		Matrix::AddMatrixAtPosition(totalTransformMatrix, step3, 1, i * 3 + 1);
	}
	*/
	MatrixOperation::AddMatrixAtPosition(totalTransformMatrix, step2, 1, 1);
	MatrixOperation::AddMatrixAtPosition(totalTransformMatrix, step3, 1, 4);

	return totalTransformMatrix;
}

//<summary>Calculates the global stiffness matrix based on the displacement based</summary>
Matrix Spring3D::GetGlobalStiffMatrixDispBased(AnalysisSpringRecorder* springRecorder) {
	Matrix transf = GetTransformationMatrix();
	Matrix transTransp = MatrixOperation::Transpose(transf);
	Matrix localStiff = GetLocalStifnessMatrixDispBased((*springRecorder->GetNewDisp())[_ID-1],
		(*springRecorder->GetNewMinDisp())[_ID - 1],
		(*springRecorder->GetNewMaxDisp())[_ID - 1],
		(*springRecorder->GetNewPlasticDisp())[_ID - 1],
		(*springRecorder->GetListStages())[_ID - 1],
		(*springRecorder->GetNewStages())[_ID - 1],
		(*springRecorder->GetNewUnlDisp())[_ID - 1],
		(*springRecorder->GetNewRelDisp())[_ID - 1]);

	Matrix mult2 = transTransp * localStiff * transf;

	return mult2;
}

//<summary>Assembles the members of the global stiffness matrix of a collection of Spring3D elements on the total stiffness matrix</summary>
//<comment>Not used</comment>
/*
void Spring3D::AssembleSpringGlobalMatrixOnComplete(std::vector<Spring3D> &vecEle, Matrix& complete, std::vector<std::vector<double>> &listOfDisp, std::vector<std::vector<double>> &listOfMinDisp, std::vector<std::vector<double>> &listOfMaxDisp, std::vector<std::vector<double>> &listOfPlasticDisp, std::vector<std::vector<std::string>> &listOfLoadStage, std::vector <std::vector<std::string>> &listOfStage, std::vector<std::vector<double>> &listOfUnlDisp, std::vector<std::vector<double>> &listOfRelDisp) {
	int DOF = 3; //# of degrees of freedom in each spring element
	for (int k = 0; k < vecEle.size(); k++) { // for each element
		Spring3D ele = vecEle[k];
		Matrix global = ele.GetGlobalStiffMatrixDispBased(listOfDisp[k], listOfMinDisp[k], listOfMaxDisp[k], listOfPlasticDisp[k], listOfLoadStage[k], listOfStage[k], listOfUnlDisp[k], listOfRelDisp[k]); //6x6 matrix returns
		std::vector<double> vec = ele.GlobalDOFVectorNotUsed();
		int runs = 2 * DOF; //this will be 6 in this case
		for (int i = 0; i < runs; i++) { //for all the lines in the local stiffness matrix
			int index1 = (vec)[i];
			for (int j = 0; j < runs; j++) { //for all the columns
				int index2 = (vec)[j];
				complete.GetMatrixDouble()[index1][index2] += global.GetMatrixDouble()[i][j];
			}
		}
	}
}
*/

//<summary>Assembles the members of the global stiffness matrix of a collection of Spring3D elements on the total stiffness matrix.</summary>
//<comment>Uses the displacement based method</comment>
Matrix Spring3D::AssembleSpringGlobalMatrixOnReducedSizedMatrix(const int* redSize, const std::map<int, Spring3D*>* vecEle, AnalysisSpringRecorder* springRecorder) {
	Matrix m(*redSize, *redSize);

	std::map<int, Spring3D*>::const_iterator it = vecEle->begin();

	while (it != vecEle->end()) { //for each spring element
		Spring3D* ele = it->second;
		Matrix global = ele->GetGlobalStiffMatrixDispBased(springRecorder); //6x6 matrix returns
		std::vector<std::vector<int>> vec = ele->GetGlobalDOFVector();
		for (int i = 0; i < vec.size(); i++) { //for all the lines in the local stiffness matrix
			int index1 = vec[i][0];
			for (int j = 0; j < vec.size(); j++) { //for all the columns
				int index2 = vec[j][0];
				m.GetMatrixDouble()[index1][index2] += global.GetMatrixDouble()[vec[i][1]][vec[j][1]];
			}
		}
		it++;
	}
	return m;
}

//<summary>Assembles the members of the global stiffness matrix of a collection of Spring3D elements on the total stiffness matrix.</summary>
//<comment>Uses the displacement based method</comment>
Matrix Spring3D::AssembleSpringGlobalRestrictedMatrixOnComplete(const PreAnalysisSetUp* setUp, const std::map<int, Spring3D*>* vecEle, AnalysisSpringRecorder* springRecorder) {
	int sizeRow = *setUp->StiffMatrixSize() - *setUp->ReducedStiffMatrixSize();
	int sizeCol = *setUp->StiffMatrixSize();
	Matrix m(sizeRow, sizeCol);
	
	std::map<int, Spring3D*>::const_iterator it = vecEle->begin();

	while (it != vecEle->end()) { //for each spring element
		Spring3D* ele = it->second;
		std::vector<std::vector<int>> vec = ele->GetGlobalRestDOFVector();
		if (vec.size() != 0) { //only calculate if there is some restricted DOF in this element
			std::vector<std::vector<int>> vec2 = ele->GetGlobalDOFVector();
			Matrix global = ele->GetGlobalStiffMatrixDispBased(springRecorder); //6x6 matrix returns
			int runs = vec.size(); //this will be 6 in this case
			for (int i = 0; i < runs; i++) { //for all the lines in the local stiffness matrix
				int index1 = vec[i][0];
				for (int j = 0; j < vec2.size(); j++) {
					int index2 = vec2[j][0];
					m.GetMatrixDouble()[index1][index2] += global.GetMatrixDouble()[vec[i][1]][vec2[j][1]];
				}
				for (int j = 0; j < runs; j++) { //for all the columns
					int index2 = vec[j][0] + *setUp->ReducedStiffMatrixSize();
					m.GetMatrixDouble()[index1][index2] += global.GetMatrixDouble()[vec[i][1]][vec[j][1]];
				}
			}
		}
		it++;
	}
	return m;
}

//<summary>Checks the material nonlinearity of Spring3D elements to decide if convergence has been achieved</summary>
//<comment>Based on the displacement based theory</comment>
bool Spring3D::CheckMaterialNonlinearityConvergenceDispBased(const StructureManager* structManager, const PreAnalysisSetUp* setUp, AnalysisSpringRecorder* springRecord, bool* breakAnalysis) {
	
	std::map<int, Spring3D*>::const_iterator it = structManager->SpringElements()->begin();

	while (it != structManager->SpringElements()->end()) {//for each Spring3D element
		Spring3D* spring = it->second;
		int springID = spring->GetID();
		std::vector<int> listOfPos = spring->GetListOfGlobalMaterialDirections();

		std::vector<SpringMaterialModels*> matList = spring->GetListOfMaterials();

		double oldStiff[3];
		double newStiff[3];

		for (int j = 0; j < matList.size(); j++) { //x then y then z direction
			oldStiff[j] = matList[listOfPos[j]]->GetSecantStiffnessFromDisplacement(*(*springRecord->GetOldDisp())[springID -1][j],
				* (*springRecord->GetOldPlasticDisp())[springID - 1][j], * (*springRecord->GetOldMaxDisp())[springID - 1][j], *(*springRecord->GetOldMinDisp())[springID - 1][j],
				* (*springRecord->GetListStages())[springID - 1][j], * (*springRecord->GetOldStages())[springID - 1][j],
				* (*springRecord->GetOldUnlDisp())[springID - 1][j], * (*springRecord->GetOldRelDisp())[springID - 1][j]); //should this last 2 be "OLD"?
			
			newStiff[j] = matList[listOfPos[j]]->GetSecantStiffnessFromDisplacement(*(*springRecord->GetNewDisp())[springID - 1][j],
				*(*springRecord->GetNewPlasticDisp())[springID - 1][j], *(*springRecord->GetNewMaxDisp())[springID - 1][j], *(*springRecord->GetNewMinDisp())[springID - 1][j],
				* (*springRecord->GetListStages())[springID - 1][j], * (*springRecord->GetNewStages())[springID - 1][j],
				* (*springRecord->GetNewUnlDisp())[springID - 1][j], * (*springRecord->GetNewRelDisp())[springID - 1][j]);
		}

		if (newStiff[0] == -1 || newStiff[1] == -1 || newStiff[2] == -1) {
			*breakAnalysis = true;
			return false;
		}

		for (int i = 0; i < matList.size(); i++) { //for each material on the element
			double convergence = abs((newStiff[i] - oldStiff[i]) / newStiff[i]);
			if (convergence > *setUp->Analysis()->ConvergenceLimit()) {
				return false;
			}
		}
		it++;
	}
	//if the next lines are executed, it means that no spring element had a convergence above the limit, therefore:
	return true;
}

std::vector<int> Spring3D::GetPlasticDispIndexes(const std::vector<std::vector<double*>> listOfPlasticDisp, const StructureManager* structManager, const int* DOF)
{
	std::vector<int> vec;
	
	for (int i = 0; i < listOfPlasticDisp.size(); i++) {//for each spring3D
		for (int j = 0; j < listOfPlasticDisp[i].size(); j++) { //for x, y, and z
			if (*listOfPlasticDisp[i][j] != 0) {
				int nodeID = *structManager->SpringElements()->find(i + 1)->second->GetNode2().GetID();
				int index = (nodeID - 1) * (*DOF) + j;
				int count = Support::NumberOfDOFBeforeDOF(&index, structManager->Supports(), DOF);
				vec.push_back(index - count);
			}
		}
	}
	return vec;
}

//<summary>Get a vector with the displacements of the Spring3D element from the total vector of displacements</summary>
Matrix Spring3D::GetElementGlobalDisplacementVector(Matrix completeD) {
	Matrix disp1 = Displacement::GetDisplacementByNodeID(*_n1->GetID(), completeD);
	Matrix disp2 = Displacement::GetDisplacementByNodeID(*_n2->GetID(), completeD);

	Matrix ans(6, 1);
	ans.GetMatrixDouble()[0][0] = disp1.GetMatrixDouble()[0][0];
	ans.GetMatrixDouble()[1][0] = disp1.GetMatrixDouble()[1][0];
	ans.GetMatrixDouble()[2][0] = disp1.GetMatrixDouble()[2][0];
	ans.GetMatrixDouble()[3][0] = disp2.GetMatrixDouble()[0][0];
	ans.GetMatrixDouble()[4][0] = disp2.GetMatrixDouble()[1][0];
	ans.GetMatrixDouble()[5][0] = disp2.GetMatrixDouble()[2][0];

	return ans;
}

void Spring3D::UpdateSpringLoadStages(const std::map<int, Spring3D*>* vecEle, AnalysisSpringRecorder* springRecord) {
	
	std::map<int, Spring3D*>::const_iterator it = vecEle->begin();

	while (it != vecEle->end()) {
		Spring3D* spring = it->second;
		int ID = spring->GetID();
		std::vector<int> pos = spring->GetListOfGlobalMaterialDirections();

		for (int j = 0; j < 3; j++) {
			(*springRecord->GetListStages())[ID - 1][j] = (*springRecord->GetNewStages())[ID - 1][j];
		}
		it++;
	}
}

Spring3D::~Spring3D()
{
}
