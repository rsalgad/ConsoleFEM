#include "pch.h"
#include "ShellElement.h"
#include "ElasticMaterial.h"
#include "OrthotropicElasticMaterial.h"
#include "Support.h"
#include "Mass.h"
#include "MatrixOperation.h"
#include "VectorOperation.h"
#include <math.h>
#include <iostream>
#include <thread>
#include <mutex>

ShellElement::ShellElement()
{
	/*
    (3)   (6)   (2)
	 4_____7_____3
	 |           |
(7) 8|     9     |6 (5)   <-- Node Configuration
	 |    (8)    |
	 |___________|
	 1     5     2   <- Normal Numbering
    (0)   (4)   (1)  <- Array Numbering
	*/
}

//<summary>Initializes a new shell element</summary>
//<ID>Initializes a new shell element</ID>
//<thick>The thickness of the element</thick>
//<layers>The number of layers of the element</layers>
//<matList>The list of materials (one for each layer)</matList>
ShellElement::ShellElement(int ID, Node* n1, Node* n2, Node* n3, Node* n4, Node* n5, Node* n6, Node* n7, Node* n8, Node* n9, double thick, int layers, std::vector<OrthotropicElasticMaterial> matList) {
	_ID = ID;
	//nodes should be numbered in sequence first with the outer nodes then the mid-nodes
	_nodeList.reserve(9);
	_nodeList.emplace_back(n1);
	_nodeList.emplace_back(n2);
	_nodeList.emplace_back(n3);
	_nodeList.emplace_back(n4);
	_nodeList.emplace_back(n5);
	_nodeList.emplace_back(n6);
	_nodeList.emplace_back(n7);
	_nodeList.emplace_back(n8);
	_nodeList.emplace_back(n9);
	_thick = thick;
	_layers = layers;
	_matList = matList;
	_v3 = Get_v3(true, 0, 0);
	_v1 = Get_v1(_v3);
	_v2 = Get_v2(_v3, _v1);
}

std::string ShellElement::ToString() const {
	std::string shell = "";
	shell += "(";
	shell += std::to_string(_ID);
	shell += ")";
	shell += "(";
	shell += "Layers: ";
	shell += std::to_string(_layers);
	shell += ", ";
	shell += "Thickness: ";
	shell += std::to_string(_thick);
	shell += ", ";
	shell += "Nodes: ";
	for (int i = 0; i < _nodeList.size(); i++){
		shell += std::to_string(_nodeList[i]->GetID());
		if (i != _nodeList.size() - 1) {
			shell += " ";
		}
	}
	shell += ")";
	return shell;
}

//<summary>Returns the ID of the element</summary>
int ShellElement::GetID() const {
	return _ID;
}

double ShellElement::GetMassRho() const {
	return _massRho;
}

//<summary>Returns the global degree of freedom vector of the element</summary>
const std::vector<std::vector<int>> ShellElement::GetGlobalDOFVector() const {
	return _globalDOFList;
}

//<summary>Returns the global degree of freedom vector of the element</summary>
const std::vector<std::vector<int>> ShellElement::GetGlobalMassDOFVector() const {
	return _globalMassDOFList;
}

const std::vector<std::vector<int>> ShellElement::GetGlobalRestDOFVector() const{
	return _globalRestDOFList;
}

//<summary>Calculates the global degree of freedom vector of the element</summary>
//<comment>This is a vector that stores the global DOF position of each term in the local stiffness matrix. This function is deprecated. </comment>
/*
std::vector<int> ShellElement::GlobalDOFVector() {
	std::vector<int> vec;
	vec.reserve(6 * 9);
	int DOF = 6;
	for (int i = 0; i < _nodeList.size(); i++) // for each node
	{
		Node n = *_nodeList[i];
		int init = (n.GetID() - 1)*DOF;//if node 1, its 'coordinate' in the stiffness matrix is 0, if 2 -> 5, if 3 -> 10;
		vec.emplace_back(init);
		vec.emplace_back(init + 1);
		vec.emplace_back(init + 2);
		vec.emplace_back(init + 3);
		vec.emplace_back(init + 4);
		vec.emplace_back(init + 5);
	}
	return vec; // every 5 terms in the vector is related to a node coordinates in the global stiff matrix
}
*/

//<summary>Calculates the global degree of freedom vector of the element</summary>
//<comment>This is a vector that stores the global DOF position of each term in the local stiffness matrix. It returns the DOFs that are not constrained by a support. </comment>
void ShellElement::CalculateGlobalDOFVector(std::vector<Support*> supList) {
	std::vector<std::vector<int>> vecFree;
	std::vector<std::vector<int>> vecRestr;

	int DOF = 6;
	for (int i = 0; i < _nodeList.size(); i++) // for each node
	{
		Node n = *_nodeList[i];
		//int count = Support::NumberOfDOFBeforeNode(n.GetID(), supList);
		int init = (n.GetID() - 1)*DOF;//if node 1, its 'coordinate' in the stiffness matrix is 0, if 2 -> 6, if 3 -> 12;
		int initLocal = i * DOF;

		if (Support::IsNodeConstrained(supList, n.GetID())) { //if true = node has some constraint
			for (int j = 0; j < DOF; j++) {
				if (!Support::IsDOFConstrained(init + j, supList)) { //if true = DOF of the node is not constrained
					int count = Support::NumberOfDOFBeforeDOF(init + j, supList);
					std::vector<int> vec1 = { init + j - count, initLocal + j};
					vecFree.push_back(vec1);
				}
				else { // DOF of the node is constrained
					int count = Support::NumberOfDOFBeforeDOF(init + j, supList);
					std::vector<int> vec1 = { count, initLocal + j };
					vecRestr.push_back(vec1);
				}
			}
		}
		else { //if node has no restraint, no need to check if 'IsDOFConstrained' all the time
			int count = Support::NumberOfDOFBeforeDOF(init, supList);
			for (int j = 0; j < DOF; j++) {
				std::vector<int> vec1 = { init + j - count, initLocal + j};
				vecFree.push_back(vec1);
			}
		}
	}
	_globalDOFList = vecFree; // every 6 terms in the vector is related to a node coordinates in the global stiff matrix
	_globalRestDOFList = vecRestr;
}


void ShellElement::CalculateGlobalMassDOFVector(std::vector<Mass*> massList, std::vector<Support*> supList) {
	std::vector<std::vector<int>> vecMass;

	int DOF = 6;
	for (int i = 0; i < _nodeList.size(); i++) // for each node
	{
		Node n = *_nodeList[i];
		//int count = Support::NumberOfDOFBeforeNode(n.GetID(), supList);
		int init = (n.GetID() - 1)*DOF;//if node 1, its 'coordinate' in the stiffness matrix is 0, if 2 -> 6, if 3 -> 12;
		int initLocal = i * DOF;

		if (Support::IsNodeConstrained(supList, n.GetID())) { //if true = node has some constraint
			for (int j = 0; j < DOF; j++) {
				if (!Support::IsDOFConstrained(init + j, supList)) { //if DOF is not constrained
					int count = Support::NumberOfDOFBeforeDOF(init + j, supList);
					int globDOF = init + j - count;
					int locDOF = initLocal + j;
					
					if (!Mass::HasDOFAppliedMass(massList, globDOF + count)) { //if DOF does not have mass
						bool rot = IsDOFRotational(locDOF);
						if (rot || _massRho == 0 || locDOF == 48 || locDOF == 49 || locDOF == 50) {//if it is not rotational nor the ninthNode
							std::vector<int> vec1 = { globDOF, locDOF };
							vecMass.push_back(vec1);
						}
					}
				}
			}
		}
		else { //if node has no restraint, no need to check if 'IsDOFConstrained' all the time
			int count = Support::NumberOfDOFBeforeDOF(init, supList);
			for (int j = 0; j < DOF; j++) {
				int globDOF = init + j - count;
				int locDOF = initLocal + j;

				if (!Mass::HasDOFAppliedMass(massList, globDOF + count)) { //if DOF does not have mass
					bool rot = IsDOFRotational(locDOF);
					if (rot || _massRho == 0 || locDOF == 48 || locDOF == 49 || locDOF == 50) {//if it is not rotational nor the ninthNode
						std::vector<int> vec1 = { globDOF, locDOF };
						vecMass.push_back(vec1);
					}
				}
			}
		}
	}
	_globalMassDOFList = vecMass; // every 6 terms in the vector is related to a node coordinates in the global stiff matrix
}

//<summary>Calculates the v1 vector used to defined the local axis of the element</summary>
//<v3>The v3 vectors that is required to defined the v1 vector </v3>
Matrix ShellElement::Get_v1(Matrix &v3) const {
	//creates a vector v1 which is used in the B Matrix

	double** v1 = Matrix::CreateMatrixDouble(3, 1);
	if (v3.GetMatrixDouble()[0][0] == 0 && v3.GetMatrixDouble()[2][0] == 0 && v3.GetMatrixDouble()[1][0] != 0) {
		v1[0][0] = -v3.GetMatrixDouble()[1][0];
	} else {
		v1[0][0] = v3.GetMatrixDouble()[2][0];
		v1[2][0] = -v3.GetMatrixDouble()[0][0];
	}
	Matrix V1(v1, 3, 1);
	Matrix ans = VectorOperation::UnitVector(V1);
	
	return ans;
}

//<summary>Calculates the v3 vector used to defined the local axis of the element</summary>
//<first>Defines if the v3 vector is being calculated at startup of the element, or by subsequent calculations </first>
//<e>The e position in natural coordinate that defines which point is being calculated </e>
//<n>The n position in natural coordinate that defines which point is being calculated </n>
Matrix ShellElement::Get_v3(bool first, double e, double n) const {
	Matrix unitV3(3, 1);
	if (!first){
		double rx = 0, ry = 0, rz = 0;
		for (int i = 0; i < _nodeList.size(); i++) { //should this include all 9 nodes or just the 8?
			double shape = ShapeFunctionN(i + 1, e, n);
			rx += shape * _nodeList[i]->GetRx();
			ry += shape * _nodeList[i]->GetRy();
			rz += shape * _nodeList[i]->GetRz();
		}
		if (rx == 0 && ry == 0 && rz == 0){
			unitV3.SetMatrixDouble(MatrixOperation::CopyMatrixDouble(_v3));
		}
		else {
			Matrix V3 = VectorOperation::RotateVectorIn3D(_v3, rx, ry, rz);
			unitV3 = VectorOperation::UnitVector(V3);
		}
	} else {
		Matrix vector1 = VectorOperation::VectorFromNodes(*_nodeList[0], *_nodeList[1]);
		Matrix vector2 = VectorOperation::VectorFromNodes(*_nodeList[0], *_nodeList[3]);
		Matrix V3 = VectorOperation::VectorCrossProduct(vector1, vector2);
		unitV3 = VectorOperation::UnitVector(V3);
	}
	return unitV3;
}

//<summary>Calculates the v2 vector used to defined the local axis of the element</summary>
//<v3>The v3 vectors that is required to defined the v2 vector </v3>
//<v1>The v3 vectors that is required to defined the v2 vector </v1>
Matrix ShellElement::Get_v2(Matrix &v3, Matrix &v1) const {
	Matrix V2 = VectorOperation::VectorCrossProduct(v3, v1);
	Matrix unitV2 = VectorOperation::UnitVector(V2);
	return unitV2;
}

//<summary>Calculates the 3x3 local-global axis transformation matrix following the formulas in the Hyrnyk PhD thesis</summary>
//<jacob>The jacobian matrix</jacob>
Matrix ShellElement::GetTransformationMatrix3x3Hyrnyk(const Matrix &jacob) const {
	double** J = jacob.GetMatrixDouble();
	double** vec1 = Matrix::CreateMatrixDouble(3, 1);
	vec1[0][0] = J[0][0];
	vec1[1][0] = J[0][1];
	vec1[2][0] = J[0][2];

	double** vec2 = Matrix::CreateMatrixDouble(3, 1);
	vec2[0][0] = J[1][0];
	vec2[1][0] = J[1][1];
	vec2[2][0] = J[1][2];

	double** vec3 = Matrix::CreateMatrixDouble(3, 1);
	/*
	vec3[0][0] = J[0][0];
	vec3[1][0] = J[1][1];
	vec3[2][0] = J[0][2];
	*/

	Matrix vec11(vec1, 3, 1);
	Matrix vec22(vec2, 3, 1);

	Matrix vecZ = VectorOperation::VectorCrossProduct(vec11, vec22);

	if (vecZ.GetMatrixDouble()[0][0] == 0 && vecZ.GetMatrixDouble()[2][0] == 0 && vecZ.GetMatrixDouble()[1][0] != 0) {
		vec3[0][0] = -vecZ.GetMatrixDouble()[1][0];
	}
	else {
		vec3[0][0] = vecZ.GetMatrixDouble()[2][0];
		vec3[2][0] = -vecZ.GetMatrixDouble()[0][0];
	}

	Matrix vecX(vec3, 3, 1);
	Matrix vecY = VectorOperation::VectorCrossProduct(vecZ, vecX);

	Matrix unitVecX = VectorOperation::UnitVector(vecX);
	Matrix unitVecY = VectorOperation::UnitVector(vecY);
	Matrix unitVecZ = VectorOperation::UnitVector(vecZ);

	Matrix matrix1 = MatrixOperation::AddMatrixRight(unitVecX, unitVecY);
	Matrix finalMatrix = MatrixOperation::AddMatrixRight(matrix1, unitVecZ);

	return finalMatrix;
}

//<summary>Calculates the 3x3 local-global axis transformation matrix following the formulas in the Hinton and Owen, 1984 book</summary>
//<jacob>The jacobian matrix</jacob>
//<comment>This function results in error. That might be some mistake with the formulations given by the book</comment>
Matrix ShellElement::GetTransformationMatrix3x3Book(const Matrix &jacob) const {
	double** J = jacob.GetMatrixDouble();

	double** vec1 = Matrix::CreateMatrixDouble(3, 1);
	vec1[0][0] = J[0][0];
	vec1[1][0] = J[1][1];
	vec1[2][0] = J[0][2];

	Matrix vecX(vec1, 3, 1);

	double** vec31 = Matrix::CreateMatrixDouble(3, 1);
	vec31[0][0] = J[0][0];
	vec31[1][0] = J[0][1];
	vec31[2][0] = J[0][2];

	double** vec32 = Matrix::CreateMatrixDouble(3, 1);
	vec32[0][0] = J[1][0];
	vec32[1][0] = J[1][1];
	vec32[2][0] = J[1][2];

	Matrix vec311(vec31, 3, 1);
	Matrix vec322(vec32, 3, 1);

	Matrix vecZ = VectorOperation::VectorCrossProduct(vec311, vec322);

	Matrix vecY = VectorOperation::VectorCrossProduct(vecZ, vecX);

	Matrix vecXUnit = VectorOperation::UnitVector(vecX);
	Matrix vecYUnit = VectorOperation::UnitVector(vecY);
	Matrix vecZUnit = VectorOperation::UnitVector(vecZ);

	Matrix matrix1 = MatrixOperation::AddMatrixRight(vecXUnit, vecYUnit);
	Matrix finalMatrix = MatrixOperation::AddMatrixRight(matrix1, vecZUnit);

	return finalMatrix;
}

//<summary>Calculates the 54x54 local-global axis transformation matrix following the formulas in the Hinton and Owen, 1984 book</summary>
//<jacob>The jacobian matrix</jacob>
Matrix ShellElement::GetTransformationMatrix(const Matrix& jacob) const {
	Matrix transform = GetTransformationMatrix3x3Hyrnyk(jacob); //For some reason the formulation contained in the book does not work

	double** totalTransform = Matrix::CreateMatrixDouble(54);
	Matrix totalTransformMatrix(totalTransform, 54, 54);
	for (int i = 0; i < 18; i++) { //each node has 2 instances of the transform matrix (1 for transl and 1 for rot). So 9 nodes require 18 instances.
		//Add the transform matrix at all the position in the 54x54 size of the transform matrix
		MatrixOperation::AddMatrixAtPosition(totalTransformMatrix, transform, i * 3 + 1, i * 3 + 1);
	}

	//For 5 DOF
	/*
	for (int i = 1; i <= 9; i++) {
		//Remove the rows and columns related to the degree of freedom that is not present in this element (i.e., rotation in z)
		totalTransformMatrix = Matrix::DeleteRowAndColumn(totalTransformMatrix, 5 * i);
	}
	*/
	return totalTransformMatrix;
}

//<summary>Calculates the result of the element's shape functions</summary>
//<i>Which node's shape function is being considered</i>
//<e>The e location in natural coordinates</e>
//<n>The n location in natural coordinates</n>
double ShellElement::ShapeFunctionN(int i, double &e, double &n) {
	
		switch (i) {
		case 1:
			return (1.0 - e)*(1.0 - n)*(-e - n - 1.0) / 4.0;
			break;
		case 2:
			return (1.0 + e)*(1.0 - n)*(e - n - 1.0) / 4.0;
			break;
		case 3:
			return (1.0 + e)*(1.0 + n)*(e + n - 1.0) / 4.0;
			break;
		case 4:
			return (1.0 - e)*(1.0 + n)*(-e + n - 1.0) / 4.0;
			break;
		case 5:
			return (1.0 + e)*(1.0 - e)*(1.0 - n) / 2.0;
			break;
		case 6:
			return (1.0 + e)*(1.0 + n)*(1.0 - n) / 2.0;
			break;
		case 7:
			return (1.0 + e)*(1.0 - e)*(1.0 + n) / 2.0;
			break;
		case 8:
			return (1.0 - e)*(1.0 + n)*(1.0 - n) / 2.0;
			break;
		case 9:
			return (1.0 - pow(e, 2.0))*(1 - pow(n, 2.0));
			break;
		default:
			//i needs to be between 1 and 8
			return 0; //If this is executed, i is either greater than 8 or less than 1. This is wrong!
		}
}

//<summary>Calculates the result of the shape functions specifically used to convert the 2x2 B matrix to 3x3 gauss points</summary>
//<i>Which node's shape function is being considered</i>
//<e>The e location in natural coordinates</e>
//<n>The n location in natural coordinates</n>
double ShellElement::ShapeFunctionNr(int i, double &e, double &n) {
	switch (i) {
	case 1:
		return (1.0 + e * (-sqrt(3)))*(1.0 + n * (-sqrt(3))) / 4.0;
		break;
	case 2:
		return (1.0 + e * (-sqrt(3)))*(1.0 + n * (sqrt(3))) / 4.0;
		break;
	case 3:
		return (1.0 + e * (sqrt(3)))*(1.0 + n * (-sqrt(3))) / 4.0;
		break;
	default:
		return (1.0 + e * (sqrt(3)))*(1.0 + n * (sqrt(3))) / 4.0;
	}
}

//<summary>Calculates the result of the element's shape functions derivatives</summary>
//<var>'e' = derivative in e; 'n' = derivative in n; 'c' = derivative in c</var>
//<i>Which node's shape function is being considered</i>
//<e>The e location in natural coordinates</e>
//<n>The n location in natural coordinates</n>
double ShellElement::ShapeFunctionNDerivative(int i, char var, double &e, double &n) {
	if (var == 'e') {
		switch (i)
		{
		case 1:
			return -(n - 1.0)*(2.0 * e + n) / (4.0);
			break;
		case 2:
			return -(n - 1.0)*(2.0 * e - n) / (4.0);
			break;
		case 3:
			return (n + 1.0)*(2.0 * e + n) / (4.0);
			break;
		case 4:
			return (n + 1.0)*(2.0 * e - n) / (4.0);
			break;
		case 5:
			return (n - 1.0)* e;
			break;
		case 6:
			return -(pow(n, 2.0) - 1.0) / 2.0;
			break;
		case 7:
			return (-n - 1.0)* e;
			break;
		case 8:
			return (pow(n,2.0) - 1.0) / 2.0;
			break;
		case 9:
			return (2.0*pow(n, 2) - 2)*e;
			break;
		default:
			return 0; //If this is executed, i is either greater than 8 or less than 1. This is wrong!
			break;
		}
	}
	else if (var == 'n') {
		switch (i)
		{
		case 1:
			return -(e - 1.0)*(2.0 * n + e) / (4.0);
			break;
		case 2:
			return (e + 1.0)*(2.0 * n - e) / (4.0);
			break;
		case 3:
			return (e + 1.0)*(2.0 * n + e) / (4.0);
			break;
		case 4:
			return -(e - 1.0)*(2.0 * n - e) / (4.0);
			break;
		case 5:
			return (pow(e, 2.0) - 1.0) / 2.0;
			break;
		case 6:
			return (-e - 1.0)* n;
			break;
		case 7:
			return -(pow(e, 2.0) - 1.0) / 2.0;
			break;
		case 8:
			return (e - 1.0)* n;
			break;
		case 9:
			return (2.0*pow(e, 2) - 2)*n;
			break;
		default:
			return 0; //If this is executed, i is either greater than 8 or less than 1. This is wrong!
			break;
		}
	}
	else if (var == 'c') {
		return 0; //The N's are not a function of z (or 'c'). Then the result is like deriving a constant, thus = 0;
	}
	else{
		return 0; //If this is executed, var is neither x or y. This is wrong!
	}
}

//<summary>Returns the width of the element</summary>
//<comment>Should only be used when the element is at its undeformed positions</comment>
double ShellElement::GetWidth() const {
	//Get the width based on the difference between nodes 1 and 2 x components.
	Matrix vec = VectorOperation::VectorFromNodes(*_nodeList[0], *_nodeList[1]); // assumes that nodes 1 and 2 are the edge nodes in the x direction
	double width = VectorOperation::VectorLength(vec);

	return width;
}

//<summary>Returns the height of the element</summary>
//<comment>Should only be used when the element is at its undeformed positions</comment>
double ShellElement::GetHeight() const {
	//Get the height based on the difference between nodes 2 and 3 y components.
	Matrix vec = VectorOperation::VectorFromNodes(*_nodeList[0], *_nodeList[3]); // assumes that nodes 1 and 4 are the edge nodes in the y direction
	double height = VectorOperation::VectorLength(vec);
	
	return height;
}

//<summary>Returns the thickness of the element</summary>
double ShellElement::GetTotalThickness() const {
	return _thick * _layers;
}

//<summary>Returns the number of layers of the element</summary>
double ShellElement::GetLayers() const {
	return _layers;
}

//<summary>Returns the thickness of each layer of the element</summary>
//<comment>Assumes uniform thickness</comment>
double ShellElement::GetLayerThickness() const {
	return GetTotalThickness() /_layers;
}

//<summary>Returns the list of nodes of the element</summary>
std::vector<Node *> ShellElement::GetNodeList() {
	return _nodeList;
}

//<summary>Defines the list of nodes of the element</summary>
void ShellElement::SetNodeList(std::vector<Node *> nodeList) {
	_nodeList = nodeList;
}

//<summary>Returns the complete (5x5) D matrix of the element</summary>
const Matrix ShellElement::GetDMatrix(OrthotropicElasticMaterial mat) const {
	double alpha = 5.0 / 6.0;
	double** D = Matrix::CreateMatrixDouble(5, 5);
	//double v21 = mat.GetStiffnessY() * mat.GetPoissonXY() / mat.GetStiffnessX();
	D[0][0] = mat.GetStiffnessX() / (1 - pow(mat.GetPoissonXY(), 2));
	D[0][1] = (mat.GetPoissonXY() * mat.GetStiffnessX()) / (1 - pow(mat.GetPoissonXY(), 2));
	D[1][0] = D[0][1];
	D[1][1] = mat.GetStiffnessY() / (1 - pow(mat.GetPoissonXY(), 2));
	D[2][2] = mat.GetShearStiffnessXY();
	D[3][3] = alpha * mat.GetShearStiffnessYZ();
	D[4][4] = alpha * mat.GetShearStiffnessXZ();

	Matrix Ds(D, 5);
	return Ds;
}

//<summary>Returns the shear (3x3) D matrix of the element</summary>
Matrix ShellElement::GetShearDMatrix(OrthotropicElasticMaterial mat) {
	double alpha = 5.0 / 6.0;
	double mult = mat.GetStiffnessX()*alpha / (2 * (1 + mat.GetPoissonXY()));
	double** D = Matrix::CreateMatrixDouble(3, 3);
	D[0][0] = 1 / alpha;
	D[1][1] = 1;
	D[2][2] = D[0][0];
	Matrix Ds(D, 3);
	Matrix ans = Ds * mult;
	return ans;
}

//<summary>Returns the membrane and flexure (2x2) D matrix of the element</summary>
Matrix ShellElement::GetMembFlexDMatrix(OrthotropicElasticMaterial mat) {
	double mult = mat.GetStiffnessX() / (1 - pow(mat.GetPoissonXY(), 2));
	double** D = Matrix::CreateMatrixDouble(2, 2);
	D[0][0] = 1;
	D[0][1] = mat.GetPoissonXY();
	D[1][0] = D[0][1];
	D[1][1] = D[0][0];
	Matrix Df(D, 2);
	Matrix ans = Df * mult;
	return ans;
}

//<summary>Loops through the diagonal terms of the stiffness matrix to find the maximum value, and returns 1/100000 of this value</summary>
//<comment>To be used in the Z rotation stiffness term of shell elements</comment>
double ShellElement::GetRotZStiffTerm(Matrix &m) {
	int dim = m.GetDimX();
	double val = 0;

	for (int i = 0; i < dim; i++) {
		if (m.GetMatrixDouble()[i][i] > val) {
			val = m.GetMatrixDouble()[i][i];
		}
	}
	return val / 100000;
}

//<summary>Calculates the jacobian matrix at a given gauss integration point coordinate</summary>
//<e>The e location in natural coordinates</e>
//<n>The n location in natural coordinates</n>
//<c>The c location in natural coordinates</c>
Matrix ShellElement::GetJacobian(double &e, double &n, double &c) const {
	double** J = Matrix::CreateMatrixDouble(3);
	double dNde, dNdn, N;
	Matrix V3 = Get_v3(false, e, n);
	double** v3 = V3.GetMatrixDouble();
	double thick = GetTotalThickness();

	for (int node = 0; node < (_nodeList.size() - 1); node++) { //I guess only the 8 serendipity nodes should be used for this
		dNde = ShellElement::ShapeFunctionNDerivative(node + 1, 'e', e, n);
		dNdn = ShellElement::ShapeFunctionNDerivative(node + 1, 'n', e, n);
		N = ShellElement::ShapeFunctionN(node + 1, e, n);
		double nodeX, nodeY, nodeZ;
		nodeX = (*_nodeList[node]).GetX();
		nodeY = (*_nodeList[node]).GetY();
		nodeZ = (*_nodeList[node]).GetZ();

		J[0][0] += dNde * nodeX + dNde*(thick / 2)*c*v3[0][0];
		J[0][1] += dNde * nodeY + dNde*(thick / 2)*c*v3[1][0];
		J[0][2] += dNde * nodeZ + dNde*(thick / 2)*c*v3[2][0];
		J[1][0] += dNdn * nodeX + dNdn*(thick / 2)*c*v3[0][0];
		J[1][1] += dNdn * nodeY + dNdn*(thick / 2)*c*v3[1][0];
		J[1][2] += dNdn * nodeZ + dNdn*(thick / 2)*c*v3[2][0];
		J[2][0] += N*(thick / 2)*v3[0][0];
		J[2][1] += N*(thick / 2)*v3[1][0];
		J[2][2] += N*(thick / 2)*v3[2][0];

	}
	return Matrix(J, 3, 3);
}

//<summary>Calculates the B matrix at a given gauss integration point coordinate</summary>
//<type>'s' = shear B matrix; 'f' = membrane flexure B matrix</type>
//<e>The e location in natural coordinates</e>
//<n>The n location in natural coordinates</n>
//<c>The c location in natural coordinates</c>
//<Jacobian>The jacobian matrix</Jacobian>
const Matrix ShellElement::GetBMatrix(char type, double &e, double &n, double &c, const Matrix &Jacobian) const {
	
	Matrix V1, V2;
	double** v1;
	double** v2;
	double** BMatrixS = Matrix::CreateMatrixDouble(3, 45);
	double** BMatrixF = Matrix::CreateMatrixDouble(2, 45);
	
	Matrix V3 = Get_v3(false, e, n);
	V1 = Get_v1(V3);
	V2 = Get_v2(V3, V1);
	v1 = V1.GetMatrixDouble();
	v2 = V2.GetMatrixDouble();
	char derivatives[] = { 'e','n','c' };
	double thick = GetTotalThickness();

	for (int node = 0; node < _nodeList.size(); node++) { //for each of the 9 nodes of the element
		double** matrix = Matrix::CreateMatrixDouble(3, 15);
		
		//the next for is to calculate the matrix that mulplies the inverse of the Jacobian
		for (int i = 0; i < 3; i++) { //for each row
			double func = ShellElement::ShapeFunctionN(node + 1, e, n);
			double funcDerivative = ShellElement::ShapeFunctionNDerivative(node + 1, derivatives[i], e, n);
			if (i < 2) { //these rows are derived based on 'e' and 'n'
				matrix[i][0] = funcDerivative;
				matrix[i][3] = funcDerivative * c*(thick / 2)*v1[0][0];
				matrix[i][4] = -funcDerivative * c*(thick / 2)*v2[0][0];
				matrix[i][6] = funcDerivative;
				matrix[i][8] = funcDerivative * c*(thick / 2)*v1[1][0];
				matrix[i][9] = -funcDerivative * c*(thick / 2)*v2[1][0];
				matrix[i][12] = funcDerivative;
				matrix[i][13] = funcDerivative * c*(thick / 2)*v1[2][0];
				matrix[i][14] = -funcDerivative * c*(thick / 2)*v2[2][0];
			}
			else { //these rows are derived based on 'c'
				matrix[i][3] = func * (thick / 2)*v1[0][0];
				matrix[i][4] = -func * (thick / 2)*v2[0][0];
				matrix[i][8] = func * (thick / 2)*v1[1][0];
				matrix[i][9] = -func * (thick / 2)*v2[1][0];
				matrix[i][13] = func *(thick / 2)*v1[2][0];
				matrix[i][14] = -func * (thick / 2)*v2[2][0];
			}
		}

		Matrix full(matrix, 3, 15);

		if (type == 's') { // Calculation of the Shear Matrix Terms
			for (int i = 0; i < 5; i++) { //for each term of the B Matrix
				Matrix m = MatrixOperation::ExtractSetOfColumns(full, 3, i, 5);
				Matrix invJacob = MatrixOperation::GetInverseWithAdjoint(Jacobian);
				Matrix transf = GetTransformationMatrix3x3Hyrnyk(Jacobian);
				Matrix transpTransf = MatrixOperation::Transpose(transf);
				Matrix m2 = transpTransf * invJacob * m * transf;
				BMatrixS[0][i + node * 5] = m2.GetMatrixDouble()[1][0] + m2.GetMatrixDouble()[0][1];
				BMatrixS[1][i + node * 5] = m2.GetMatrixDouble()[2][0] + m2.GetMatrixDouble()[0][2];
				BMatrixS[2][i + node * 5] = m2.GetMatrixDouble()[2][1] + m2.GetMatrixDouble()[1][2];
			}
		}
		else { // Calculation of the Flexural Matrix Terms
			for (int i = 0; i < 5; i++) { //for each term of the B Matrix
				Matrix m = MatrixOperation::ExtractSetOfColumns(full, 3, i, 5);
				Matrix invJacob = MatrixOperation::GetInverseWithAdjoint(Jacobian);
				Matrix transf = GetTransformationMatrix3x3Hyrnyk(Jacobian);
				Matrix transpTransf = MatrixOperation::Transpose(transf);
				Matrix m2 = transpTransf * invJacob * m * transf;
				BMatrixF[0][i + node * 5] = m2.GetMatrixDouble()[0][0];
				BMatrixF[1][i + node * 5] = m2.GetMatrixDouble()[1][1];
			}
		}
	}
	if (type == 's') {
		Matrix::DestroyMatrixDouble(BMatrixF, 2);
		return Matrix(BMatrixS, 3, 45);
	}
	else {
		Matrix::DestroyMatrixDouble(BMatrixS, 3);
		return Matrix(BMatrixF, 2, 45);
	}
}

//<summary>Calculates the B matrix at a given gauss integration point coordinate</summary>
//<e>The e location in natural coordinates</e>
//<n>The n location in natural coordinates</n>
//<c>The c location in natural coordinates</c>
//<Jacobian>The jacobian matrix</Jacobian>
//<comment>To be used when shear and membrane-flexure matrices are calculated together</comment>
const Matrix ShellElement::GetTotalBMatrix(double &e, double &n, double &c, const Matrix &Jacobian) const {

	Matrix V1, V2;
	double** v1;
	double** v2;
	double** BMatrix = Matrix::CreateMatrixDouble(5, 45);

	Matrix V3 = Get_v3(false, e, n);
	V1 = Get_v1(V3);
	V2 = Get_v2(V3, V1);
	v1 = V1.GetMatrixDouble();
	v2 = V2.GetMatrixDouble();
	char derivatives[] = { 'e','n','c' };
	double thick = GetTotalThickness();

	for (int node = 0; node < _nodeList.size(); node++) { //for each of the 9 nodes of the element
		double** matrix = Matrix::CreateMatrixDouble(3, 15);

		//the next for is to calculate the matrix that mulplies the inverse of the Jacobian
		for (int i = 0; i < 3; i++) { //for each row
			double func = ShellElement::ShapeFunctionN(node + 1, e, n);
			double funcDerivative = ShellElement::ShapeFunctionNDerivative(node + 1, derivatives[i], e, n);
			if (i < 2) { //these rows are derived based on 'e' and 'n'
				matrix[i][0] = funcDerivative;
				matrix[i][3] = funcDerivative * c*(thick / 2)*v1[0][0];
				matrix[i][4] = -funcDerivative * c*(thick / 2)*v2[0][0];
				matrix[i][6] = funcDerivative;
				matrix[i][8] = funcDerivative * c*(thick / 2)*v1[1][0];
				matrix[i][9] = -funcDerivative * c*(thick / 2)*v2[1][0];
				matrix[i][12] = funcDerivative;
				matrix[i][13] = funcDerivative * c*(thick / 2)*v1[2][0];
				matrix[i][14] = -funcDerivative * c*(thick / 2)*v2[2][0];
			}
			else { //these rows are derived based on 'c'
				matrix[i][3] = func * (thick / 2)*v1[0][0];
				matrix[i][4] = -func * (thick / 2)*v2[0][0];
				matrix[i][8] = func * (thick / 2)*v1[1][0];
				matrix[i][9] = -func * (thick / 2)*v2[1][0];
				matrix[i][13] = func * (thick / 2)*v1[2][0];
				matrix[i][14] = -func * (thick / 2)*v2[2][0];
			}
		}

		Matrix full(matrix, 3, 15);


		for (int i = 0; i < 5; i++) { //for each term of the B Matrix
			Matrix m = MatrixOperation::ExtractSetOfColumns(full, 3, i, 5);
			Matrix invJacob = MatrixOperation::GetInverseWithAdjoint(Jacobian);
			//is this where the transformation is performed in the B matrix and, consequently, on the stiffness matrix?
			Matrix transf = GetTransformationMatrix3x3Hyrnyk(Jacobian);
			Matrix transpTransf = MatrixOperation::Transpose(transf);
			Matrix m2 = transpTransf * invJacob * m * transf;

			BMatrix[0][i + node * 5] = m2.GetMatrixDouble()[0][0];
			BMatrix[1][i + node * 5] = m2.GetMatrixDouble()[1][1];
			BMatrix[2][i + node * 5] = m2.GetMatrixDouble()[1][0] + m2.GetMatrixDouble()[0][1];
			BMatrix[3][i + node * 5] = m2.GetMatrixDouble()[2][0] + m2.GetMatrixDouble()[0][2];
			BMatrix[4][i + node * 5] = m2.GetMatrixDouble()[2][1] + m2.GetMatrixDouble()[1][2];
		}
	}

	return Matrix(BMatrix, 5, 45);
}

Matrix ShellElement::GetMassMatrixTheory() {
	double DOF = 6;

	double gaussPoints[3][3] = { {-sqrt(0.6), 0, sqrt(0.6)},{-sqrt(0.6), 0, sqrt(0.6)}, {-sqrt(0.6), 0, sqrt(0.6)} };
	double w[3] = { 5.0 / 9.0,8.0 / 9.0,5.0 / 9.0 };

	Matrix mass(_nodeList.size()*DOF, _nodeList.size()*DOF);
	
	if (_massRho != 0) {
		
		int count = 0;
		for (int node = 0; node < _nodeList.size() - 1; node++) { //node 9 doesn't have because all DOFs are rotational
			double mii = 0;
			double firstIntegral = 0;
			double secondIntegral = 0;
			double thirdIntegral = 0;

			for (int i = 0; i < 3; i++) { // for the first gauss integration
				for (int j = 0; j < 3; j++) { // for the second gauss integration
					double mult = 0;
					for (int node2 = 0; node2 < _nodeList.size() - 1; node2++) { //for each of the 9 nodes of the element. Not depending on 'k' yet.
						double Nk = ShapeFunctionN(node2 + 1, gaussPoints[0][i], gaussPoints[1][j]);
						mult += Nk * _massRho * Nk;
					}
					for (int k = 0; k < 3; k++) { // for the third gauss integration this should be relative to the number of layers the element has
						Matrix jacob3x3 = GetJacobian(gaussPoints[0][i], gaussPoints[1][j], gaussPoints[2][k]);
						double detJacob = MatrixOperation::DeterminantWithCofactor(jacob3x3, 3);

						double Ni = ShapeFunctionN(node + 1, gaussPoints[0][i], gaussPoints[1][j]);

						firstIntegral += Ni * _massRho * Ni * detJacob * w[i] * w[j] * w[k];
						secondIntegral += _massRho * detJacob * w[i] * w[j] * w[k];
						thirdIntegral += mult * detJacob * w[i] * w[j] * w[k];
					}
				}
			}
			mii += firstIntegral * secondIntegral / thirdIntegral;

			int index = node * DOF;
			mass.GetMatrixDouble()[index][index] = mii;
			mass.GetMatrixDouble()[index + 1][index + 1] = mii;
			mass.GetMatrixDouble()[index + 2][index + 2] = mii;
		}
	}

	return mass;
}

//this can be used if the elements are quadratic (rectangles or squares). If the element is of 'general shape' then the theoretical formulation should be used and integrated
Matrix ShellElement::GetMassMatrix() { 
	double totalMass = GetWidth() * GetHeight() * _thick * _massRho;
	double extremeNodes = 0.0395;
	double midNodes = 0.215;
	double DOF = 6;

	Matrix mass(_nodeList.size()*DOF, _nodeList.size()*DOF);
	if (_massRho != 0) {
		for (int node = 0; node < _nodeList.size() - 1; node++) { //node 9 doesn't have because all DOFs are rotational

			int index = node * DOF;
			if (node < 4) {
				mass.GetMatrixDouble()[index][index] = totalMass * extremeNodes;
				mass.GetMatrixDouble()[index + 1][index + 1] = totalMass * extremeNodes;
				mass.GetMatrixDouble()[index + 2][index + 2] = totalMass * extremeNodes;
			}
			else {
				mass.GetMatrixDouble()[index][index] = totalMass * midNodes;
				mass.GetMatrixDouble()[index + 1][index + 1] = totalMass * midNodes;
				mass.GetMatrixDouble()[index + 2][index + 2] = totalMass * midNodes;
			}
		}
	}
	return mass;
}

//<summary>Calculates the global stiffness matrix for the element following the procedure described in the Hinton and Owen, 1984 book</summary>
//<comment>Not used because. It considers different gauss integration points for membrane and shear matrices.</comment>
Matrix ShellElement::GetGlobalStiffMatrixFollowingBook() { 
	Matrix shearK(54, 54);
	Matrix transfK(54, 54);
	Matrix transfFlex(54, 54);
	Matrix shearK3x3(54, 54);
	Matrix membFlexK(54, 54);

	//Calculate the gauss points associate with the thickness of each layer
	/*
	std::vector<double> gaussPointsThick;
	gaussPointsThick.reserve(_layers);
	double step = (layerThick / (eleThick / 2));
	double initPos = -1 - step / 2;
	for (int i = 0; i < _layers; i++) {
		gaussPointsThick.emplace_back(initPos + (i + 1) * step);
	}
	*/

	double gaussPoints[3][3] = { {-sqrt(0.6), 0, sqrt(0.6)},{-sqrt(0.6), 0, sqrt(0.6)}, {-sqrt(0.6), 0, sqrt(0.6)} };
	//double gaussPoints[3][2] = { {-1 / sqrt(3), 1 / sqrt(3)},{-1 / sqrt(3), 1 / sqrt(3)}, {-1 / sqrt(3), 1 / sqrt(3)} };
	//int w = 1;
	double w[3] = { 5.0 / 9.0,8.0 / 9.0,5.0 / 9.0 };
	//2x2 gauss integration
	for (int i = 0; i < 3; i++) { // for the first gauss integration
		for (int j = 0; j < 3; j++) { // for the second gauss integration
			for (int k = 0; k < 3; k++) { // for the third gauss integration this should be relative to the number of layers the element has
				//For each point in the Gauss Integration
				Matrix shearDMatrix = ShellElement::GetShearDMatrix(_matList[0]);
				Matrix jacob2x2 = ShellElement::GetJacobian(gaussPoints[0][i], gaussPoints[1][j], gaussPoints[2][k]);
				Matrix BMatrix = GetBMatrix('s', gaussPoints[0][i], gaussPoints[1][j], gaussPoints[2][k], jacob2x2);
				Matrix transpB = MatrixOperation::Transpose(BMatrix);
				double detJacob = MatrixOperation::DeterminantWithCofactor(jacob2x2, 3);
				double constmult = w[i] * w[j] * w[k];
				//double constmult = w * w;
				Matrix localShearB = transpB * shearDMatrix * BMatrix * detJacob * constmult;

				double zTerm = ShellElement::GetRotZStiffTerm(localShearB);
				for (int i = 0; i < 9; i++) {
					localShearB = MatrixOperation::AddLineAndColumnWithTermAtPosition(localShearB, (i + 1) * 6, (i + 1) * 6, zTerm);
				}

				Matrix transf = ShellElement::GetTransformationMatrix(jacob2x2);
				Matrix transpTransf = MatrixOperation::Transpose(transf);
				Matrix final = transpTransf * localShearB * transf;

				shearK += final; //already at 54x54
			}
		}
	}

	double w2[3] = { 5.0 / 9.0,8.0 / 9.0,5.0 / 9.0 };
	//double gaussPoints3[2][3] = { {-sqrt(0.6), 0, sqrt(0.6)},{-sqrt(0.6), 0, sqrt(0.6)}};
	double gaussPoints3[3][3] = { {-sqrt(0.6), 0, sqrt(0.6)},{-sqrt(0.6), 0, sqrt(0.6)}, {-sqrt(0.6), 0, sqrt(0.6)} };

	//3x3 gauss integration
	for (int i = 0; i < 3; i++) { // for the first gauss integration
		for (int j = 0; j < 3; j++) { // for the second gauss integration
			//Tranforming the shear B matrix to 3x3
			/*
			for (int l = 0; l < 4; l++) {
				Matrix mult = shearK * ShellElement::ShapeFunctionNr(l, gaussPoints3[0][i], gaussPoints3[1][j]);
				shearK3x3 += mult;
			}
			*/
			for (int k = 0; k < 3; k++) {  // for the third gauss integration
				//For each point in the Gauss Integration
				Matrix flexureDMatrix = ShellElement::GetMembFlexDMatrix(_matList[0]);
				Matrix jacob3x3 = ShellElement::GetJacobian(gaussPoints3[0][i], gaussPoints3[1][j], gaussPoints3[2][k]);
				Matrix flexureB = GetBMatrix('f', gaussPoints3[0][i], gaussPoints3[1][j], gaussPoints3[2][k], jacob3x3);
				Matrix transpB = MatrixOperation::Transpose(flexureB);
				double detJacob = MatrixOperation::DeterminantWithCofactor(jacob3x3, 3);
				double constmult = w2[i] * w2[j] * w2[k];
				Matrix localFlexureB = transpB * flexureDMatrix * flexureB * detJacob * constmult;

				double zTerm = ShellElement::GetRotZStiffTerm(localFlexureB);
				for (int i = 0; i < 9; i++) {
					localFlexureB = MatrixOperation::AddLineAndColumnWithTermAtPosition(localFlexureB, (i + 1) * 6, (i + 1) * 6, zTerm);
				}

				Matrix transf = ShellElement::GetTransformationMatrix(jacob3x3);
				Matrix transpTransf = MatrixOperation::Transpose(transf);
				Matrix final = transpTransf * localFlexureB * transf;
				membFlexK += final;
			}
		}
	}

	Matrix stiffMatrix = shearK + membFlexK;
	//Matrix stiffMatrix = shearK3x3 + membFlexK;

	return stiffMatrix;
}

//<summary>Calculates the global stiffness matrix for the element</summary>
//<comment>Not used because it considers shear and membrane matrices separately. It considers the same number of gauss integration points for membrane and shear matrices.</comment>
Matrix ShellElement::GetGlobalStiffMatrixSeparate() {
	Matrix stiffMatrix(54, 54);


	//Calculate the gauss points associate with the thickness of each layer
	/*
	std::vector<double> gaussPointsThick;
	double step = (layerThick / (eleThick / 2));
	double initPos = -1 - step / 2;
	for (int i = 0; i < _layers; i++) {
		gaussPointsThick.push_back(initPos + (i + 1) * step);
	}
	*/

	double gaussPoints[3][3] = { {-sqrt(0.6), 0, sqrt(0.6)},{-sqrt(0.6), 0, sqrt(0.6)}, {-sqrt(0.6), 0, sqrt(0.6)} };
	//double gaussPoints[3][2] = { {-1 / sqrt(3), 1 / sqrt(3)},{-1 / sqrt(3), 1 / sqrt(3)}, {-1 / sqrt(3), 1 / sqrt(3)} };
	//int w = 1;
	double w[3] = { 5.0 / 9.0,8.0 / 9.0,5.0 / 9.0 };
	//2x2 gauss integration
	for (int i = 0; i < 3; i++) { // for the first gauss integration
		for (int j = 0; j < 3; j++) { // for the second gauss integration
			for (int k = 0; k < 3; k++) { // for the third gauss integration this should be relative to the number of layers the element has
				//For each point in the Gauss Integration
				Matrix shearDMatrix = ShellElement::GetShearDMatrix(_matList[0]);
				Matrix jacob3x3 = ShellElement::GetJacobian(gaussPoints[0][i], gaussPoints[1][j], gaussPoints[2][k]);
				Matrix BMatrix = GetBMatrix('s', gaussPoints[0][i], gaussPoints[1][j], gaussPoints[2][k], jacob3x3);
				Matrix transpB = MatrixOperation::Transpose(BMatrix);
				double detJacob = MatrixOperation::DeterminantWithCofactor(jacob3x3, 3);
				double constmult = w[i] * w[j] * w[k];
				//double constmult = w * w;
				Matrix localShearB = transpB * shearDMatrix * BMatrix * detJacob * constmult;

				Matrix flexureDMatrix = ShellElement::GetMembFlexDMatrix(_matList[0]);
				Matrix flexureB = GetBMatrix('f', gaussPoints[0][i], gaussPoints[1][j], gaussPoints[2][k], jacob3x3);
				Matrix transpBFlex = MatrixOperation::Transpose(flexureB);
				Matrix localFlexureB = transpBFlex * flexureDMatrix * flexureB * detJacob * constmult;

				double zTermShear = ShellElement::GetRotZStiffTerm(localShearB);
				double zTermFlexure = ShellElement::GetRotZStiffTerm(localFlexureB);
				for (int i = 0; i < 9; i++) {
					localShearB = MatrixOperation::AddLineAndColumnWithTermAtPosition(localShearB, (i + 1) * 6, (i + 1) * 6, zTermShear);
					localFlexureB = MatrixOperation::AddLineAndColumnWithTermAtPosition(localFlexureB, (i + 1) * 6, (i + 1) * 6, zTermFlexure);
				}

				Matrix add = localShearB + localFlexureB;

				Matrix transf = ShellElement::GetTransformationMatrix(jacob3x3);
				Matrix transpTransf = MatrixOperation::Transpose(transf);

				Matrix final = transpTransf * add * transf;

				stiffMatrix += final;
			}
		}
	}

	return stiffMatrix;
}

//<summary>Calculates the global stiffness matrix for the element</summary>
//<comment>It considers the same number of gauss integration points for membrane and shear matrices.</comment>
const Matrix ShellElement::GetGlobalStiffMatrix() const {
	Matrix stiffMatrix(54, 54);

	int pointLayer = 3; //number of integration points per layer

	//Calculate the gauss points associated with the thickness of each layer
	/*
	std::vector<std::vector<double>> gaussThick;
	for (int i = 0; i < _layers; i++) { //starting from the bottom layer
		double bot, top;
		bot = -1 + i * _thick / ((_thick*_layers) / 2);
		top = bot + _thick / ((_thick*_layers) / 2);
		std::vector<double> points;
		double size = (top - bot);
		points.push_back(bot + size/2 - (sqrt(0.6)*size/2));
		points.push_back(bot + size / 2);
		points.push_back(top - size/2 + sqrt(0.6)*size/2);
		gaussThick.push_back(points);
	}
	*/
	//For some reason, the layer approach doesn't work. It always gives 1/3 of the answer it should give. I can't figure out why.

	double gaussPoints[3][3] = { {-sqrt(0.6), 0, sqrt(0.6)},{-sqrt(0.6), 0, sqrt(0.6)}, {-sqrt(0.6), 0, sqrt(0.6)} };
	//int w = 1;
	double w[3] = { 5.0 / 9.0,8.0 / 9.0,5.0 / 9.0 };
	//3x3 gauss integration
	int cCounter = 0;

	for (int i = 0; i < 3; i++) { // for the first gauss integration
		for (int j = 0; j < 3; j++) { // for the second gauss integration
			//for (int layer = 0; layer < _layers; layer++) {
				for (int k = 0; k < 3; k++) { // for the third gauss integration this should be relative to the number of layers the element has
					//For each point in the Gauss Integration
					//Matrix* DMatrix = &ShellElement::GetDMatrix(_matList[floor(k / 3)]); //to be used with layers
					const Matrix DMatrix = ShellElement::GetDMatrix(_matList[0]);
					const Matrix jacob3x3 = GetJacobian(gaussPoints[0][i], gaussPoints[1][j], gaussPoints[2][k]);
					const Matrix BMatrix = GetTotalBMatrix(gaussPoints[0][i], gaussPoints[1][j], gaussPoints[2][k], jacob3x3);
					Matrix transpB = MatrixOperation::Transpose(BMatrix);
					double detJacob = MatrixOperation::DeterminantWithCofactor(jacob3x3, 3);
					double constmult = w[i] * w[j] * w[k];
					Matrix localB = (transpB) * (DMatrix) * (BMatrix) * (detJacob) * (constmult);

					double zTerm = ShellElement::GetRotZStiffTerm(localB);
					for (int i = 0; i < 9; i++) {
						localB = MatrixOperation::AddLineAndColumnWithTermAtPosition(localB, (i + 1) * 6, (i + 1) * 6, zTerm);
					}

					//If this is used to transform the matrix, the results are wrong. Somewhere else in the procedure the transformation is performed. I believe it is done in the B matrix.
					//Matrix* transf = &ShellElement::GetTransformationMatrix(*jacob3x3);
					//Matrix* transpTransf = &Matrix::Transpose(*transf);

					//Matrix firstMultFinal = (*transpTransf) * localB;
					//Matrix final = firstMultFinal * (*transf);

					stiffMatrix += localB;
				}
			//}
		}
	}

	return stiffMatrix;
}

//<summary>Assembles the complete global matrix of a list of elements using threads</summary>
//<vecEle>The list of elements used</vecEle>
//<complete>The complete matrix, defined somewhere else, used to store the values</complete>
//<sup>The list of supports in the structure</sup>
//<mu>A way of protecting the different threads of acessing the same information at the same time</mu>
void ShellElement::AssembleCompleteGlobalMatrixThreads(const std::vector<ShellElement*>* vecEle, Matrix& complete, std::mutex& mu) {
	for (int k = 0; k < vecEle->size(); k++) { // for each element
		const ShellElement* ele = (*vecEle)[k];
		Matrix global = ele->GetGlobalStiffMatrix(); //45x45 matrix returns
		std::vector<std::vector<int>> vec = ele->GetGlobalDOFVector();
		int runs = vec.size(); //number of DOFs not restrained
		for (int i = 0; i < runs; i++) { //for all the lines in the local stiffness matrix
			int index1 = vec[i][0];
			for (int j = 0; j < runs; j++) { //for all the columns
				int index2 = vec[j][0];
				std::lock_guard<std::mutex> lock(mu);
				complete.GetMatrixDouble()[index1][index2] += global.GetMatrixDouble()[vec[i][1]][vec[j][1]];
			}
		}
	}
}

void ShellElement::AssembleCompleteRestrictedGlobalMatrixThreads(const std::vector<ShellElement*>* vecEle, Matrix& complete, const int* unrestrictDOFs, std::mutex& mu) {
	for (int k = 0; k < vecEle->size(); k++) { // for each element
		const ShellElement* ele = (*vecEle)[k];
		std::vector<std::vector<int>> vecRest = ele->GetGlobalRestDOFVector();
		if (vecRest.size() != 0) {
			std::vector<std::vector<int>> vecGlobal = ele->GetGlobalDOFVector();
			Matrix global = ele->GetGlobalStiffMatrix(); //45x45 matrix returns
			int runs = vecRest.size(); //number of DOFs not restrained
			for (int i = 0; i < runs; i++) { //for all the lines in the local stiffness matrix
				int index1 = vecRest[i][0];
				for (int j = 0; j < vecGlobal.size(); j++) {
					int index2 = vecGlobal[j][0];
					std::lock_guard<std::mutex> lock(mu);
					complete.GetMatrixDouble()[index1][index2] += global.GetMatrixDouble()[vecRest[i][1]][vecGlobal[j][1]];
				}
				for (int j = 0; j < runs; j++) { //for all the columns
					int index2 = vecRest[j][0] + *unrestrictDOFs;
					std::lock_guard<std::mutex> lock(mu);
					complete.GetMatrixDouble()[index1][index2] += global.GetMatrixDouble()[vecRest[i][1]][vecRest[j][1]];
				}
			}
		}
	}
}

//<summary>Calculates the displacement of the ninth node of the element</summary>
void ShellElement::GetNinthNodeDisplacement(Matrix &totalDisplacementMatrix, std::vector<ShellElement> &listElements) {
	int DOF = 6;
	double** matrix = totalDisplacementMatrix.GetMatrixDouble();
	for (int i = 0; i < listElements.size(); i++) { //for each element
		double x = 0, y = 0, z = 0, rx = 0, ry = 0;
		double e = 0, n = 0;
		std::vector<Node *> nodeList = listElements[i].GetNodeList();
		int ninthNodeID = (*nodeList[8]).GetID();

		for (int j = 0; j < 8; j++) { //use the other 8 nodes to calculate the displacement of the ninth node
			double shape = ShellElement::ShapeFunctionN(j + 1, e, n);
			int nodeID = (*nodeList[j]).GetID();
			x += shape * matrix[(nodeID - 1)*DOF][0];
			y += shape * matrix[(nodeID - 1)*DOF + 1][0];
			z += shape * matrix[(nodeID - 1)*DOF + 2][0];
			rx += shape * matrix[(nodeID - 1)*DOF + 3][0];
			ry += shape * matrix[(nodeID - 1)*DOF + 4][0];
		}

		totalDisplacementMatrix.GetMatrixDouble()[(ninthNodeID - 1)*DOF][0] = x;
		totalDisplacementMatrix.GetMatrixDouble()[(ninthNodeID - 1)*DOF + 1][0] = y;
		totalDisplacementMatrix.GetMatrixDouble()[(ninthNodeID - 1)*DOF + 2][0] = z;
		totalDisplacementMatrix.GetMatrixDouble()[(ninthNodeID - 1)*DOF + 3][0] += rx;
		totalDisplacementMatrix.GetMatrixDouble()[(ninthNodeID - 1)*DOF + 4][0] += ry;
	}
}

// see if the local DOF is a rotational DOF in the local coordinate system
bool ShellElement::IsDOFRotational(int &DOF) {
	for (int i = 0; i < _nodeList.size(); i++) {
		double val1, val2, val3;
		val1 = 3 + 6 * i;
		val2 = 4 + 6 * i;
		val3 = 5 + 6 * i;
		if (DOF == val1 || DOF == val2 || DOF == val3) {
			return true;
		}
	}
	return false;
}

std::vector<std::vector<int>> ShellElement::CrescentOrderDOFVector(std::vector<std::vector<int>>& oriVec) {
	std::vector<std::vector<int>> newVec;
	newVec.reserve(oriVec.size());
	newVec = oriVec;

	int i, j;
	for (i = 0; i < newVec.size() - 1; i++) {
		// Last i elements are already in place    
		for (j = 0; j < newVec.size() - i - 1; j++) {
			if (newVec[j][0] > newVec[j + 1][0]) {
				int temp1 = newVec[j][0];
				int temp2 = newVec[j][1];

				newVec[j][0] = newVec[j + 1][0];
				newVec[j][1] = newVec[j + 1][1];
				newVec[j + 1][0] = temp1;
				newVec[j + 1][1] = temp2;
			}
		}
	}
	return newVec;
}


Matrix ShellElement::ConvertAccFromReducedToTotal(Matrix &acc, std::vector<std::vector<int>> &totalMassDOFVec, int size) {
	int DOF = 6;
	int count = 0;

	Matrix ans(size, 1);

	bool hasMass = true;
	for (int k = 0; k < size; k++) {
		/*
		for (int i = 0; i < vecEle.size(); i++) { //for each shell element
			ShellElement* ele = &vecEle[i];
			std::vector<std::vector<int>> vec = ele->GetGlobalMassDOFVector();
			std::vector<std::vector<int>> vec1 = ShellElement::CrescentOrderDOFVector(vec); //DOFs without mass
			if (vec1.size() != 0) { //if there are any DOF not restricted in this element (possible)
				for (int j = 0; j < vec1.size(); j++) { //for each not-restricted DOF of the element
					if (k == vec1[j][0]) {
						hasMass = false;
					}
				}
			}
		}
		*/
		for (int i = 0; i < totalMassDOFVec.size(); i++) { //for each shell element
			if (k == totalMassDOFVec[i][0]) {
				hasMass = false;
			}
		}

		if (hasMass) {
			ans.GetMatrixDouble()[k][0] = acc.GetMatrixDouble()[count][0];
			count++;
		}
		else {
			hasMass = true;
		}
	}

	return ans;
}

std::vector<std::vector<int>> ShellElement::GetTotalGlobalMassDOFVector(std::vector<ShellElement> &vecEle) {
	std::vector<std::vector<int>> ans;
	for (int i = 0; i < vecEle.size(); i++) { //for each shell element
		ShellElement* ele = &vecEle[i];
		std::vector<std::vector<int>> vec = ele->GetGlobalMassDOFVector();
		ans.insert(ans.end(), vec.begin(), vec.end());
	}

	std::vector<std::vector<int>> organized = ShellElement::CrescentOrderDOFVector(ans);
	std::vector<std::vector<int>> answer;
	int count = 0;
	for (int i = 0; i < organized.size(); i++) {
		bool test;
		if (i == organized.size() - 1) {//if last term
			test = false;
		}
		else {
			test = (organized[i][0] == organized[i + 1][0]);
		}
		if (!test) {
			answer.push_back(organized[i]);
		}
	}
	return answer;
}

Matrix ShellElement::CondensedReducedStiffMatrixForModal(Matrix &m, std::vector<std::vector<int>> &totalMassDOFVec) {

	int DOF = 6;
	int count = 0;

	Matrix matrix(MatrixOperation::CopyMatrixDouble(m), m.GetDimX(), m.GetDimY());
	/*
	for (int i = 0; i < vecEle.size(); i++) { //for each shell element
		ShellElement* ele = &vecEle[i];
		std::vector<std::vector<int>> vec = ele->GetGlobalMassDOFVector();
		std::vector<std::vector<int>> vec1 = ShellElement::CrescentOrderDOFVector(vec);
		if (vec1.size() != 0) { //if there are any DOF not restricted in this element (possible)
			for (int j = 0; j < vec1.size(); j++) { //for each not-restricted DOF of the element
				int DOF = vec1[j][0]; //index on the global reduced matrix
				Matrix::MoveLineAndColumnToEnd(matrix, DOF - count); //move the line and column from the original DOF in the global reduced matrix to the last line and column
				count++;
			}
		}
	}
	*/
	for (int i = 0; i < totalMassDOFVec.size(); i++) { //for each shell element
		int DOF = totalMassDOFVec[i][0]; //index on the global reduced matrix
		MatrixOperation::MoveLineAndColumnToEnd(matrix, DOF - count); //move the line and column from the original DOF in the global reduced matrix to the last line and column
		count++;
	}
	Matrix k11 = MatrixOperation::ExtractMatrixBasedOnLineAndColumns(matrix, 0, matrix.GetDimX() - count, 0, matrix.GetDimY() - count);
	Matrix k12 = MatrixOperation::ExtractMatrixBasedOnLineAndColumns(matrix, 0, matrix.GetDimX() - count, matrix.GetDimY() - count, count);
	Matrix k21 = MatrixOperation::Transpose(k12);
	Matrix k22 = MatrixOperation::ExtractMatrixBasedOnLineAndColumns(matrix, matrix.GetDimX() - count, count, matrix.GetDimY() - count, count);
	Matrix k22Inv = MatrixOperation::GetInverse(k22);

	Matrix ans = k11 - k12 * k22Inv * k21;

	return ans;
}

Matrix ShellElement::ReducedAccelerationForceMatrix(Matrix &m, std::vector<std::vector<int>>& totalMassDOFVec) {
	int DOF = 6;
	int count = 0;

	Matrix matrix(MatrixOperation::CopyMatrixDouble(m), m.GetDimX(), m.GetDimY());

	/*
	for (int i = 0; i < vecEle.size(); i++) { //for each shell element
		ShellElement* ele = &vecEle[i];
		std::vector<std::vector<int>> vec = ele->GetGlobalMassDOFVector();
		std::vector<std::vector<int>> vec1 = ShellElement::CrescentOrderDOFVector(vec);
		if (vec1.size() != 0) { //if there are any DOF not restricted in this element (possible)
			for (int j = 0; j < vec1.size(); j++) { //for each not-restricted DOF of the element
				int DOF = vec1[j][0]; //index on the global reduced matrix
				Matrix::MoveLineToEnd(matrix, DOF - count); //move the line and column from the original DOF in the global reduced matrix to the last line and column
				count++;
			}
		}
	}
	*/

	for (int i = 0; i < totalMassDOFVec.size(); i++) { //for each shell element
		int DOF = totalMassDOFVec[i][0]; //index on the global reduced matrix
		MatrixOperation::MoveLineToEnd(matrix, DOF - count); //move the line and column from the original DOF in the global reduced matrix to the last line and column
		count++;
	}

	Matrix k11 = MatrixOperation::ExtractMatrixBasedOnLineAndColumns(matrix, 0, matrix.GetDimX() - count, 0, 1);

	return k11;
}

void ShellElement::AssembleCompleteGlobalMassMatrixThreads(std::vector<ShellElement> &vecEle, Matrix& complete, std::vector<Support> &sup, std::mutex& mu) {
	int DOF = 6; //# of degrees of freedom in each shell element

	for (int k = 0; k < vecEle.size(); k++) { // for each element
		ShellElement* ele = &vecEle[k];
		Matrix global = (ele->GetMassMatrixTheory()); //45x45 matrix returns
		std::vector<std::vector<int>> vec = ele->GetGlobalDOFVector();
		int runs = vec.size(); //number of DOFs not restrained
		for (int i = 0; i < runs; i++) { //for all the lines in the local stiffness matrix
			int index1 = vec[i][0];
			for (int j = 0; j < runs; j++) { //for all the columns
				int index2 = vec[j][0];
				std::lock_guard<std::mutex> lock(mu);
				complete.GetMatrixDouble()[index1][index2] += global.GetMatrixDouble()[vec[i][1]][vec[j][1]];
			}
		}
	}
}

Matrix ShellElement::GetMassMatrixNonZeroMassOnly(Matrix &m, std::vector<std::vector<int>> &totalMassDOFVec) {

	int DOF = 6;
	int count = 0;
	int massCount = 0;
	Matrix k11(m.GetDimX() - totalMassDOFVec.size());
	Matrix matrix(MatrixOperation::CopyMatrixDouble(m), m.GetDimX(), m.GetDimY());

	/*
	for (int i = 0; i < vecEle.size(); i++) { //for each shell element
		ShellElement* ele = &vecEle[i];
		std::vector<std::vector<int>> vec = ele->GetGlobalMassDOFVector();
		std::vector<std::vector<int>> vec1 = ShellElement::CrescentOrderDOFVector(vec);
		if (vec1.size() != 0) { //if there are any DOF not restricted in this element (possible)
			for (int j = 0; j < vec1.size(); j++) { //for each not-restricted DOF of the element
				int DOF = vec1[j][0]; //index on the global reduced matrix
				Matrix::MoveLineAndColumnToEnd(matrix, DOF - count); //move the line and column from the original DOF in the global reduced matrix to the last line and column
				count++;
			}
		}
	}
	*/


	for (int i = 0; i < m.GetDimX(); i++) { //for each shell element
		if (m.GetMatrixDouble()[i][i] != 0) {
			k11.GetMatrixDouble()[massCount][massCount] = m.GetMatrixDouble()[i][i];
			massCount++;
		}
	}

		/*
	for (int i = 0; i < totalMassDOFVec.size(); i++) { //for each shell element
		int DOF = totalMassDOFVec[i][0]; //index on the global reduced matrix
		MatrixOperation::MoveLineAndColumnToEnd(matrix, DOF - count); //move the line and column from the original DOF in the global reduced matrix to the last line and column
		count++;
	}
		*/

	//Matrix k11 = MatrixOperation::ExtractMatrixBasedOnLineAndColumns(matrix, 0, matrix.GetDimX() - count, 0, matrix.GetDimY() - count);
	return k11;
}

bool ShellElement::IsDOFInTheNinthNode(int DOF, std::vector<ShellElement> &vecEle) {
	//TODO:GlobalDOFVector is deprecated. Need to investigate why is this used here.
	/*
	for (int i = 0; i < vecEle.size(); i++) {
		std::vector<int> vec = vecEle[i].GlobalDOFVector();
		for (int j = 0; j < 3; j++) {
			if (vec[vec.size() - 1 - 3 - j] == DOF) {
				return true;
			}
		}
	}
	*/
	return false;
}

ShellElement::~ShellElement()
{
}