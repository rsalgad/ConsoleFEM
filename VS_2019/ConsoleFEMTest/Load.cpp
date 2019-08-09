#include "pch.h"
#include "Load.h"
#include "Support.h"
#include "MatrixOperation.h"
#include "StructureManager.h"
#include "PreAnalysisSetUp.h"
#include <iostream>


Load::Load()
{
}

Load::Load(int ID, int nodeID)
{
	_ID = ID;
	_nodeID = nodeID;
}

Load::Load(int ID, int nodeID, std::string status)
{
	_ID = ID;
	_nodeID = nodeID;
	_status = status;
}

std::string Load::ToString()
{
	std::string load = "";
	load += "(";
	load += std::to_string(_ID);
	load += ")";
	load += "(";
	load += "Node: ";
	load += std::to_string(_nodeID);
	load += ", ";
	load += "Load Type: ";
	load += _status;
	load += ", ";
	for (int i = 0; i < _load[0].size(); i++) {
		
		if (i != 0)
		{
			load += ", ";
		}
		
		load += "Component ";
		load += (i + 1);
		load += ": ";
		int dir = _load[i][0];
		
		switch (dir)
		{
			case 1:
				load += "Fx = ";
				break;
			case 2: 
				load += "Fy = ";
				break;
			case 3:
				load += "Fz = ";
				break;
			case 4:
				load += "Mx = ";
				break;
			case 5:
				load += "My = ";
				break;
			case 6:
				load += "Mz = ";
				break;
			default:
				break;
		}

		load += std::to_string(_load[i][1]);

	}
	
	load += ")";

	return load;
}

int Load::GetID()
{
	return _ID;
}

int Load::GetNode()
{
	return _nodeID;
}

std::string Load::GetStatus() {
	return _status;
}

void Load::SetFx(double val) {
	std::vector<double> vec = { 1, val };
	_load.push_back(vec);
}

void Load::SetFy(double val) {
	std::vector<double> vec = { 2, val };
	_load.push_back(vec);
}

void Load::SetFz(double val) {
	std::vector<double> vec = { 3, val };
	_load.push_back(vec);
}

void Load::SetMx(double val) {
	std::vector<double> vec = { 4, val };
	_load.push_back(vec);
}

void Load::SetMy(double val) {
	std::vector<double> vec = { 5, val };
	_load.push_back(vec);
}

void Load::SetMz(double val) {
	std::vector<double> vec = { 6, val };
	_load.push_back(vec);
}

std::vector<std::vector<double>> Load::GetLoadVector() {
	return _load;
}

void Load::SetLoadVector(std::vector<std::vector<double>> vec) {
	_load = vec;
}

bool Load::operator ==(Load const &l2) {
	bool stat = false;
	if (_nodeID == l2._nodeID) {
		stat = true;
	}
	return stat;
}

void Load::SortByNodeID(std::vector<Load> &load) {
	std::vector<Load> load1;
	int index = 0, count = 0;

	for (int j = 0; j < load.size(); j++) { //for each support
		int min = INT32_MAX;
		for (int i = 0; i < load.size(); i++) { //for each support
			int nodeID = load[i].GetNode();
			if (count == 0) {
				if (nodeID < min) {
					min = nodeID;
					index = i;
				}
			}
			else {
				if (!(std::find(load1.begin(), load1.end(), load[i]) != load1.end())) {
					//if is not inside the sup1 vector
					if (nodeID < min) {
						min = nodeID;
						index = i;
					}
				}
			}
		}
		load1.push_back(load[index]);
		count++;
	}
	load = load1;
}

Matrix Load::AssembleLoadMatrix(const StructureManager* structManager, const PreAnalysisSetUp* setUp) {
	// Assemble the total load matrix. The vector of loads MUST be already sorted when it is passed.
	const int* DOF = setUp->DOF();
	double** complete = Matrix::CreateMatrixDouble(*setUp->StiffMatrixSize(), 1);
	

	std::map<int, Load*>::const_iterator it = structManager->Loads()->begin();		
	//for each Load instance
	while (it != structManager->Loads()->end()) {
		int node = it->second->GetNode(); //get node ID
		
		 //for each load vector
		for (int j = 0; j < it->second->GetLoadVector().size(); j++) {
			int dir = it->second->GetLoadVector()[j][0];
			double load = it->second->GetLoadVector()[j][1];
			complete[(node - 1) * (*DOF) + (dir - 1)][0] = load;
		}
		it++;
	}
	return Matrix(complete, *setUp->StiffMatrixSize(), 1);
}

Matrix Load::AssembleLoadMatrixWithFlag(std::vector<Node> &vecNode, std::vector<Load> &vecLoad, std::string flag) {
	// Assemble the total load matrix. The vector of loads MUST be already sorted when it is passed.
	int DOF = 6; //For Shell Elements
	int size = vecNode.size() * DOF;
	double** complete = Matrix::CreateMatrixDouble(size, 1);
	for (int i = 0; i < vecLoad.size(); i++) { //for each Load instance
		int node = vecLoad[i].GetNode(); //get node ID
		for (int j = 0; j < vecLoad[i].GetLoadVector().size(); j++) { //for each load vector
			int dir = vecLoad[i].GetLoadVector()[j][0];
			double load = vecLoad[i].GetLoadVector()[j][1];
			if (vecLoad[i].GetStatus() == flag) {
				complete[(node - 1) * DOF + (dir - 1)][0] = load;
			}
		}
	}
	return Matrix(complete, size, 1);
}

Matrix Load::AssembleDispLoadMatrix(std::vector<Node> &vecNode, std::vector<Support> &vecSup) {
	// Assemble the total load matrix. The vector of loads MUST be already sorted when it is passed.
	int DOF = 6; //For Shell Elements
	int size = vecNode.size() * DOF;
	double** complete = Matrix::CreateMatrixDouble(size, 1);
	for (int i = 0; i < vecSup.size(); i++) { //for each support instance
		int node = vecSup[i].GetNode(); //get node ID
		for (int j = 0; j < vecSup[i].GetSupportVector().size(); j++) { //for each sup vector
			int dir = vecSup[i].GetSupportVector()[j][0];
			double load = vecSup[i].GetSupportVector()[j][1];
			if (load != 0) {
				complete[(node - 1) * DOF + (dir - 1)][0] = load;
			}
		}
	}
	return Matrix(complete, size, 1);
}

Matrix Load::GetReducedLoadMatrix(Matrix &loadMatrix, const std::map<int, Support*>* mapSup, const int* DOF) {
	// Returns the load matrix minus the rows associated with support conditions.
	// The vector of loads MUST be already sorted when it is passed.

	int remove = 0; //keeps track of the amount of loads that were removed from the matrix
	Matrix reduced(MatrixOperation::CopyMatrixDouble(loadMatrix), loadMatrix.GetDimX(), loadMatrix.GetDimY());
	
	std::map<int, Support*>::const_iterator it = mapSup->begin();
	
	//for each Support
	while (it != mapSup->end()) {
		int node = it->second->GetNode(); //node ID

		// for each direction of support
		for (int j = 0; j < it->second->GetSupportVector().size(); j++) { 
			
			//only do this if the type of support is really a constraint (i.e., = 0), not a displacement load
			if (it->second->GetSupportVector()[j][1] == 0) { 
				int dir = it->second->GetSupportVector()[j][0];
				reduced = MatrixOperation::DeleteRow(reduced, (node - 1) * *DOF + (dir - 1) - remove);
				remove++;
			}
		}
		it++;
	}
	return reduced;
}

std::vector<int> Load::IdentifyIncrementalLoads(std::vector<Load> &vecLoad) {

	std::vector<int> indexes;
	int DOF = 6;

	for (int i = 0; i < vecLoad.size(); i++) { //for each load
		if (vecLoad[i].GetStatus() == "increment") {
			for (int j = 0; j < vecLoad[i].GetLoadVector().size(); j++) {
				indexes.emplace_back((vecLoad[i].GetNode() - 1)*DOF + (vecLoad[i].GetLoadVector()[j][0] - 1));
			}
		}
	}

	return indexes;
}

Matrix Load::MultiplyIncrementalTerms(Matrix &redLoadMatrix, std::vector<int> &incIndex, std::vector<Support> &vecSups, double mult) {

	Matrix newMatrix(redLoadMatrix.GetDimX(), 1);
	bool multiplied = false;

	for (int i = 0; i < redLoadMatrix.GetDimX(); i++) {
		for (int j = 0; j < incIndex.size(); j++) { //for each index that can be multiplied

			int index = incIndex[j] - Support::NumberOfDOFBeforeDOF(incIndex[j], vecSups);

			if (index == i) { //if the matrix term is one of the multipliable ones
				newMatrix.GetMatrixDouble()[i][0] = redLoadMatrix.GetMatrixDouble()[i][0] * mult;
				multiplied = true; 
			}
		}
		//this wil run if no load is to be multiplied OR if the current index is not to be multiplied
		if (!multiplied) {
			newMatrix.GetMatrixDouble()[i][0] = redLoadMatrix.GetMatrixDouble()[i][0];
		}
		multiplied = false;
	}
	return newMatrix;
}

Matrix Load::GetTotalForceNotOrganized(Matrix &m, Matrix &m2, std::vector<Support> &vecSup, std::vector<Node> &vecNode) {
	int DOF = 6;
	Matrix d(vecNode.size()*DOF, 1);

	for (int i = 0; i < m.GetDimX(); i++) {
		d.GetMatrixDouble()[i][0] = m.GetMatrixDouble()[i][0];
	}

	
	for (int i = 0; i < d.GetDimX() - m.GetDimX(); i++) {
		d.GetMatrixDouble()[m.GetDimX() + i][0] = m2.GetMatrixDouble()[i][0];
	}
	
	return d;
}

Matrix Load::GetTotalForceMatrix(Matrix &m, std::vector<Support> &vecSup, std::vector<Node> &vecNode, double& biggest, double loadFraction) {
	//This function returns the full displacement vector from the reduced version
	int DOF = 6; // for shell elements
	int size = vecNode.size() * DOF;

	double** disp = Matrix::CreateMatrixDouble(size, 1);

	std::vector<int> vecPos;
	std::vector<int> globalPosOfDispLoad;
	std::vector<double> dispLoad;
	std::vector<double> vecForce;
	int startOfRestricted = vecNode.size()*DOF - Support::TotalDOFsRestrained(vecSup);
	int count = 0;

	for (int i = 0; i < vecSup.size(); i++) //for each support
	{
		for (int j = 0; j < vecSup[i].GetSupportVector().size(); j++) { //this will store all the global positions and displacement values in vectors
			int nodeID = vecSup[i].GetNode();
			int dir = vecSup[i].GetSupportVector()[j][0];
			int pos = (nodeID - 1)*DOF + (dir - 1);
			if (vecSup[i].GetSupportVector()[j][1] == 0) { //only do this if the support is a boundary condition
				double force = m.GetMatrixDouble()[startOfRestricted + count][0]; //get the reaction on the support
				vecPos.emplace_back(pos);
				vecForce.emplace_back(force);
				count++;
			}
			else {
				dispLoad.push_back(vecSup[i].GetSupportVector()[j][1]); //if support is not a fixed boundary, get the applied displacement
				globalPosOfDispLoad.push_back(pos);
			}
		}
	}

	int index = 0;
	for (int i = 0; i < size; i++) { //for each term of the global D matrix
		if (std::find(vecPos.begin(), vecPos.end(), i) != vecPos.end()) {//if DOF = i a fixed boundary support type
			disp[i][0] = vecForce[index]; //add the reaction
			index++;
		}
		else {
			disp[i][0] = m.GetMatrixDouble()[i - index][0]; // get the force from the force vector
		}
	}

	
	double highStiff = biggest * pow(10, 6);
	for (int i = 0; i < globalPosOfDispLoad.size(); i++) { //remove the fictitious force applied due to displacement load
		disp[globalPosOfDispLoad[i]][0] -= highStiff * dispLoad[i] * loadFraction;
	}
	
	return Matrix(disp, size, 1);
}

double Load::DefineTotalLoadSteps(std::string type, int &specStep, int cyclicRepeat, double totalTime, double deltaT) {

	if (type == "monotonic") {
		return specStep;
	}
	else if (type == "cyclic" || type == "reverse-cyclic") {
		return specStep * cyclicRepeat;
	}
	else if (type == "dynamic") {
		return totalTime / deltaT;
	}
	else if (type == "seismic") {
		return totalTime / deltaT;
	}
	else if (type == "impulse") {
		return totalTime / deltaT;
	}
	else {
		return 0;
	}

}

double Load::DefineLoadFractionAtLoadStep(std::string type, int& step, int& totalSteps, int stepsPerPeak, double peakInc, int cyclesPerPeak, double iniPeak) {

	if (type == "monotonic") {
		return (step + 1)*(1.0 / totalSteps);
	}
	else if (type == "cyclic") {
		int stepsPerCycle = stepsPerPeak * 2;//return the total number of steps inside each full cycle
		int cycle = step / stepsPerCycle; //return the number of cycles already covered
		int aux = step + 1 - cycle * stepsPerCycle; //returns the loadstep inside the current cycle
		int aux2 = aux / stepsPerPeak; // returns the number of peaks of the current cycle it already covered

		double currentPeak = iniPeak + (cycle / cyclesPerPeak) * peakInc;
		double increment = currentPeak / stepsPerPeak;

		if (aux2 == 0) {
			return (aux) * increment;
		}
		else {
			return currentPeak - (aux - stepsPerPeak)*increment;
		}
	}
	else { //this runs if 'reverse cyclic'
		
		int stepsPerCycle = stepsPerPeak * 4;//return the total number of steps inside each full cycle
		int cycle = step / stepsPerCycle; //return the number of cycles already covered
		int aux = step + 1 - cycle * stepsPerCycle; //returns the loadstep inside the current cycle
		int aux2 = aux / stepsPerPeak; // returns the number of peaks of the current cycle it already covered

		double currentPeak = iniPeak + (cycle / cyclesPerPeak) * peakInc;
		double increment = currentPeak / stepsPerPeak;

		if (aux2 == 0) {
			return aux * increment;
		} else if (aux2 == 1) {
			return currentPeak - (aux - stepsPerPeak)*increment;
		}
		else if (aux2 == 2) {
			return -(aux - aux2 * stepsPerPeak) * increment;
		}
		else if (aux2 == 3) {
			return -currentPeak + (aux - aux2 * stepsPerPeak) * increment;
		}
		else {
			return 0;
		}
	}
}

double Load::SampleDynamicForceFunction(double amplitude, double period, double phase, double t) {
	return amplitude * sin(period * t + phase);
}

Load::~Load()
{
}