#include "pch.h"
#include "StructureManager.h"
#include "Node.h"
#include "ShellElement.h"
#include "Spring3D.h"

StructureManager::StructureManager()
{
}

StructureManager::~StructureManager()
{
}

void StructureManager::AddNode(Node* node) {
	_strucNodes.push_back(node);
}

std::vector<Node*> StructureManager::Nodes()
{
	return _strucNodes;
}

Node* StructureManager::FindNodeByID(int ID) {
	//this should be optimized
	for (int i = 0; i < _strucNodes.size(); i++) {
		if (_strucNodes[i]->GetID() == ID) {
			return _strucNodes[i];
		}
	}
}

Node* StructureManager::FindNodeByCoordinates(double x, double y, double z) {
	//this should be optimized
	for (int i = 0; i < _strucNodes.size(); i++) {
		if (_strucNodes[i]->GetX() == x && _strucNodes[i]->GetY() == y && _strucNodes[i]->GetZ() == z) {
			return _strucNodes[i];
		}
	}
}

void StructureManager::AddShellElement(ShellElement* shell) {
	_strucShells.push_back(shell);
}

std::vector<ShellElement*> StructureManager::ShellElements()
{
	return _strucShells;
}

void StructureManager::AddSpringElement(Spring3D* spring) {
	_strucSprings.push_back(spring);
}

std::vector<Spring3D*> StructureManager::SpringElements()
{
	return _strucSprings;
}

void StructureManager::AddLoad(Load* load) {
	_strucLoads.push_back(load);
}

void StructureManager::AddSupport(Support* sup) {
	_strucSupports.push_back(sup);
}

void StructureManager::AddMass(Mass* mass)
{
	_strucMasses.push_back(mass);
}

void StructureManager::AddMaterial(MaterialModel* mat)
{
	_strucMaterials.push_back(mat);
}

MaterialModel* StructureManager::FindMaterialByID(int ID)
{
	for (int i = 0; i < _strucMaterials.size(); i++) {
		if (_strucMaterials[i]->GetID() == ID) {
			return _strucMaterials[i];
		}
	}
}

void StructureManager::SortSupportsByNodeID()
{
	//IF I USE A MAP, NONE OF THIS IS REQUIRED
	std::vector<Support*> sup1;
	int index = 0, count = 0;

	for (int j = 0; j < _strucSupports.size(); j++) { //for each support
		int min = INT32_MAX;
		for (int i = 0; i < _strucSupports.size(); i++) { //for each support
			int nodeID = _strucSupports[i]->GetNode();
			if (count == 0) {
				if (nodeID < min) {
					min = nodeID;
					index = i;
				}
			}
			else {
				if (!(std::find(sup1.begin(), sup1.end(), _strucSupports[i]) != sup1.end())) {
					//if is not inside the sup1 vector
					if (nodeID < min) {
						min = nodeID;
						index = i;
					}
				}
			}
		}
		sup1.push_back(_strucSupports[index]);
		count++;
	}

	//Need to destroy the previous _strucSupports first before assigning it, otherwise memory leak
	_strucSupports = sup1; 
}

void StructureManager::SortLoadsByNodeID()
{
	//IF I USE A MAP, NONE OF THIS IS REQUIRED
	std::vector<Load*> load1;
	int index = 0, count = 0;

	for (int j = 0; j < _strucLoads.size(); j++) { //for each support
		int min = INT32_MAX;
		for (int i = 0; i < _strucLoads.size(); i++) { //for each support
			int nodeID = _strucLoads[i]->GetNode();
			if (count == 0) {
				if (nodeID < min) {
					min = nodeID;
					index = i;
				}
			}
			else {
				if (!(std::find(load1.begin(), load1.end(), _strucLoads[i]) != load1.end())) {
					//if is not inside the sup1 vector
					if (nodeID < min) {
						min = nodeID;
						index = i;
					}
				}
			}
		}
		load1.push_back(_strucLoads[index]);
		count++;
	}

	//Need to destroy the previous _strucSupports first before assigning it, otherwise memory leak
	_strucLoads = load1;
}

void StructureManager::PrintNodes() {
	for (int i = 0; i < _strucNodes.size(); i++) {
		std::cout << _strucNodes[i]->ToString() << std::endl;
	}
}

void StructureManager::PrintShellElement() {
	for (int i = 0; i < _strucShells.size(); i++) {
		std::cout << _strucShells[i]->ToString() << std::endl;
	}
}

void StructureManager::PrintSpringElement() {
	for (int i = 0; i < _strucSprings.size(); i++) {
		std::cout << _strucSprings[i]->ToString() << std::endl;
	}
}

void StructureManager::PrintLoads() {
	for (int i = 0; i < _strucLoads.size(); i++) {
		std::cout << _strucLoads[i]->ToString() << std::endl;
	}
}

void StructureManager::PrintSupports() {
	for (int i = 0; i < _strucSupports.size(); i++) {
		std::cout << _strucSupports[i]->ToString() << std::endl;
	}
}

void StructureManager::PrintMasses()
{
	for (int i = 0; i < _strucMasses.size(); i++) {
		std::cout << _strucMasses[i]->ToString() << std::endl;
	}
}

void StructureManager::CalculateShellsGlobalDOFVector() {
	if (_strucShells.size() != 0) {
		for (int i = 0; i < _strucShells.size(); i++) {
			_strucShells[i]->CalculateGlobalDOFVector(_strucSupports);
			_strucShells[i]->CalculateGlobalMassDOFVector(_strucMasses, _strucSupports);
		}
	}
}

void StructureManager::CalculateSpringsGlobalDOFVector() {
	if (_strucSprings.size() != 0) {
		for (int i = 0; i < _strucSprings.size(); i++) {
			_strucSprings[i]->CalculateGlobalDOFVector(_strucSupports);
		}
	}
}
