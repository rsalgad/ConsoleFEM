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

Node* StructureManager::FindNodeByID(int& ID) {
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
