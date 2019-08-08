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
	_strucNodes.insert(std::make_pair(node->GetID(), node));
}

const std::map<int, Node*> StructureManager::Nodes()
{
	return _strucNodes;
}

Node* StructureManager::FindNodeByID(int ID) {
	std::map<int, Node*>::iterator it = _strucNodes.find(ID);
	if (it != _strucNodes.end()) {
		return it->second;
	}
	else {
		return nullptr;
	}
}

Node* StructureManager::FindNodeByCoordinates(double x, double y, double z) {
	std::map<int, Node*>::iterator it = _strucNodes.begin();
	while (it != _strucNodes.end())
	{
		if (it->second->GetX() == x && it->second->GetY() == y && it->second->GetZ() == z)
			return it->second;
		it++;
	}
	return nullptr;
}

void StructureManager::AddShellElement(ShellElement* shell) {
	_strucShells.insert(std::make_pair(shell->GetID(), shell));
}

const std::map<int, ShellElement*> StructureManager::ShellElements()
{
	return _strucShells;
}

void StructureManager::AddSpringElement(Spring3D* spring) {
	_strucSprings.insert(std::make_pair(spring->GetID(), spring));
}

const std::map<int, Spring3D*> StructureManager::SpringElements()
{
	return _strucSprings;
}

void StructureManager::AddLoad(Load* load) {
	_strucLoads.insert(std::make_pair(load->GetID(), load));
}

const std::map<int, Load*> StructureManager::Loads()
{
	return _strucLoads;
}

void StructureManager::AddSupport(Support* sup) {
	_strucSupports.insert(std::make_pair(sup->GetID(), sup));
}

const std::map<int, Support*> StructureManager::Supports()
{
	return _strucSupports;
}

void StructureManager::AddMass(Mass* mass)
{
	_strucMasses.insert(std::make_pair(mass->GetID(), mass));
}

const std::map<int, Mass*> StructureManager::Masses()
{
	return _strucMasses;
}

void StructureManager::AddMaterial(MaterialModel* mat)
{
	_strucMaterials.insert(std::make_pair(mat->GetID(), mat));
}

const std::map<int, MaterialModel*> StructureManager::Materials()
{
	return _strucMaterials;
}

MaterialModel* StructureManager::FindMaterialByID(int ID)
{
	std::map<int, MaterialModel*>::iterator it = _strucMaterials.find(ID);
	if (it != _strucMaterials.end()) {
		return it->second;
	}
	else {
		return nullptr;
	}
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
