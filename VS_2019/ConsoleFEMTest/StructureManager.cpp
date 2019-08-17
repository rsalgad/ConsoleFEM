#include "pch.h"


StructureManager::StructureManager()
{
}

StructureManager::~StructureManager()
{
}

void StructureManager::AddNode(Node* node) {
	_strucNodes.insert(std::make_pair(*node->GetID(), node));
}

//the first 'const' means the returned map is constant (can't be changed)
//the last 'const' means everything that will be done in this function
//do not change the values of any class parameters.
const std::map<int, Node*>* StructureManager::Nodes() const
{
	return &_strucNodes;
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
		if (*it->second->GetX() == x && *it->second->GetY() == y && *it->second->GetZ() == z)
			return it->second;
		it++;
	}
	return nullptr;
}

void StructureManager::AddShellElement(ShellElement* shell) {
	_strucShells.insert(std::make_pair(shell->GetID(), shell));
}

const std::map<int, ShellElement*>* StructureManager::ShellElements() const
{
	return &_strucShells;
}

void StructureManager::AddSpringElement(Spring3D* spring) {
	_strucSprings.insert(std::make_pair(spring->GetID(), spring));
}

const std::map<int, Spring3D*>* StructureManager::SpringElements() const
{
	return &_strucSprings;
}

void StructureManager::AddLoad(Load* load) {
	_strucLoads.insert(std::make_pair(load->GetID(), load));
}

const std::map<int, Load*>* StructureManager::Loads() const
{
	return &_strucLoads;
}

void StructureManager::AddSupport(Support* sup) {
	_strucSupports.insert(std::make_pair(sup->GetID(), sup));
}

const std::map<int, Support*>* StructureManager::Supports() const
{
	return &_strucSupports;
}

void StructureManager::AddMass(Mass* mass)
{
	_strucMasses.insert(std::make_pair(mass->GetID(), mass));
}

const std::map<int, Mass*>* StructureManager::Masses() const
{
	return &_strucMasses;
}

void StructureManager::AddMaterial(MaterialModel* mat)
{
	_strucMaterials.insert(std::make_pair(mat->GetID(), mat));
}

const std::map<int, MaterialModel*>* StructureManager::Materials() const
{
	return &_strucMaterials;
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

Support* StructureManager::FindSupportByNodeID(int nodeID) const {
	std::map<int, Support*>::const_iterator it = _strucSupports.begin();
	while (it != _strucSupports.end()) {
		if (it->second->GetNode() == nodeID) {
			return it->second;
		}
		it++;
	}
	return nullptr;
}

void StructureManager::SortSupportsByNodeID()
{
	std::map<int, Support*>::const_iterator it = _strucSupports.begin();
	std::map<int, Support*> newMap;
	std::vector<int> IDs;
	IDs.reserve(_strucSupports.size());
	while (it != _strucSupports.end()) {
		IDs.push_back(it->second->GetNode());
		it++;
	}

	std::sort(IDs.begin(), IDs.end());

	for (int i = 0; i < IDs.size(); i++) {
		newMap.insert(std::pair<int, Support*>(i+1, FindSupportByNodeID(IDs[i])));
	}

	_strucSupports = newMap;

}
