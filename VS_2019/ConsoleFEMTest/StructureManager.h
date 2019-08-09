#pragma once
#include "pch.h"
#include "Node.h"
#include "Load.h"
#include "Support.h"
#include "ShellElement.h"
#include "Spring3D.h"
#include "Mass.h"
#include <map>

class StructureManager
{
public:
	StructureManager();
	~StructureManager();
	void AddNode(Node* node);
	const std::map<int, Node*>* Nodes() const;
	Node* FindNodeByID(int ID);
	Node* FindNodeByCoordinates(double x, double y, double z);
	void AddShellElement(ShellElement* shell);
	const std::map<int, ShellElement*>* ShellElements() const;
	void AddSpringElement(Spring3D* spring);
	const std::map<int, Spring3D*>* SpringElements() const;
	void AddLoad(Load* load);
	const std::map<int, Load*>* Loads() const;
	void AddSupport(Support* sup);
	const std::map<int, Support*>* Supports() const;
	void AddMass(Mass* mass);
	const std::map<int, Mass*>* Masses() const;
	void AddMaterial(MaterialModel* mat);
	const std::map<int, MaterialModel*>* Materials() const;
	MaterialModel* FindMaterialByID(int ID);

	void PrintNodes();
	void PrintShellElement();
	void PrintSpringElement();
	void PrintLoads();
	void PrintSupports();
	void PrintMasses();

	void CalculateShellsGlobalDOFVector();
	void CalculateSpringsGlobalDOFVector();

private:
	std::map<int, Node*> _strucNodes;
	std::map<int, ShellElement*> _strucShells;
	std::map<int, Spring3D*> _strucSprings;
	std::map<int, Load*> _strucLoads;
	std::map<int, Support*> _strucSupports;
	std::map<int, Mass*> _strucMasses;
	std::map<int, MaterialModel*> _strucMaterials;

};

