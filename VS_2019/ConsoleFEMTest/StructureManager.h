#pragma once
#include "pch.h"
#include "Node.h"
#include "Load.h"
#include "Support.h"
#include "ShellElement.h"
#include "Spring3D.h"
#include "Mass.h"

class StructureManager
{
public:
	StructureManager();
	~StructureManager();
	void AddNode(Node* node);
	std::vector<Node*> Nodes();
	Node* FindNodeByID(int ID);
	Node* FindNodeByCoordinates(double x, double y, double z);
	void AddShellElement(ShellElement* shell);
	std::vector<ShellElement*> ShellElements();
	void AddSpringElement(Spring3D* spring);
	std::vector<Spring3D*> SpringElements();
	void AddLoad(Load* load);
	void AddSupport(Support* sup);
	void AddMass(Mass* mass);
	void AddMaterial(MaterialModel* mat);
	MaterialModel* FindMaterialByID(int ID);

	void SortSupportsByNodeID();
	void SortLoadsByNodeID();

	void PrintNodes();
	void PrintShellElement();
	void PrintSpringElement();
	void PrintLoads();
	void PrintSupports();
	void PrintMasses();

	void CalculateShellsGlobalDOFVector();
	void CalculateSpringsGlobalDOFVector();

private:
	std::vector<Node*> _strucNodes;
	std::vector<ShellElement*> _strucShells;
	std::vector<Spring3D*> _strucSprings;
	std::vector<Load*> _strucLoads;
	std::vector<Support*> _strucSupports;
	std::vector<Mass*> _strucMasses;
	std::vector<MaterialModel*> _strucMaterials;

};

