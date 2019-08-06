#pragma once
#include "pch.h"
#include "Node.h"
#include "ShellElement.h"
#include "Spring3D.h"

class StructureManager
{
public:
	StructureManager();
	~StructureManager();
	void AddNode(Node* node);
	std::vector<Node*> Nodes();
	Node* FindNodeByID(int& ID);
	Node* FindNodeByCoordinates(double x, double y, double z);
	void AddShellElement(ShellElement* shell);
	std::vector<ShellElement*> ShellElements();
	void AddSpringElement(Spring3D* spring);
	std::vector<Spring3D*> SpringElements();
	void PrintNodes();
	void PrintShellElement();
	void PrintSpringElement();

private:
	std::vector<Node*> _strucNodes;
	std::vector<ShellElement*> _strucShells;
	std::vector<Spring3D*> _strucSprings;
};

