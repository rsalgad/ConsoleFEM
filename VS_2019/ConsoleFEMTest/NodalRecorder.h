#pragma once
//#include "pch.h"

#include <map>
#include <vector>
#include "Node.h"

template <typename T>
class NodalRecorder
{
private:
	std::vector<int> _nodesID;
	std::map<int, std::map<int, T*>> _records;
public:
	NodalRecorder();
	NodalRecorder(std::vector<int> nodesID);
	NodalRecorder(const std::map<int, Node*>* allNodes);
	void Add(std::map<int, T*> vec);
	const std::map<int, std::map<int, T*>>* GetRecord() const;
	~NodalRecorder();
};

