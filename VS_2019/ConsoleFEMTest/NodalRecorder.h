#pragma once
#include <map>
#include <vector>
#include "Node.h"

template <typename T>
class NodalRecorder
{
private:
	std::vector<int> _nodesID;
	std::map<int, std::map<int, T>> _records;
public:
	NodalRecorder(){}

	NodalRecorder(std::vector<int> nodesID) {
		_nodesID = nodesID;
	}

	NodalRecorder(const std::map<int, Node*>* allNodes) {

		_nodesID.reserve(allNodes->size());
		for (int i = 0; i < allNodes->size(); i++) {
			_nodesID.push_back(i + 1);
		}

	}

	void Add(std::map<int, T> map)
	{
		std::map<int, T> values;

		//finds all the nodes that were specified in the recorder
		for (int i = 0; i < _nodesID.size(); i++) {
			std::map<int, T>::iterator it = map.find(_nodesID[i]);
			if (it != map.end()) {
				values.insert(std::pair<int, T>(it->first, it->second));
			}
		}
		_records.insert(std::pair<int, std::map<int, T>>(_records.size() + 1, values));
	}

	const std::map<int, std::map<int, T>>* GetRecord() const {
		return &_records;
	}

	const std::vector<int>* Nodes() const {
		return &_nodesID;
	}

	~NodalRecorder() {}
};

