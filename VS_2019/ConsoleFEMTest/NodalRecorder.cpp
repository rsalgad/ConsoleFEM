#include "pch.h"

template<typename T>
NodalRecorder<T>::~NodalRecorder()
{
}

template<typename T>
NodalRecorder<T>::NodalRecorder()
{
}

template<typename T>
NodalRecorder<T>::NodalRecorder(std::vector<int> nodesID)
{
	_nodesID = nodesID;
}

template<typename T>
NodalRecorder<T>::NodalRecorder(const std::map<int, Node*>* allNodes)
{
	std::map<int, Node*>::iterator it;
	_nodesID.reserve(allNodes->size());
	for (int i = 0; i < allNodes->size(); i++) {
		_nodesID.push_back(i + 1);
	}
}

template<typename T>
void NodalRecorder<T>::Add(std::map<int, T*> map)
{
	std::map<int, T*> values;

	//finds all the nodes that were specified in the recorder
	for (int i = 0; i < _nodesID.size(); i++) {
		std::map<int, T*>::iterator it = map.find(_nodesID[i]);
		if (it != map.end()) {
			values.insert(std::pair<int,T*>(it->first, it->second));
		}
	}
	_records.insert(std::pair<int, std::map<int, T*>>(_records.size() + 1, values));
}

template<typename T>
const std::map<int, std::map<int, T*>>* NodalRecorder<T>::GetRecord() const
{
	return &_records;
}


