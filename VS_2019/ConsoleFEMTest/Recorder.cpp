#include "pch.h"
#include "Recorder.h"

template<typename T>
Recorder<T>::~Recorder()
{
}

template<typename T>
Recorder<T>::Recorder()
{
}

template<typename T>
Recorder<T>::Recorder(std::vector<int> nodesID)
{
	_nodesID = nodesID;
}

template<typename T>
Recorder<T>::Recorder(char c, std::map<int, Node*> allNodes)
{
	if (c == 'a') {
		std::map<int, Node*>::iterator it;
		_nodesID.reserve(allNodes.size());
		for (int i = 0; i < allNodes.size(); i++) {
			_nodesID.push_back(i + 1);
		}
	}
	else {
		std::cout << "Error initializing recorder" << std::endl;
	}
}

template<typename T>
void Recorder<T>::Add(std::map<int, T> map)
{
	std::map<int, T> values;

	//finds all the nodes that were specified in the recorder
	for (int i = 0; i < _nodesID.size(); i++) {
		std::map<int, T>::iterator it = map.find(_nodesID[i]);
		if (it != map.end()) {
			values.insert(std::pair<int,T>(it->first, it->second));
		}
	}
	_records.insert(std::pair<int, std::map<int, T>>(_records.size() + 1, values));
}

template<typename T>
const std::map<int, std::map<int, T>> Recorder<T>::GetRecord() const
{
	return _records;
}


