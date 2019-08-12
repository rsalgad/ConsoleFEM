#pragma once
#include <map>

template <typename T>
class NodalRecorder
{
private:
	std::vector<int> _nodesID;
	std::map<int, std::map<int, T>> _records;
public:
	NodalRecorder();
	NodalRecorder(std::vector<int> nodesID);
	NodalRecorder(char c, std::map<int, Node*>* allNodes);
	void Add(std::map<int, T> vec);
	const std::map<int, std::map<int, T>> GetRecord() const;
	~NodalRecorder();
};

