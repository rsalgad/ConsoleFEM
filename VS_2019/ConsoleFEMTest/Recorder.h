#pragma once
#include <map>

template <typename T>
class Recorder
{
private:
	std::vector<int> _nodesID;
	std::map<int, std::map<int, T>> _records;
public:
	Recorder();
	Recorder(std::vector<int> nodesID);
	Recorder(char c, std::map<int, Node*> allNodes);
	void Add(std::map<int, T> vec);
	const std::map<int, std::map<int, T>> GetRecord() const;
	~Recorder();
};

