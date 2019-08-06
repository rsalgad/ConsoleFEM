#pragma once
#include <vector>
#include <string>

class MaterialModel
{
public:
	virtual int GetID() = 0;
	virtual std::string GetType() = 0;
	virtual ~MaterialModel();
};

