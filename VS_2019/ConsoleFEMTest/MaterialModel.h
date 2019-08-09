#pragma once

class MaterialModel
{
public:
	virtual int GetID() = 0;
	virtual std::string GetType() = 0;
	virtual std::string ToString() = 0;
	virtual ~MaterialModel();
};

