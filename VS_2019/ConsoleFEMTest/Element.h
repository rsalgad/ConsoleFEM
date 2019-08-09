#pragma once

class Element
{
public:
	virtual int GetID() = 0;
	virtual std::string ToString() = 0;
	virtual ~Element();
};

