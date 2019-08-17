#pragma once
//#include "pch.h"
#include <string>

class Element
{
public:
	virtual int GetID() const = 0;
	virtual std::string ToString() const = 0;
	virtual ~Element();
};

