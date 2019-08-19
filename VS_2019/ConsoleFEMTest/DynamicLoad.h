#pragma once
//#include "pch.h"
#include <string>

class DynamicLoad
{
public:
	DynamicLoad();
	virtual std::string GetType() = 0;
};

