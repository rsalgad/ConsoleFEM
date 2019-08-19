#pragma once
//#include "pch.h"
#include "DynamicLoad.h"
#include <vector>
#include <string>

class ImpulseLoad : public DynamicLoad
{
public:
	ImpulseLoad();
	std::vector<std::vector<double>> GetPoints() const;
	void SetPoints(std::vector<std::vector<double>> points);
	std::string GetType() override;

	static double LoadFromTime(const ImpulseLoad* impLoad, const double* t);
	~ImpulseLoad();

private:
	std::vector<std::vector<double>> _points;
};

