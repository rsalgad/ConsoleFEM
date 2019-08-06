#pragma once
#include <vector>

class ImpulseLoad
{
public:
	ImpulseLoad();
	std::vector<std::vector<double>> GetPoints();
	void SetPoints(std::vector<std::vector<double>> points);
	
	static double LoadFromTime(ImpulseLoad &impLoad, double t);
	~ImpulseLoad();

private:

	std::vector<std::vector<double>> _points;
};

