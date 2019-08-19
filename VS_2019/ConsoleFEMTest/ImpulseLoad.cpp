#include "pch.h"


ImpulseLoad::ImpulseLoad()
{

}

std::vector<std::vector<double>> ImpulseLoad::GetPoints() const
{
	return _points;
}

void ImpulseLoad::SetPoints(std::vector<std::vector<double>> points)
{
	_points = points;
}

std::string ImpulseLoad::GetType()
{
	return std::string();
}

double ImpulseLoad::LoadFromTime(const ImpulseLoad* impLoad,const double* t) {
	
	if (*t <= impLoad->GetPoints()[0][0]) {
		return impLoad->GetPoints()[0][1] * (*t) / impLoad->GetPoints()[0][0];
	}
	
	for (int i = 0; i < impLoad->GetPoints().size() - 1; i++) {
		if (*t > impLoad->GetPoints()[i][0] && *t <= impLoad->GetPoints()[i + 1][0]) {
			return (impLoad->GetPoints()[i][1] * (impLoad->GetPoints()[i + 1][0] - *t) + impLoad->GetPoints()[i + 1][1] * (*t - impLoad->GetPoints()[i][0])) / (impLoad->GetPoints()[i + 1][0] - impLoad->GetPoints()[i][0]);
		}
	}

	return 0;
	
	/*
	if (t <= impLoad.GetPoints()[0][0]) {
		return impLoad.GetPoints()[0][1] * t / impLoad.GetPoints()[0][0];
	}
	else if (t > impLoad._points[0][0] && t <= impLoad._points[1][0]) {
		return (impLoad.GetPoints()[0][1] * (impLoad.GetPoints()[1][0] - t) + impLoad.GetPoints()[1][1] *(t - impLoad.GetPoints()[0][0])) / (impLoad.GetPoints()[1][0] - impLoad.GetPoints()[0][0]);
	}
	else {
		if (t <= impLoad._points[2][0]) {
			return (impLoad.GetPoints()[1][1] * (impLoad.GetPoints()[2][0] - t) + impLoad.GetPoints()[2][1] * (t - impLoad.GetPoints()[1][0])) / (impLoad.GetPoints()[2][0] - impLoad.GetPoints()[1][0]);
		}
		else {
			return 0;
		}
	}
	*/
}

ImpulseLoad::~ImpulseLoad()
{
}
