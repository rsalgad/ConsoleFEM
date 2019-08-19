#pragma once
//#include "pch.h"
#include <vector>
#include <string>
#include "Matrix.h"


class SeismicLoad : public DynamicLoad
{
public:
	SeismicLoad();
	std::vector<double> GetTime();
	void SetTime(std::vector<double> time);
	std::vector<std::vector<double>> GetRecords();
	std::vector<char> GetDirections();
	void SetRecordX(std::vector<double> record);
	void SetRecordY(std::vector<double> record);
	void SetRecordZ(std::vector<double> record);
	double LoadFromTime(double t, char dir) const;
	int GetIndexOfDirection(char dir) const;
	std::string GetType() override;

	static Matrix GetSeismicLoadVector(const SeismicLoad& sLoad, Matrix& FInc, const double* t);
	
	~SeismicLoad();

private:
	std::vector<std::vector<double>> _records;
	std::vector<char> _directions;
	std::vector<double> _time;
};

