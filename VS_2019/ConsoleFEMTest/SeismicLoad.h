#pragma once
#include <vector>
#include "Matrix.h"

class SeismicLoad
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
	double LoadFromTime(double t, char dir);
	int GetIndexOfDirection(char dir);

	static Matrix GetSeismicLoadVector(SeismicLoad &sLoad, Matrix &FInc, double t);
	
	~SeismicLoad();

private:
	std::vector<std::vector<double>> _records;
	std::vector<char> _directions;
	std::vector<double> _time;
};

