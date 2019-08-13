#include "pch.h"
#include "DynamicAnalysis.h"

DynamicAnalysis::DynamicAnalysis(int totalTime, int deltaT)
{
	_totalTime = totalTime;
	_deltaT = deltaT;
}

const int* DynamicAnalysis::LoadSteps() const
{
	int* loadSteps = new int;
	*loadSteps = _totalTime / _deltaT;
	return loadSteps;
}

const int* DynamicAnalysis::Iterations() const
{
	return 0;
}

const double* DynamicAnalysis::TotalTime() const
{
	return &_totalTime;
}

const double* DynamicAnalysis::DeltaT() const
{
	return &_deltaT;
}
