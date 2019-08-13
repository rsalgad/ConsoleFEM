#pragma once
#include "AnalysisMethod.h"
class DynamicAnalysis :
	public AnalysisMethod
{
public:
	DynamicAnalysis(int totalTime, int deltaT); //TODO Parameters
	const int* LoadSteps() const override;
	const int* Iterations() const override;
	const double* TotalTime() const;
	const double* DeltaT() const;
	const AnalysisTypes Type() const override;

private: //TODO Parameters
	double _totalTime;
	double _deltaT;
};

