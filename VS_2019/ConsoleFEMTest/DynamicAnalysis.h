#pragma once
//#include "pch.h"
#include "AnalysisMethod.h"
#include "TimeIntegrationMethod.h"
#include "DynamicLoad.h"

class DynamicAnalysis :
	public AnalysisMethod
{
public:
	DynamicAnalysis(double totalTime, double deltaT, AnalysisTypes type, TimeIntegrationMethod intMethod, DynamicLoad* dynLoad, int iterations, double convLimit); //TODO Parameters
	DynamicAnalysis(double totalTime, double deltaT, AnalysisTypes type, TimeIntegrationMethod intMethod, DynamicLoad* dynLoad, int dampMode1, int dampMode2, double damp1, double damp2, int iterations, double convLimit);
	const int* LoadSteps() const override;
	const double* TotalTime() const;
	const double* DeltaT() const;
	const int* DampMode1() const;
	const int* DampMode2() const;
	const double* Damp1() const;
	const double* Damp2() const;
	DynamicLoad* Load() const;
	const TimeIntegrationMethod* IntegrationMethod() const;
	const AnalysisTypes Type() const override;

private: //TODO Parameters
	double CalculateTotalTime();
	DynamicLoad* _dynLoad;
	AnalysisTypes _type;
	double _totalTime;
	double _deltaT;
	int _dampMode1 = 1, _dampMode2 = 2;
	double _damp1 = 0, _damp2 = 0;
	TimeIntegrationMethod _intMethod;
};

