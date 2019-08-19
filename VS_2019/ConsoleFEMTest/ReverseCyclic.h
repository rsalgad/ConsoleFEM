#pragma once
//#include "pch.h"
#include "AnalysisMethod.h"

class ReverseCyclic :
	public AnalysisMethod
{
public:
	ReverseCyclic(int nLoadSteps, int nIterations, double convLimit, int stepsPerPeak, double iniPeak, int cyclicRepeat, int cyclesPerPeak, double peakInc); //TODO Parameters
	const int* LoadSteps() const override;
	const AnalysisTypes Type() const override;
	const int* StepsPerPeak() const;
	const int* CyclesPerPeak() const;
	const int* CyclicRepeat() const;
	const double* IniPeak() const;
	const double* PeakInc() const;

private: //TODO Additional Parameters
	int _stepsPerPeak, _cyclesPerPeak, _cyclicRepeat;
	double _iniPeak, _peakInc;
};

