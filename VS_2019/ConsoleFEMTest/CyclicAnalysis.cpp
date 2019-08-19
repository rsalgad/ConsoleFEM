#include "pch.h"
#include "CyclicAnalysis.h"

CyclicAnalysis::CyclicAnalysis(int nLoadSteps, int nIterations, double convLimit, int stepsPerPeak, double iniPeak, int cyclicRepeat, int cyclesPerPeak, double peakInc) : AnalysisMethod(nLoadSteps, nIterations, convLimit)
{
	_stepsPerPeak = stepsPerPeak;
	_iniPeak = iniPeak;
	_cyclicRepeat = cyclicRepeat;
	_cyclesPerPeak = cyclesPerPeak;
	_peakInc = peakInc;
}

const int* CyclicAnalysis::LoadSteps() const
{
	int* ans = new int();
	*ans = _nLoadSteps * _cyclicRepeat;
	return ans;
}

const AnalysisTypes CyclicAnalysis::Type() const
{
	return AnalysisTypes::Cyclic;
}

const int* CyclicAnalysis::StepsPerPeak() const
{
	return &_stepsPerPeak;
}

const int* CyclicAnalysis::CyclesPerPeak() const
{
	return &_cyclesPerPeak;
}

const int* CyclicAnalysis::CyclicRepeat() const
{
	return &_cyclicRepeat;
}

const double* CyclicAnalysis::IniPeak() const
{
	return &_iniPeak;
}

const double* CyclicAnalysis::PeakInc() const
{
	return &_peakInc;
}
