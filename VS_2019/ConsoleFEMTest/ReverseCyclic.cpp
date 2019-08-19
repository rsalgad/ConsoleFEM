#include "pch.h"
#include "ReverseCyclic.h"

ReverseCyclic::ReverseCyclic(int nLoadSteps, int nIterations, double convLimit, int stepsPerPeak, double iniPeak, int cyclicRepeat, int cyclesPerPeak, double peakInc) : AnalysisMethod(nLoadSteps, nIterations, convLimit)
{
	_stepsPerPeak = stepsPerPeak;
	_iniPeak = iniPeak;
	_cyclicRepeat = cyclicRepeat;
	_cyclesPerPeak = cyclesPerPeak;
	_peakInc = peakInc;
}

const int* ReverseCyclic::LoadSteps() const
{
	int* ans = new int();
	*ans = _nLoadSteps * _cyclicRepeat;
	return ans;
}

const AnalysisTypes ReverseCyclic::Type() const
{
	return AnalysisTypes::Reverse_Cyclic;
}

const int* ReverseCyclic::StepsPerPeak() const
{
	return &_stepsPerPeak;
}

const int* ReverseCyclic::CyclesPerPeak() const
{
	return &_cyclesPerPeak;
}

const int* ReverseCyclic::CyclicRepeat() const
{
	return &_cyclicRepeat;
}

const double* ReverseCyclic::IniPeak() const
{
	return &_iniPeak;
}

const double* ReverseCyclic::PeakInc() const
{
	return &_peakInc;
}
