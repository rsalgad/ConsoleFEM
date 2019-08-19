#include "pch.h"
#include "MonotonicAnalysis.h"

MonotonicAnalysis::MonotonicAnalysis(int nLoadSteps, int nIterations, double convLimit) : AnalysisMethod(nLoadSteps, nIterations, convLimit)
{
}

const AnalysisTypes MonotonicAnalysis::Type() const
{
	return AnalysisTypes::Monotonic;
}
