#pragma once
#include "AnalysisMethod.h"
class MonotonicAnalysis :
	public AnalysisMethod
{
	MonotonicAnalysis(int nLoadSteps, int nIterations);
	const AnalysisTypes Type() const override;
};

