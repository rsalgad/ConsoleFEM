#pragma once
//#include "pch.h"
#include "AnalysisMethod.h"

class MonotonicAnalysis :
	public AnalysisMethod
{
public:
	MonotonicAnalysis(int nLoadSteps, int nIterations, double convLimit);
	const AnalysisTypes Type() const override;
};

