#pragma once
#include "AnalysisMethod.h"
class ElasticAnalysis :
	public AnalysisMethod
{
public:
	ElasticAnalysis (int nLoadStep, int nIterations);
	const AnalysisTypes Type() const override;
};

