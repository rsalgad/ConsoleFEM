#include "pch.h"
#include "ElasticAnalysis.h"

ElasticAnalysis::ElasticAnalysis(int nLoadStep, int nIterations) : AnalysisMethod(nLoadStep, nIterations)
{
}

const AnalysisTypes ElasticAnalysis::Type() const
{
	return AnalysisTypes::Elastic;
}
