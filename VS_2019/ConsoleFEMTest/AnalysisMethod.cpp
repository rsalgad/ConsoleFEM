#include "pch.h"
#include "AnalysisMethod.h"

AnalysisMethod::AnalysisMethod()
{
}

AnalysisMethod::AnalysisMethod(int nLoadSteps, int nIterations)
{
	_nLoadSteps = nLoadSteps;
	_niterations = nIterations;
}

AnalysisMethod::~AnalysisMethod()
{
}

const int AnalysisMethod::LoadSteps() const
{
	return _nLoadSteps;
}

const int AnalysisMethod::Iterations() const
{
	return _niterations;
}
