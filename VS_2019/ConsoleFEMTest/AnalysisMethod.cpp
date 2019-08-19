#include "pch.h"

AnalysisMethod::AnalysisMethod(int nLoadSteps, int nIterations)
{
	_nLoadSteps = nLoadSteps;
	_niterations = nIterations;
	_convLimit = -1;
}

AnalysisMethod::AnalysisMethod(int nLoadSteps, int nIterations, double convLimit)
{
	_nLoadSteps = nLoadSteps;
	_niterations = nIterations;
	_convLimit = convLimit;
}

AnalysisMethod::~AnalysisMethod()
{
}

const int* AnalysisMethod::LoadSteps() const
{
	return &_nLoadSteps;
}

const int* AnalysisMethod::Iterations() const
{
	return &_niterations;
}

const double* AnalysisMethod::ConvergenceLimit() const
{
	return &_convLimit;
}
