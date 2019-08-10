#pragma once
#include "AnalysisMethod.h"

class CyclicAnalysis :
	public AnalysisMethod
{
public:
	CyclicAnalysis(); //TODO Parameters
	const int LoadSteps() const override;
	const int Iterations() const override;
	const AnalysisTypes Type() const override;

private: //TODO Additional Parameters
};

