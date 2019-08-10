#pragma once
#include "AnalysisMethod.h"
class ReverseCyclic :
	public AnalysisMethod
{
public:
	ReverseCyclic(); ; //TODO Parameters
	const int LoadSteps() const override;
	const int Iterations() const override;
	const AnalysisTypes Type() const override;

private: //TODO Additional Parameters
};

