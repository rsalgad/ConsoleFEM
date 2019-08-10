#pragma once
#include "AnalysisMethod.h"
class DynamicAnalysis :
	public AnalysisMethod
{
public:
	DynamicAnalysis(); //TODO Parameters
	const int LoadSteps() const override;
	const int Iterations() const override;
	const AnalysisTypes Type() const override;

private: //TODO Parameters

};

