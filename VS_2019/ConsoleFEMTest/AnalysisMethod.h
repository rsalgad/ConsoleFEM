#pragma once

enum class AnalysisTypes {
	Elastic, Monotonic, Cyclic, Reverse_Cyclic, Seismic, Impulse
};

class AnalysisMethod
{
public:
	AnalysisMethod();
	AnalysisMethod(int nLoadSteps, int nIterations);
	virtual ~AnalysisMethod();

	virtual const int* LoadSteps() const;
	virtual const int* Iterations() const;
	virtual const AnalysisTypes Type() const = 0;

protected:
	int _nLoadSteps;
	int _niterations;

};

