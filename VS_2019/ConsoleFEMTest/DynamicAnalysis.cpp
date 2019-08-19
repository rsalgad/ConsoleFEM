#include "pch.h"

DynamicAnalysis::DynamicAnalysis(double totalTime, double deltaT, AnalysisTypes type, TimeIntegrationMethod intMethod, DynamicLoad* dynLoad, int iterations, double convLimit) : AnalysisMethod(-1, iterations)
{
	_totalTime = totalTime;
	_deltaT = deltaT;
	_type = type;
	_dynLoad = dynLoad;
	_intMethod = intMethod;
	_convLimit = convLimit;
}

DynamicAnalysis::DynamicAnalysis(double totalTime, double deltaT, AnalysisTypes type, TimeIntegrationMethod intMethod, DynamicLoad* dynLoad, int dampMode1, int dampMode2, double damp1, double damp2, int iterations, double convLimit) : AnalysisMethod(-1, iterations)
{
	_totalTime = totalTime;
	_deltaT = deltaT;
	_type = type;
	_dynLoad = dynLoad;
	_intMethod = intMethod;
	_dampMode1 = dampMode1;
	_dampMode2 = dampMode2;
	_damp1 = damp1;
	_damp2 = damp2;
	_convLimit = convLimit;
}


const int* DynamicAnalysis::LoadSteps() const
{
	int* loadSteps = new int;
	*loadSteps = _totalTime / _deltaT;
	return loadSteps;
}

const double* DynamicAnalysis::TotalTime() const
{
	return &_totalTime;
}

const double* DynamicAnalysis::DeltaT() const
{
	return &_deltaT;
}

const int* DynamicAnalysis::DampMode1() const
{
	return &_dampMode1;
}

const int* DynamicAnalysis::DampMode2() const
{
	return &_dampMode2;
}

const double* DynamicAnalysis::Damp1() const
{
	return &_damp1;
}

const double* DynamicAnalysis::Damp2() const
{
	return &_damp2;
}

DynamicLoad* DynamicAnalysis::Load() const
{
	return _dynLoad;
}

const TimeIntegrationMethod* DynamicAnalysis::IntegrationMethod() const
{
	return &_intMethod;
}

const AnalysisTypes DynamicAnalysis::Type() const
{
	return _type;
}

double DynamicAnalysis::CalculateTotalTime()
{
	switch (_type)
	{
	case AnalysisTypes::Seismic:
		{
			SeismicLoad* seisLoad = static_cast<SeismicLoad*>(_dynLoad);
			return seisLoad->GetTime()[0] * seisLoad->GetRecords()[0].size();
		}
		break;
	case AnalysisTypes::Impulse:
		{
			ImpulseLoad* impLoad = static_cast<ImpulseLoad*>(_dynLoad);
			return impLoad->GetPoints()[impLoad->GetPoints().size() - 1][0] + 1.93; //1.93 needs to be changed
		}
		break;
	default:
		{
		}
		break;
	}
	return 0.0;
}
