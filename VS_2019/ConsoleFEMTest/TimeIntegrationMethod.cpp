#include "pch.h"

TimeIntegrationMethod::TimeIntegrationMethod()
{
}

TimeIntegrationMethod::TimeIntegrationMethod(IntegrationMethod method)
{
	switch (method) {
	case IntegrationMethod::AverageNewmark:
		_newmarkGama = 1.0 / 2;
		_newmarkBeta = 1.0 / 4;
		_wilsonTheta = 1.0;
		_alphaF = 0;
		_alphaM = 0;
		break;
	case IntegrationMethod::LinearNewmark:
		_newmarkGama = 1.0 / 2;
		_newmarkBeta = 1.0 / 6;
		_wilsonTheta = 1.0;
		_alphaF = 0;
		_alphaM = 0;
		break;
	case IntegrationMethod::WilsonTheta:
		_newmarkGama = 1.0 / 2;
		_newmarkBeta = 1.0 / 6;
		_wilsonTheta = 1.42;
		_alphaF = 0;
		_alphaM = 0;
		break;
	case IntegrationMethod::HHTalpha:
		_newmarkGama = 0.6;
		_newmarkBeta = 0.3025;
		_wilsonTheta = 1.0;
		_alphaF = 0.1;
		_alphaM = 0;
		break;
	}
}

double TimeIntegrationMethod::GetNewmarkGama() const
{
	return _newmarkGama;
}

double TimeIntegrationMethod::GetNewmarkBeta() const
{
	return _newmarkBeta;
}

double TimeIntegrationMethod::GetWilsonTheta() const
{
	return _wilsonTheta;
}

double TimeIntegrationMethod::GetAlphaF() const
{
	return _alphaF;
}

double TimeIntegrationMethod::GetAlphaM() const
{
	return _alphaM;
}


