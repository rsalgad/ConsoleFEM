#pragma once

enum class IntegrationMethod {
	AverageNewmark, LinearNewmark, WilsonTheta, HHTalpha
};

class TimeIntegrationMethod
{
public:
	TimeIntegrationMethod(IntegrationMethod method);
	double GetNewmarkGama();
	double GetNewmarkBeta();
	double GetWilsonTheta();
	double GetAlphaF();
	double GetAlphaM();

private:
	TimeIntegrationMethod();
	double _newmarkGama;
	double _newmarkBeta;
	double _wilsonTheta;
	double _alphaF;
	double _alphaM;
};



