#pragma once

enum class IntegrationMethod {
	AverageNewmark, LinearNewmark, WilsonTheta, HHTalpha
};

class TimeIntegrationMethod
{
	//Wilson-Theta: Gama = 1/2, Beta = 1/6, Theta = 1.42, alpha = 0; HTT-alpha: Gama = 0.6, Beta 0.3025, theta = 1, alpha = 0.1
public:
	TimeIntegrationMethod();
	TimeIntegrationMethod(IntegrationMethod method);
	double GetNewmarkGama() const;
	double GetNewmarkBeta() const;;
	double GetWilsonTheta() const;;
	double GetAlphaF() const;;
	double GetAlphaM() const;;

private:
	double _newmarkGama;
	double _newmarkBeta;
	double _wilsonTheta;
	double _alphaF;
	double _alphaM;
};



