#include "pch.h"
#include "SpringGeneralModel.h"
#include <math.h>

// This Material Model is said to be "General" because it is assumed that both shear and out-of-plane will have similar shapes and parameters
SpringGeneralModel::SpringGeneralModel(int ID, double iniStiff, double dMax, double fMax, double degStiff, double fRes, double dUlt, double unlStiff, double fUnl, double conStiff, double relStiff)
{
	_ID = ID;
	_iniStiff = iniStiff;
	_dMax = dMax;
	_fMax = fMax;
	_degStiff = degStiff;
	_fRes = fRes;
	_dUlt = dUlt;
	_unlStiff = unlStiff;
	_fUnl = fUnl;
	_conStiff = conStiff;
	_relStiff = relStiff;
}

int SpringGeneralModel::GetID() {
	return _ID;
}

double SpringGeneralModel::GetSecantStiffnessFromDisplacement(double disp, double plasticDisp, double maxD, double minD, std::string prevStage, std::string stage, double minUnlDisp, double maxRelDisp) {
	double n = GetNValue();
	double realDisp = fabs(disp - plasticDisp);

	//std::string stage = GetLoadingStage(disp, maxD, minD, prevStage, minUnlDisp, maxRelDisp);

	if (stage == "backbone-positive" || stage == "backbone-negative") {
		double dAbs = fabs(disp); //even if displacement is negative, stiffness is positive
		if (disp == 0)
		{
			return _iniStiff;
		}

		if (dAbs <= _dMax)
		{
			double force = _fMax * (dAbs / _dMax)*(n / (n - 1 + pow((dAbs / _dMax), n)));
			return abs(force / realDisp);
		}
		else if (dAbs > _dMax & dAbs <= _dUlt)
		{
			double force = (_fMax - (dAbs - _dMax) * _degStiff);
			return abs(force / realDisp);
		}
		else
		{
			double force = _fRes;
			return abs(force / realDisp);
		}
	}
	else if (stage == "unloading-from-positive") {
		double dAbs = fabs(disp); //even if displacement is negative, stiffness is positive
		double maxForce = 0;
		if (maxD < _dMax) { //below maximum but above 80% of it
			maxForce = _fMax * (maxD / _dMax)*(n / (n - 1 + pow((maxD / _dMax), n)));
		}
		else {
			maxForce = (_fMax - (maxD - _dMax) * _degStiff);
		}

		double unlDisp = maxD - (maxForce - (-_fUnl)) / _unlStiff;
		if (disp > unlDisp) {//spring is in the first branch of the unloading
			return abs((maxForce - (maxD - disp)*_unlStiff) / realDisp);
		}
		else { //spring is in the connection brach OR the reloading branch
			double relDisp = (maxForce + (-_fUnl) - _conStiff * unlDisp - _relStiff * maxD) / (_relStiff - _conStiff);
			double relForce = _conStiff * (relDisp - unlDisp) + (-_fUnl);
			if (disp > relDisp) { //spring is in the connection branch
				double force = (-_fUnl) - (unlDisp - disp)*_conStiff;
				return abs(force / realDisp);
			}
			else { //spring is in the reloading branch
				double force = relForce - (relDisp - disp)*_relStiff;
				return abs(force / realDisp);
			}
		}

	}
	else if (stage == "reloading-from-negative") {
		double maxForce = 0;
		if (minD > -_dMax) { //below maximum but above 80% of it
			maxForce = (-_fMax) * (minD / (-_dMax))*(n / (n - 1 + pow((minD / (-_dMax)), n)));
		}
		else {
			maxForce = ((-_fMax) - (minD - (-_dMax)) * _degStiff);
		}
		double unlDisp = minD - (maxForce - _fUnl) / _unlStiff;
		if (disp < unlDisp) {//spring is in the first branch of the unloading
			return abs((maxForce - (minD - disp)*_unlStiff) / realDisp);
		}
		else { //spring is in the connection brach OR the reloading branch
			double relDisp = (maxForce + _fUnl - _conStiff * unlDisp - _relStiff * minD) / (_relStiff - _conStiff);
			double relForce = _conStiff * (relDisp - unlDisp) + _fUnl;
			if (disp < relDisp) { //spring is in the connection branch
				double force = _fUnl - (unlDisp - disp)*_conStiff;
				return abs(force / realDisp);
			}
			else { //spring is in the reloading branch
				double force = relForce - (relDisp - disp)*_relStiff;
				return abs(force / realDisp);
			}
		}
	}
	else if (stage == "unload-reload-connection") {
		double prevForce = GetUnloadForce(-(minUnlDisp - maxD), maxD);
		return abs((_unlStiff * (minUnlDisp - (maxD - disp)) + prevForce) / realDisp);
	}
	else if (stage == "reload-unload-connection") {
		double prevForce = GetReloadForce((maxRelDisp + minD), minD);
		return abs((prevForce - (_unlStiff * (maxRelDisp - (disp - minD)))) / realDisp);
	}
	else {
		return -1; //some problem, unpredicted behavior
	}

double dAbs = fabs(disp); //even if displacement is negative, stiffness is positive
if (dAbs == 0)
{
	return _iniStiff;
}

if (dAbs <= _dMax)
{
	double force = _fMax * (dAbs / _dMax)*(n / (n - 1 + pow((dAbs / _dMax), n)));
	return force / dAbs;
}
else if (dAbs > _dMax & dAbs <= _dUlt)
{
	double force = (_fMax - (dAbs - _dMax) * _degStiff);
	return force / dAbs;
}
else
{
	double force = _fRes;
	return force / dAbs;
}
}

double SpringGeneralModel::GetSecantStiffnessFromForce(double force) {
	double fAbs = fabs(force); //even if the force is negative, stiffness is positive
	double n = GetNValue();

	if (fAbs == 0)
	{
		return _iniStiff;
	}

	if (fAbs <= _fMax)
	{
		double curveForce = 0;
		double disp = _dMax;
		do {
			force = _fMax * (disp / _dMax)*(n / (n - 1 + pow((disp / _dMax), n)));
			double ratio = fAbs / force;
			disp *= ratio;
		} while ((force - fAbs) > 0.001);
		return fAbs / disp;
	}
	else {
		//if the force is greater than the peak force, the analysis by force cannot proceed
		return -1;
	}
}

double SpringGeneralModel::GetYieldDisp() {
	double n = GetNValue();
	double force = 0;
	double elaFract = 0.8;
	double fAbs = elaFract * _fMax;
	double disp = _dMax;
	do {
		force = _fMax * (disp / _dMax)*(n / (n - 1 + pow((disp / _dMax), n)));
		double ratio = fAbs / force;
		disp *= ratio;
	} while ((force - fAbs) > 0.001);
	return disp;
}

double SpringGeneralModel::GetInitialStiffness() {
	return _iniStiff;
}

double SpringGeneralModel::GetPlasticDisplacement(double disp, double maxD, double minD, std::string prevStage, std::string stage, double minUnlDisp, double maxRelDisp)
{
	double endElas = GetEndElasticDisplacement(disp, maxD, minD); //this can be positive or negative
	double elastic = GetYieldDisp();

	//added the next two else ifs in order to try to remove the else ifs after that
	if ((prevStage == "backbone-positive" && stage == "unloading-from-positive") || stage == "unloading-from-positive") {
		double plasticDisp = maxD - GetEndElasticDisplacement(maxD, maxD, minD);
		if (disp >= plasticDisp) {
			return plasticDisp;
		}
		else {
			double force = -GetUnloadForce(disp, maxD);
			return disp + force / _unlStiff;
		}
		//return maxD - GetEndElasticDisplacement(maxD, maxD, minD);
	}
	else if ((prevStage == "backbone-negative" && stage == "reloading-from-negative") || stage == "reloading-from-negative") {
		double plasticDisp = minD - GetEndElasticDisplacement(minD, maxD, minD);
		if (disp <= plasticDisp) {
			return plasticDisp;
		}
		else {
			double force = GetReloadForce(disp, minD);
			return disp - force / _unlStiff;
		}
		//return minD - GetEndElasticDisplacement(minD, maxD, minD);
	}
	//end of added else ifs
	else {
		if (stage == "reloading-from-negative") {
			return minD - GetEndElasticDisplacement(minD, maxD, minD);
		}
		else if (stage == "unload-reload-connection") {
			//return -(minUnlDisp - maxD);
			double force = -GetUnloadForce(-(minUnlDisp - maxD), maxD);
			return -(minUnlDisp - maxD) + force / (_unlStiff);
		}
		else if (stage == "reload-unload-connection") {
			//return maxRelDisp + minD;
			double force = GetReloadForce(maxRelDisp + minD, minD);
			return maxRelDisp + minD - force / (_unlStiff);
		}
		else if (stage == "backbone-negative") {
			return 0; //test
			if (disp < -elastic) {
				return disp - endElas;
			}
			else if (disp > -elastic) {
				return 0;
			}
			else {
				return minD - GetEndElasticDisplacement(minD, maxD, minD);
			}
		}
		else if (stage == "backbone-positive") {
			return 0; //test
			if (disp > elastic) {
				return disp - endElas;
			}
			else if (disp < elastic) {
				return 0;
			}
			else {
				return maxD - GetEndElasticDisplacement(maxD, maxD, minD);
			}
		}
		else {
			if (prevStage == "backbone-positive") {
				return maxD - GetEndElasticDisplacement(maxD, maxD, minD);
			}
			else {
				return minD - GetEndElasticDisplacement(minD, maxD, minD);
			}
		}
	}
}

double SpringGeneralModel::GetEndElasticDisplacement(double disp, double maxD, double minD) {
	double force = GetForceFromDisplacement(disp, maxD, minD);
	return force / _unlStiff;
}

double SpringGeneralModel::GetForceFromDisplacement(double disp, double maxD, double minD)
{
	double dAbs = fabs(disp);
	double sign = 1;
	if (disp < 0) {
		sign = -1;
	}
	double n = GetNValue();

	if (dAbs == 0)
	{
		return 0;
	}

	if (dAbs <= _dMax)
	{
		double force = _fMax * (dAbs / _dMax)*(n / (n - 1 + pow((dAbs / _dMax), n)));
		return force * sign;
	}
	else if (dAbs > _dMax & dAbs <= _dUlt)
	{
		double force = (_fMax - (dAbs - _dMax) * _degStiff);
		return force * sign;
	}
	else
	{
		double force = _fRes;
		return force * sign;
	}
}

std::string SpringGeneralModel::GetType()
{
	return "Spring-General";
}

std::string SpringGeneralModel::ToString()
{
	std::string str = "";
	str += "(";
	str += std::to_string(_ID);
	str += ")";
	str += "(";
	str += "Initial Stiff = ";
	str += std::to_string(_iniStiff);
	str += ", ";
	str += "Peak Displacement = ";
	str += std::to_string(_dMax);
	str += ", ";
	str += "Peak Force = ";
	str += std::to_string(_fMax);
	str += ", ";
	str += "Residual Force = ";
	str += std::to_string(_fRes);
	str += ", ";
	str += "Ultimate Displacement = ";
	str += std::to_string(_dUlt);
	str += ", ";
	str += "Unload Stiffness = ";
	str += std::to_string(_unlStiff);
	str += ", ";
	str += "Unload Force = ";
	str += std::to_string(_fUnl);
	str += ", ";
	str += "Connecting Stiffness = ";
	str += std::to_string(_conStiff);
	str += ", ";
	str += "Reload Stiffness = ";
	str += std::to_string(_relStiff);
	str += ")";
	return str;
}

bool SpringGeneralModel::IsConnectingFromUnload(double disp, double maxD, double prevDisp) {
	double maxForce = 0;
	double n = GetNValue();
	if (maxD < _dMax) { //below maximum but above yield of it
		maxForce = (-_fMax) * ((-maxD) / (-_dMax))*(n / (n - 1 + pow(((-maxD) / (-_dMax)), n)));
	}
	else {
		maxForce = ((-_fMax) - ((-maxD) - (-_dMax)) * _degStiff);
	}
	double unlDisp = (-maxD) - (maxForce - _fUnl) / _unlStiff;
	double relDisp = (maxForce + _fUnl - _conStiff * unlDisp - _relStiff * (-maxD)) / (_relStiff - _conStiff);
	double relForce = _conStiff * (relDisp - unlDisp) + _fUnl;

	double force;
	if (disp > (-maxD) && disp < unlDisp) {
		force = _unlStiff * (disp + maxD) - maxForce;
	}
	else if (disp > unlDisp && disp < relDisp) {
		force = _conStiff * (disp - unlDisp) + _fUnl;
	}
	else if (disp > relDisp && disp < maxD){
		force = _relStiff * (disp - relDisp) + relForce;
	}
	else {
		force = -1;//some problem
	}

	double prevForce = GetUnloadForce(-(prevDisp - maxD), maxD);
	if (_unlStiff*(prevDisp - (maxD - disp)) + prevForce < force) {
		return true; //it is in the connecting branch between unload and reload
	}
	else {
		return false; //it is already in the oposite branch
	}
}

double SpringGeneralModel::GetUnloadForce(double disp, double maxD) {
	double maxForce = 0;
	double n = GetNValue();
	if (maxD < _dMax) { //below maximum but above yield of it
		maxForce = _fMax * (maxD / _dMax)*(n / (n - 1 + pow((maxD / _dMax), n)));
	}
	else {
		maxForce = (_fMax - (maxD - _dMax) * _degStiff);
	}

	double unlDisp = maxD - (maxForce - (-_fUnl)) / _unlStiff;
	if (disp > unlDisp) {//spring is in the first branch of the unloading
		return (maxForce - (maxD - disp)*_unlStiff);
	}
	else { //spring is in the connection brach OR the reloading branch
		double relDisp = (maxForce + (-_fUnl) - _conStiff * unlDisp - _relStiff * maxD) / (_relStiff - _conStiff);
		double relForce = _conStiff * (relDisp - unlDisp) + (-_fUnl);
		if (disp > relDisp) { //spring is in the connection branch
			return (-_fUnl) - (unlDisp - disp)*_conStiff;
		}
		else { //spring is in the reloading branch
			return relForce - (relDisp - disp)*_relStiff;
		}
	}
}

double SpringGeneralModel::GetReloadForce(double disp, double minD) {
	double maxForce = 0;
	double n = GetNValue();
	if (minD > -_dMax) { //below maximum but above yield of it
		maxForce = (-_fMax) * (minD / (-_dMax))*(n / (n - 1 + pow((minD / (-_dMax)), n)));
	}
	else {
		maxForce = ((-_fMax) - (minD - (-_dMax)) * _degStiff);
	}
	double unlDisp = minD - (maxForce - _fUnl) / _unlStiff;
	if (disp < unlDisp) {//spring is in the first branch of the unloading
		return (maxForce - (minD - disp)*_unlStiff);
	}
	else { //spring is in the connection brach OR the reloading branch
		double relDisp = (maxForce + _fUnl - _conStiff * unlDisp - _relStiff * minD) / (_relStiff - _conStiff);
		double relForce = _conStiff * (relDisp - unlDisp) + _fUnl;
		if (disp < relDisp) { //spring is in the connection branch
			return _fUnl - (unlDisp - disp)*_conStiff;
		}
		else { //spring is in the reloading branch
			return relForce - (relDisp - disp)*_relStiff;
		}
	}
}

bool SpringGeneralModel::IsConnectingFromReload(double disp, double minD, double prevDisp) {
	double maxForce;
	double n = GetNValue();
	if ((-minD) < _dMax) { //below maximum but above yield of it
		maxForce = _fMax * ((-minD) / _dMax)*(n / (n - 1 + pow(((-minD) / _dMax), n)));
	}
	else {
		maxForce = (_fMax - ((-minD) - _dMax) * _degStiff);
	}

	double unlDisp = (-minD) - (maxForce - (-_fUnl)) / _unlStiff;
	double relDisp = (maxForce + (-_fUnl) - _conStiff * unlDisp - _relStiff * (-minD)) / (_relStiff - _conStiff);
	double relForce = _conStiff * (relDisp - unlDisp) + (-_fUnl);

	double force;
	if (disp < (-minD) && disp > unlDisp) {
		force = maxForce - _unlStiff * ((-minD) - disp);
	}
	else if (disp < unlDisp && disp > relDisp) {
		force = (-_fUnl) - _conStiff * (unlDisp - disp);
	}
	else if (disp < relDisp && disp > minD) {
		force = relForce - _relStiff * (relDisp - disp);
	}
	else {
		force = -1;//some problem
	}

	double prevForce = GetReloadForce((prevDisp + minD), minD);
	if (prevForce - (_unlStiff*(prevDisp - (disp - minD))) > force) {
		return true; //it is in the connecting branch between unload and reload
	}
	else {
		return false; //it is already in the oposite branch
	}

}

std::string SpringGeneralModel::GetLoadingStage(double disp, double maxD, double minD, std::string prevStage, double minUnlDisp, double maxRelDisp)
{
	double yieldDisp = GetYieldDisp();
	if (prevStage == "initial") { //those are the stages the spring can go after the initial
		if (disp >= 0) {
			return "backbone-positive";
		}
		else if (disp < 0) {
			return "backbone-negative";
		}
		else { return "somewhere"; }
	}
	else if (prevStage == "backbone-positive") {
		double ratio = maxD / yieldDisp;
		//if (disp == maxD || (maxD < yieldDisp && disp >= 0 && (ratio < 0.98 || ratio > 1.02))) { //80% of the max is the 'yield/elastic' criteria
		if (disp == maxD || (maxD < yieldDisp && disp >= 0)) { //80% of the max is the 'yield/elastic' criteria	
			return "backbone-positive";
		}
		if (disp < 0 && maxD < yieldDisp) {
			return "backbone-negative";
		}
		else if (disp < maxD && disp > -maxD) {
			return "unloading-from-positive";
		}
		else { return "somewhere"; }
	}
	else if (prevStage == "backbone-negative") {
		double ratio = -minD / yieldDisp;
		//if (disp == minD || (-minD < yieldDisp && disp <= 0 && (ratio < 0.98 || ratio > 1.02))) { //80% of the max is the 'yield/elastic' criteria
		if (disp == minD || (-minD < yieldDisp && disp <= 0)) { //80% of the max is the 'yield/elastic' criteria	
			return "backbone-negative";
		}
		else if (disp > 0 && -minD < yieldDisp) {
			return "backbone-positive";
		}
		else if (disp > minD && disp < -minD) {
			return "reloading-from-negative";
		}
		else { return "somewhere"; }
	}
	else if (prevStage == "unloading-from-positive") {
		double ratio = abs(disp / maxD);
		if (disp < maxD && disp > -maxD) {
			//double ratio = (maxD - disp) / minUnlDisp;
			double ratio = 2;
			if ((ratio < 0.99 || ratio > 1.01) && (maxD - disp) < minUnlDisp && minUnlDisp != 0) {
				if (IsConnectingFromUnload(disp, maxD, minUnlDisp)) {
					return "unload-reload-connection";
				}
				else {
					return "reloading-from-negative";
				}
			}
			else {
				return "unloading-from-positive";
			}
		}
		else if (disp >= maxD)
		{
			return "backbone-positive";
		}
		else if (disp == minD || disp < -maxD) {
			return "backbone-negative";
		}
		else { return "somewhere"; }
	}
	else if (prevStage == "reloading-from-negative") {
		double ratio = abs(disp / minD);
		if ((disp > minD && disp < -minD)) {
			//double ratio = (disp - minD) / maxRelDisp;
			double ratio = 2;
			if ((ratio < 0.99 || ratio > 1.1) && (disp - minD) < maxRelDisp && maxRelDisp != 0) {
				if (IsConnectingFromReload(disp, minD, maxRelDisp)) {
					return "reload-unload-connection";
				}
				else {
					return "unloading-from-positive";
				}
			}
			else {
				return "reloading-from-negative";
			}
		}
		else if (disp <= minD) {
			return "backbone-negative";
		}
		else if (disp == maxD || disp > -minD) {
			return "backbone-positive";
		}
		else { return "somewhere"; }
	}
	else if (prevStage == "unload-reload-connection") {
		double relForce = GetReloadForce(disp, -maxD);
		double prevForce = GetUnloadForce(-(minUnlDisp - maxD), maxD);

		if ((maxD - disp) < minUnlDisp) {
			if (_unlStiff*(minUnlDisp - (maxD - disp)) + prevForce > relForce) {
				return "reloading-from-negative";
			}
			else {
				return "unload-reload-connection";
			}
		}
		else {
			return "unloading-from-positive";
		}
	}
	else if (prevStage == "reload-unload-connection") {
		double unlForce = GetUnloadForce(disp, -minD);
		double prevForce = GetReloadForce((maxRelDisp + minD), minD);

		if ((disp - minD) < maxRelDisp) {
			if (prevForce - _unlStiff*(maxRelDisp - (disp - minD)) < unlForce) {
				return "unloading-from-positive";
			}
			else {
				return "reload-unload-connection";
			}
		}
		else {
			return "reloading-from-negative";
		}
	}
	else {
		return "somewhere";
	}
}

void SpringGeneralModel::UpdateUnlAndRelDisps(std::string stage, std::string prevStage, double disp, double& maxD, double& minD, double& unlDisp, double& relDisp, double prevUnlDisp, double prevRelDisp) {
	if (stage == "unloading-from-positive" && (prevStage == "reloading-from-negative" || prevStage == "reload-unload-connection")) {
		maxD = -minD; //if I was reloading or in the connection branch and am now unloading, I need to set new values for max values so it calcualtes the branches correctly
		unlDisp = maxD - disp;
		relDisp = 0;
	}
	if (stage == "reloading-from-negative" && (prevStage == "unloading-from-positive" || prevStage == "unload-reload-connection")) {
		minD = -maxD; //if I was reloading or in the connection branch and am now unloading, I need to set new values for max values so it calcualtes the branches correctly
		relDisp = disp - minD;
		unlDisp = 0;
	}
	
	double newUnlDisp = maxD - disp;
	double newRelDisp = disp - minD;
	if (stage == "unloading-from-positive" && newUnlDisp >= prevUnlDisp) {
		unlDisp = newUnlDisp;
	}
	if (stage == "reloading-from-negative" && newRelDisp >= prevRelDisp) {
		relDisp = newRelDisp;
	}
	if (stage == "backbone-negative" || stage == "backbone-positive") {
		relDisp = 0;
		unlDisp = 0;
	}
	if (stage == "reload-unload-connection" && prevStage == "reload-unload-connection") {
		relDisp = prevRelDisp;
	}
	if (stage == "unload-reload-connection" && prevStage == "unload-reload-connection") {
		unlDisp = prevUnlDisp;
	}
}

SpringGeneralModel::~SpringGeneralModel()
{
}


