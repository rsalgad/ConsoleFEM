#include "pch.h"
#include "SpringAxialModel.h"
#include <math.h>
#include <iostream>


SpringAxialModel::SpringAxialModel(int ID, double iniStiff, double dMax, double fMax, double degStiff, double fRes, double dUlt, double compStiff, double unlStiff, double fUnl, double conStiff, double relStiff)
{
	_ID = ID;
	_iniStiff = iniStiff;
	_dMax = dMax;
	_fMax = fMax;
	_degStiff = degStiff;
	_fRes = fRes;
	_dUlt = dUlt;
	_compStiff = compStiff;
	_unlStiff = unlStiff;
	_fUnl = fUnl;
	_conStiff = conStiff;
	_relStiff = relStiff;
}

int SpringAxialModel::GetID() {
	return _ID;
}

double SpringAxialModel::GetSecantStiffnessFromDisplacement(double disp, double plasticDisp, double maxD, double minD, std::string prevStage, std::string stage, double minUnlDisp, double maxRelDisp) {
	double n = GetNValue();
	double realDisp = fabs(disp - plasticDisp);

	//std::string status = GetLoadingStage(disp, maxD, minD);
	//std::cout << status << std::endl;

	//std::string stage = GetLoadingStage(disp, maxD, minD, prevStage, minUnlDisp, maxRelDisp);

	if (stage == "backbone-positive") {
		double dAbs = fabs(disp); //even if displacement is negative, stiffness is positive
		if (disp == 0) 
		{
			return _iniStiff;
		}

		if (dAbs <= _dMax)
		{
			double force = _fMax * (dAbs / _dMax)*(n / (n - 1 + pow((dAbs / _dMax), n)));
			return force / realDisp;
		}
		else if (dAbs > _dMax & dAbs <= _dUlt)
		{
			double force = (_fMax - (dAbs - _dMax) * _degStiff);
			return force / realDisp;
		}
		else
		{
			double force = _fRes;
			return force / realDisp;
		}
	}
	else if (stage == "backbone-negative") {
		double force = _compStiff * disp;
		return abs(force / realDisp);
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
			double relDisp = (_conStiff * unlDisp - (-_fUnl)) / (_conStiff - _compStiff);
			double relForce = _compStiff * relDisp;
			if (disp > relDisp) { //spring is in the connection branch
				double force = (-_fUnl) - (unlDisp - disp)*_conStiff;
				return abs(force / realDisp);
			}
			else { //spring is in the reloading branch
				double force = _compStiff * disp;
				return abs(force / realDisp);
			}
		}
	}
	else if (stage == "reloading-from-negative") {
		double relPoint = _fUnl / _unlStiff;
		double maxForce = _fMax * (maxD / _dMax)*(n / (n - 1 + pow((maxD / _dMax), n)));
		double conPoint = (maxForce - _relStiff * maxD - _fUnl + (_conStiff*_fUnl / _unlStiff)) / (_conStiff - _relStiff);
		if (disp < relPoint) {
			if (disp != 0) {
				double force = _unlStiff * disp;
				return abs(force / realDisp);
			}
			else {
				return _unlStiff;
			}
		}
		else if (disp > relPoint && disp < conPoint){
			double force = _fUnl + (disp - relPoint) * _conStiff;
			return abs(force / realDisp);
		}
		else {
			double force = maxForce - (maxD - disp) * _relStiff;
			return abs(force / realDisp);
		}
	}
	else if (stage == "unload-reload-connection") {
		double prevForce = GetUnloadForce(-(minUnlDisp - maxD), maxD);
		return abs((_unlStiff * (minUnlDisp - (maxD - disp)) + prevForce) / realDisp);
	}
	else if (stage == "reload-unload-connection") {
		double prevForce = GetReloadForce((maxRelDisp + minD), maxD);
		return abs((prevForce - (_unlStiff * (maxRelDisp - (disp - minD)))) / realDisp);
	}
	else {
		return -1; //some problem, unpredicted behavior
	}
}

double SpringAxialModel::GetSecantStiffnessFromForce(double force) {
	double fAbs = fabs(force); //even if the force is negative, stiffness is positive
	double n = GetNValue();

	if (force < 0) //if it is being compressed instead of tensioned
	{
		return _compStiff;
	}

	if (fAbs == 0)
	{
		return _iniStiff;
	}

	if (fAbs <= _fMax)
	{
		double force = 0;
		double disp = _dMax;
		do {
			force = _fMax * (disp / _dMax)*(n / (n - 1 + pow((disp / _dMax), n)));
			double ratio = fAbs / force;
			disp *= ratio;
		} while ((force - fAbs) > 0.001);
		return fAbs / disp;
	}
	else 
	{
		//if the force is greater than the peak force, the analysis cannot proceed
		return -1;
	}
}

double SpringAxialModel::GetYieldDisp() {
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

std::string SpringAxialModel::GetType()
{
	return "Spring-Axial";
}

std::string SpringAxialModel::ToString()
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
	str += ", ";
	str += "Compressive Stiffness = ";
	str += std::to_string(_compStiff);
	str += ")";
	return str;
}

double SpringAxialModel::GetEndElasticDisplacement(double disp, double maxD, double minD) {
	double force = GetForceFromDisplacement(disp, maxD, minD);
	return force / _unlStiff;
}

double SpringAxialModel::GetForceFromDisplacement(double disp, double maxD, double minD)
{
	double realDisp = fabs(disp);
	double n = GetNValue();
	double dAbs = fabs(disp); //even if displacement is negative, stiffness is positive
	double force = _fMax * (dAbs / _dMax)*(n / (n - 1 + pow((dAbs / _dMax), n)));
	return force; //returns the stiffness to the point of maximum force back to the origin
	double yieldDisp = GetYieldDisp();

	if (disp == maxD) { //spring is at backbone
		double dAbs = fabs(disp); //even if displacement is negative, stiffness is positive

		if (disp < 0) //if it is being compressed instead of tensioned
		{
			return -_compStiff * disp;
		}

		if (disp == 0) //if it is being compressed instead of tensioned
		{
			return 0;
		}

		if (dAbs <= _dMax)
		{
			double force = _fMax * (dAbs / _dMax)*(n / (n - 1 + pow((dAbs / _dMax), n)));
			return force;
		}
		else if (dAbs > _dMax & dAbs <= _dUlt)
		{
			double force = (_fMax - (dAbs - _dMax) * _degStiff);
			return force;
		}
		else
		{
			double force = _fRes;
			return force;
		}
	}
	else if (disp < 0 && maxD == 0) { //meaning that the spring is under pure compression
		return -_compStiff * disp;
	}
	else if (disp < maxD) //spring is not at backbone (either unloading or reloading)
	{
		double dAbs = fabs(disp); //even if displacement is negative, stiffness is positive
		if (maxD < _dMax && maxD < yieldDisp) { //if the maximum dispalcement was before the maximum load and before 40% the maximum, it is still elastic
			if (disp > 0) {
				double force = _fMax * (dAbs / _dMax)*(n / (n - 1 + pow((dAbs / _dMax), n)));
				return force; //returns the stiffness to the point of maximum force back to the origin
			}
			else {
				return -_compStiff * disp; //compression stiffness once it reaches zero disp coming from the elastic range
			}
		}
		else if (maxD < _dMax) { //below maximum but above 40% of it
			double maxForce = _fMax * (maxD / _dMax)*(n / (n - 1 + pow((maxD / _dMax), n)));
			double unlDisp = maxD - (maxForce - _fUnl) / _unlStiff;
			if (disp > unlDisp) {//spring is in the first branch of the unloading
				return abs((maxForce - (maxD - disp)*_unlStiff));
			}
			else { //spring is in the connection brach OR the reloading branch
				double relDisp = (_conStiff * unlDisp - _fUnl) / (_conStiff - _compStiff);
				double relForce = _compStiff * relDisp;
				if (disp > relDisp) { //spring is in the connection branch
					double force = _fUnl - (unlDisp - disp)*_conStiff;
					return force;
				}
				else { //spring is in the reloading branch
					return -_compStiff * disp;
				}
			}
		}
	}
	else {
		double dAbs = fabs(disp); //even if displacement is negative, stiffness is positive
		if (maxD < _dMax && maxD < yieldDisp) { //if the maximum dispalcement was before the maximum load and before 40% the maximum, it is still elastic
			double force = _fMax * (dAbs / _dMax)*(n / (n - 1 + pow((dAbs / _dMax), n)));
			return force; //returns the stiffness to the point of maximum force back to the origin
		}
		else if (maxD < _dMax) { //below maximum but above 40% of it
			double maxForce = _fMax * (maxD / _dMax)*(n / (n - 1 + pow((maxD / _dMax), n)));
			double unlDisp = maxD - (maxForce - _fUnl) / _unlStiff;
			if (disp > unlDisp) {//spring is in the first branch of the unloading
				return (maxForce - (maxD - disp)*_unlStiff);
			}
			else { //spring is in the connection brach OR the reloading branch
				double relDisp = (_conStiff * unlDisp - _fUnl) / (_conStiff - _compStiff);
				double relForce = _compStiff * relDisp;
				if (disp > relDisp) { //spring is in the connection branch
					double force = _fUnl - (unlDisp - disp)*_conStiff;
					return force;
				}
				else { //spring is in the reloading branch
					return -_compStiff * disp;
				}
			}
		}
	}
}

double SpringAxialModel::DisplacementAtCompressiveCurve(double maxD) {
	double n = GetNValue();
	double maxForce = _fMax * (maxD / _dMax)*(n / (n - 1 + pow((maxD / _dMax), n)));
	double unlDisp = maxD - (maxForce - _fUnl) / _unlStiff;
	return (_conStiff * unlDisp - _fUnl) / (_conStiff - _compStiff);
	
}

std::string SpringAxialModel::GetLoadingStage(double disp, double maxD, double minD, std::string prevStage, double minUnlDisp, double maxRelDisp)
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
		if (disp == maxD || (maxD < yieldDisp && disp >= 0)) { //80% of the max is the 'yield/elastic' criteria
			return "backbone-positive";
		}
		if (disp < 0 && maxD < yieldDisp) {
			return "backbone-negative";
		}
		else if (disp < maxD && disp > DisplacementAtCompressiveCurve(maxD)) {
			return "unloading-from-positive";
		}
		else { return "somewhere"; }
	}
	else if (prevStage == "backbone-negative") {
		if (disp == minD || (disp > minD && disp < 0)) {
			return "backbone-negative";
		}
		else if (disp > minD && (disp == maxD || maxD < yieldDisp)) {
			return "backbone-positive";
		}
		else if (disp >= 0 && maxD > yieldDisp) {
			return "reloading-from-negative";
		} else { return "somewhere"; }
	}
	else if (prevStage == "unloading-from-positive") {
		if (disp < maxD && disp > DisplacementAtCompressiveCurve(maxD)) {
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
		else if (disp < maxD && disp < DisplacementAtCompressiveCurve(maxD)) {
			return "backbone-negative";
		}
		else { return "somewhere"; }
	}
	else if (prevStage == "reloading-from-negative") {
		if (minD < DisplacementAtCompressiveCurve(maxD) && disp < maxD) {
			//double ratio = (disp - minD) / maxRelDisp;
			double ratio = 2;
			if ((ratio < 0.99 || ratio > 1.01) && (disp - minD) < maxRelDisp && maxRelDisp != 0) {
				if (IsConnectingFromReload(disp, maxD, minD, maxRelDisp)) {
					return "reload-unload-connection";
				}
				else {
					return "unloading-from-positive";
				}
			}
			else {
				return "reloading-from-negative";
			}
			return "reloading-from-negative";
		}
		else if (disp == maxD) {
			return "backbone-positive";
		}
		else { return "somewhere"; }
	}
	else if (prevStage == "unload-reload-connection") {
		double relForce = GetReloadForce(disp, maxD); //it is disp here because I want to see what is the force at the reload at this disp. Not the 'prevForce'
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
		double unlForce = GetUnloadForce(disp, maxD);
		double prevForce = GetReloadForce((maxRelDisp + minD), maxD);

		if (disp < DisplacementAtCompressiveCurve(maxD)) {
			return "backbone-negative";
		}
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

bool SpringAxialModel::IsConnectingFromReload(double disp, double maxD, double minD, double prevDisp) {
	//maxD is used because the is no negative minD to account for in the backbone of the axial version
	double maxForce = 0;
	double n = GetNValue();
	if (maxD < _dMax) { //below maximum but above 80% of it
		maxForce = _fMax * (maxD / _dMax)*(n / (n - 1 + pow((maxD / _dMax), n)));
	}
	else {
		maxForce = (_fMax - (maxD - _dMax) * _degStiff);
	}
	double unlDisp = maxD - (maxForce - (-_fUnl)) / _unlStiff;
	double force;
	if (disp > unlDisp) {//spring is in the first branch of the unloading
		force = (maxForce - (maxD - disp)*_unlStiff);
	}
	else { //spring is in the connection brach OR the reloading branch
		double relDisp = (_conStiff * unlDisp - (-_fUnl)) / (_conStiff - _compStiff);
		double relForce = _compStiff * relDisp;
		if (disp > relDisp) { //spring is in the connection branch
			force = (-_fUnl) - (unlDisp - disp)*_conStiff;
		}
		else { //spring is in the reloading branch
			force = _compStiff * disp;
		}
	}

	double prevForce = GetReloadForce(prevDisp + minD, maxD);
	if (prevForce - (_unlStiff*(prevDisp - (disp - minD))) > force) {
		return true; //it is in the connecting branch between unload and reload
	}
	else {
		return false; //it is already in the oposite branch
	}

}

bool SpringAxialModel::IsConnectingFromUnload(double disp, double maxD, double prevDisp) {
	double relPoint = _fUnl / _unlStiff;
	double n = GetNValue();
	double maxForce = _fMax * (maxD / _dMax)*(n / (n - 1 + pow((maxD / _dMax), n)));
	double conPoint = (maxForce - _relStiff * maxD - _fUnl + (_conStiff*_fUnl / _unlStiff)) / (_conStiff - _relStiff);
	double force;
	if (disp < relPoint) {
		force =  _unlStiff * disp;
	}
	else if (disp > relPoint && disp < conPoint) {
		force =  _fUnl + (disp - relPoint) * _conStiff;
	}
	else {
		force =  maxForce - (maxD - disp) * _relStiff;
	}

	double prevForce = GetUnloadForce(disp, maxD);
	if (_unlStiff*(prevDisp - (maxD - disp)) + prevForce < force) {
		return true; //it is in the connecting branch between unload and reload
	}
	else {
		return false; //it is already in the oposite branch
	}
}

double SpringAxialModel::GetUnloadForce(double disp, double maxD) {
	double maxForce = 0;
	double n = GetNValue();
	if (maxD < _dMax) { //below maximum but above 80% of it
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
		double relDisp = (_conStiff * unlDisp - (-_fUnl)) / (_conStiff - _compStiff);
		double relForce = _compStiff * relDisp;
		if (disp > relDisp) { //spring is in the connection branch
			return (-_fUnl) - (unlDisp - disp)*_conStiff;
		}
		else { //spring is in the reloading branch
			return _compStiff * disp;
		}
	}
}

double SpringAxialModel:: GetReloadForce(double disp, double maxD) {
	double relPoint = _fUnl / _unlStiff;
	double n = GetNValue();
	double maxForce = _fMax * (maxD / _dMax)*(n / (n - 1 + pow((maxD / _dMax), n)));
	double conPoint = (maxForce - _relStiff * maxD - _fUnl + (_conStiff*_fUnl / _unlStiff)) / (_conStiff - _relStiff);
	double force;
	if (disp < relPoint) {
		return _unlStiff * disp;
	}
	else if (disp > relPoint && disp < conPoint) {
		return _fUnl + (disp - relPoint) * _conStiff;
	}
	else {
		return maxForce - (maxD - disp) * _relStiff;
	}
}

double SpringAxialModel::GetInitialStiffness() {
	return _iniStiff;
}

double SpringAxialModel::GetPlasticDisplacement(double disp, double maxD, double minD, std::string prevStage, std::string stage, double minUnlDisp, double maxRelDisp)
{
	double endElas = GetEndElasticDisplacement(disp, maxD, minD);
	double elastic = GetYieldDisp();
	//std::string stage = GetLoadingStage(disp, maxD, minD, prevStage, minUnlDisp, maxRelDisp);
	/*
	if ((prevStage == "backbone-negative" && stage == "reloading-from-negative")) {
		return 0;
	}
	*/
	if (stage == "reloading-from-negative") {
		if (disp < (_fUnl / _unlStiff)) {
			return 0;
		}
		else {
			double force = GetReloadForce(disp, maxD);
			return disp - force / _unlStiff;
		}
	}
	else {
		if (stage == "unload-reload-connection") {
			//return -(minUnlDisp - maxD);
			double force = -GetUnloadForce(-(minUnlDisp - maxD), maxD);
			return -(minUnlDisp - maxD) + force / (_unlStiff);
		}
		else if (stage == "reload-unload-connection") {
			//return maxRelDisp + minD;
			double force = GetReloadForce(maxRelDisp + minD, maxD);
			return maxRelDisp + minD - force / (_unlStiff);
		}
		//added the next two else ifs to try to remove all the rest after these
		else if ((prevStage == "backbone-positive" && stage == "unloading-from-positive") || stage == "unloading-from-positive") {
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
		else if (stage == "backbone-positive") {
			return 0;
		}
		else if (stage == "backbone-negative") {
			return 0;
		}
		//end of added else ifs
		else if (disp > elastic && stage == "backbone-positive") {
			return disp - endElas;
		}
		else if (disp < elastic && stage == "backbone-positive") {
			return 0;
		}
		else {
			return maxD - GetEndElasticDisplacement(maxD, maxD, minD);
		}
	}
}

void SpringAxialModel::UpdateUnlAndRelDisps(std::string stage, std::string prevStage, double disp, double& maxD, double& minD, double& unlDisp, double& relDisp, double prevUnlDisp, double prevRelDisp) {
	if (stage == "unloading-from-positive" && (prevStage == "reloading-from-negative" || prevStage == "reload-unload-connection")) {
		unlDisp = maxD - disp;
		relDisp = 0;
	}
	if (stage == "reloading-from-negative" && (prevStage == "unloading-from-positive" || prevStage == "unload-reload-connection")) {
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

SpringAxialModel::~SpringAxialModel()
{
}
