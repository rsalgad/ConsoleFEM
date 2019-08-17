#pragma once
//#include "pch.h"
#include <string>

class SpringAxialModel : public SpringMaterialModels
{
private:
	double _compStiff;
	int _ID;
	double DisplacementAtCompressiveCurve(double maxD);
	bool IsConnectingFromReload(double disp, double maxD, double minD, double prevDisp);
	bool IsConnectingFromUnload(double disp, double minD, double prevDisp);
	double GetUnloadForce(double disp, double maxD);
	double GetReloadForce(double disp, double minD);


public:
	SpringAxialModel(int ID, double iniStiff, double dMax, double fMax, double degStiff, double fRes, double dUlt, double compStiff, double unlStiff, double fUnl, double conStiff, double relStiff);
	int GetID();
	double GetSecantStiffnessFromDisplacement(double disp, double plasticDisp, double maxD, double minD, std::string prevStage, std::string stage, double minUnlDisp, double maxRelDisp) override;
	double GetSecantStiffnessFromForce(double force) override;
	double GetInitialStiffness() override;
	double GetPlasticDisplacement(double disp, double maxD, double minD, std::string prevStage, std::string stage, double minUnlDisp, double maxRelDisp) override;
	double GetEndElasticDisplacement(double disp, double maxD, double minD) override;
	double GetForceFromDisplacement(double disp, double maxD, double minD) override;
	std::string GetLoadingStage(double disp, double maxD, double minD, std::string prevStage, double minUnlDisp, double maxRelDisp) override;
	double GetYieldDisp();
	void UpdateUnlAndRelDisps(std::string stage, std::string prevStage, double disp, double& maxD, double& minD, double& unlDisp, double& relDisp, double prevUnlDisp, double prevRelDisp) override;
	std::string GetType() override;
	std::string ToString() override;
	~SpringAxialModel();
};

