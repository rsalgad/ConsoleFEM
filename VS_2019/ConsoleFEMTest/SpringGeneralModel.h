#pragma once
//#include "pch.h"
#include <string>

class SpringGeneralModel : public SpringMaterialModels
{
private:
	int _ID;
	bool IsConnectingFromUnload(double disp, double maxD, double prevDisp);
	bool IsConnectingFromReload(double disp, double minD, double prevDisp);
	double GetUnloadForce(double disp, double maxD);
	double GetReloadForce(double disp, double minD);

public:
	SpringGeneralModel(int ID, double iniStiff, double dMax, double fMax, double degStiff, double fRes, double dUlt, double unlStiff, double fUnl, double conStiff, double relStiff);
	int GetID();
	double GetSecantStiffnessFromDisplacement(double disp, double plasticDisp, double maxD, double minD, std::string prevStage, std::string stage, double minUnlDisp, double maxRelDisp) override;
	double GetSecantStiffnessFromForce(double force) override;
	double GetInitialStiffness() override;
	double GetPlasticDisplacement(double disp, double maxD, double minD, std::string prevStage, std::string stage, double minUnlDisp, double maxRelDisp) override;
	double GetEndElasticDisplacement(double disp, double maxD, double minD) override;
	double GetForceFromDisplacement(double disp, double maxD, double minD);
	void UpdateUnlAndRelDisps(std::string stage, std::string prevStage, double disp, double*& maxD, double*& minD, double*& unlDisp, double*& relDisp, double prevUnlDisp, double prevRelDisp) override;
	std::string GetType() override;
	std::string ToString() override;

	std::string GetLoadingStage(double disp, double maxD, double minD, std::string prevStage, double minUnlDisp, double maxRelDisp) override;
	double GetYieldDisp();
	~SpringGeneralModel();
};

