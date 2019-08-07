#pragma once
#include "MaterialModel.h"
#include <string>

class SpringMaterialModels : public MaterialModel
{
public:
	double _iniStiff, _dMax, _fMax, _degStiff, _fRes, _dUlt, _unlStiff, _fUnl, _conStiff, _relStiff;
	SpringMaterialModels();
	virtual double GetSecantStiffnessFromDisplacement(double disp, double plasticDisp, double maxD, double minD, std::string prevStage, std::string stage, double minUnlDisp, double maxRelDisp) = 0;
	virtual double GetSecantStiffnessFromForce(double force) = 0;
	virtual double GetEndElasticDisplacement(double disp, double maxD, double minD) = 0;
	virtual int GetID() = 0;
	virtual double GetSecantStiff();
	virtual double GetNValue();
	virtual double GetForceFromDisplacement(double disp, double maxD, double minD) = 0;
	virtual double GetYieldDisp() = 0;
	virtual std::string GetType() = 0;
	virtual std::string ToString() = 0;
	virtual std::string GetLoadingStage(double disp, double maxD, double minD, std::string prevStage, double minUnlDisp, double maxRelDisp) = 0;
	virtual double GetPlasticDisplacement(double disp, double maxD, double minD, std::string prevStage, std::string stage, double minUnlDisp, double maxRelDisp) = 0;
	static SpringMaterialModels* FindSpringMaterialByID(std::vector<MaterialModel*> listOfMaterials, int ID);
	virtual double GetInitialStiffness() = 0;
	virtual void UpdateUnlAndRelDisps(std::string stage, std::string prevStage, double disp, double& maxD, double& minD, double& unlDisp, double& relDisp, double prevUnlDisp, double prevRelDisp) = 0;

	~SpringMaterialModels();
};

