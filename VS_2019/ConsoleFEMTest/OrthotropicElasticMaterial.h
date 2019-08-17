#pragma once
//#include "pch.h"
#include <string>
#include <vector>

class OrthotropicElasticMaterial : public MaterialModel
{
public:
	OrthotropicElasticMaterial();
	OrthotropicElasticMaterial(int ID, double Ex, double Ey);
	OrthotropicElasticMaterial(int ID, double Ex, double Ey, double vxy, double Gxy, double Gyz, double Gxz);
	int GetID();
	double GetStiffnessX();
	double GetStiffnessY();
	double GetShearStiffnessXY();
	double GetShearStiffnessYZ();
	double GetShearStiffnessXZ();
	double GetPoissonXY();
	std::string GetType() override;
	std::string ToString() override;
	static OrthotropicElasticMaterial* FindElasticMaterialByID(const std::map<int, MaterialModel*>* listOfMaterials, int ID);
	~OrthotropicElasticMaterial();

private:
	int _ID;
	double _Ex, _Ey, _vxy, _Gxy, _Gyz, _Gxz;
};

