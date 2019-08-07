#pragma once
#include "MaterialModel.h"
#include <string>

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
	static OrthotropicElasticMaterial FindElasticMaterialByID(std::vector<OrthotropicElasticMaterial> listOfShellMaterials, int ID);
	~OrthotropicElasticMaterial();

private:
	int _ID;
	double _Ex, _Ey, _vxy, _Gxy, _Gyz, _Gxz;
};

