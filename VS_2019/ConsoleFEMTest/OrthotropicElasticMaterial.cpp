#include "pch.h"
#include <iostream>
#include "OrthotropicElasticMaterial.h"


OrthotropicElasticMaterial::OrthotropicElasticMaterial()
{
}

OrthotropicElasticMaterial::OrthotropicElasticMaterial(int ID, double Ex, double Ey, double vxy, double Gxy, double Gyz, double Gxz)
{
	_ID = ID;
	_Ex = Ex;
	_Ey = Ey;
	_vxy = vxy;
	_Gxy = Gxy;
	_Gyz = Gyz;
	_Gxz = Gxz;
}

OrthotropicElasticMaterial::OrthotropicElasticMaterial(int ID, double Ex, double Ey)
{
	_ID = ID;
	_Ex = Ex;
	_Ey = Ey;
}

int OrthotropicElasticMaterial::GetID() {
	return _ID;
}

double OrthotropicElasticMaterial::GetStiffnessX()
{
	return _Ex;
}

double OrthotropicElasticMaterial::GetStiffnessY()
{
	return _Ey;
}

double OrthotropicElasticMaterial::GetShearStiffnessXY()
{
	return _Gxy;
}

double OrthotropicElasticMaterial::GetShearStiffnessYZ()
{
	return _Gyz;
}

double OrthotropicElasticMaterial::GetShearStiffnessXZ()
{
	return _Gxz;
}

double OrthotropicElasticMaterial::GetPoissonXY()
{
	return _vxy;
}

std::string OrthotropicElasticMaterial::GetType()
{
	return "Orthotropic-Elastic";
}

OrthotropicElasticMaterial OrthotropicElasticMaterial::FindElasticMaterialByID(std::vector<OrthotropicElasticMaterial> listOfShellMaterials, int ID) {
	for (int i = 0; i < listOfShellMaterials.size(); i++) {
		if (listOfShellMaterials[i].GetID() == ID) {
			return listOfShellMaterials[i];
		}
	}
}

OrthotropicElasticMaterial::~OrthotropicElasticMaterial()
{
}
