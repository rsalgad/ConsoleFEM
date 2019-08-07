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

std::string OrthotropicElasticMaterial::ToString()
{
	std::string str = "";
	str += "(";
	str += std::to_string(_ID);
	str += ")";
	str += "(";
	str += "Ex = ";
	str += _Ex;
	str += ", ";
	str += "Ey = ";
	str += _Ey;
	str += ", ";
	str += "Ex = ";
	str += _Ex;
	str += ", ";
	str += "Vxy = ";
	str += _vxy;
	str += ", ";
	str += "Gxy = ";
	str += _Gxy;
	str += ", ";
	str += "Gyz = ";
	str += _Gyz;
	str += ", ";
	str += "Gxz = ";
	str += _Gxz;
	str += ")";
	return str;
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
