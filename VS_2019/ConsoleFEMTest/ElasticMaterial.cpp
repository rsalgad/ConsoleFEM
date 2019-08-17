#include "pch.h"

ElasticMaterial::ElasticMaterial()
{
}

ElasticMaterial::ElasticMaterial(int ID, double E, double v)
{
	_ID = ID;
	_E = E;
	_v = v;
}

int ElasticMaterial::GetID() {
	return _ID;
}

double ElasticMaterial::GetStiffness()
{
	return _E;
}

double ElasticMaterial::GetPoisson()
{
	return _v;
}

std::string ElasticMaterial::GetType()
{
	return "Elastic";
}

std::string ElasticMaterial::ToString()
{
	std::string str = "";
	str += "(";
	str += std::to_string(_ID);
	str += ")";
	str += "(";
	str += "Ex = ";
	str += _E;
	str += ", ";
	str += "Poisson = ";
	str += _v;
	str += ")";
	return str;
}

ElasticMaterial ElasticMaterial::FindElasticMaterialByID(std::vector<ElasticMaterial> listOfShellMaterials, int ID) {
	for (int i = 0; i < listOfShellMaterials.size(); i++) {
		if (listOfShellMaterials[i].GetID() == ID) {
			return listOfShellMaterials[i];
		}
	}
}

ElasticMaterial::~ElasticMaterial()
{
}
