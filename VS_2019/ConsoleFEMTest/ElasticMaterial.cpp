#include "pch.h"
#include "ElasticMaterial.h"
#include <vector>
#include <string>


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
