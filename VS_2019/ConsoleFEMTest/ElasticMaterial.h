#pragma once
//#include "pch.h"
#include "MaterialModel.h"
#include <string>
#include <vector>

class ElasticMaterial : public MaterialModel
{
public:
	ElasticMaterial();
	ElasticMaterial(int ID, double E, double u);
	int GetID();
	double GetStiffness();
	double GetPoisson();
	std::string GetType() override;
	std::string ToString() override;
	static ElasticMaterial FindElasticMaterialByID(std::vector<ElasticMaterial> listOfShellMaterials, int ID);
	~ElasticMaterial();

private:
	int _ID;
	double _E, _v;
};

