#include "pch.h"


SpringMaterialModels::SpringMaterialModels()
{
}

double SpringMaterialModels::GetSecantStiff() {
	return _fMax / _dMax;
}

double SpringMaterialModels::GetNValue() {
	return _iniStiff / (_iniStiff - GetSecantStiff());
}

SpringMaterialModels* SpringMaterialModels::FindSpringMaterialByID(std::vector<MaterialModel*> listOfMaterials, int ID) {
	for (int i = 0; i < listOfMaterials.size(); i++) {
		if (listOfMaterials[i]->GetID() == ID) {
			return (SpringMaterialModels*) listOfMaterials[i];
		}
	}
}

SpringMaterialModels::~SpringMaterialModels()
{
}
