#include "pch.h"
#include "FileOperation.h"
#include "tinyxml2.h"
#include "MaterialModel.h"
#include "ElasticMaterial.h"
#include "SpringAxialModel.h"
#include "SpringGeneralModel.h"
#include "Load.h"
#include "Mass.h"
#include "Support.h"
#include "SeismicLoad.h"
#include "ImpulseLoad.h"
#include "Node.h"
#include "ShellElement.h"
#include "Spring3D.h"
#include "OrthotropicElasticMaterial.h"
#include "StructureManager.h"
#include "NodalRecorder.h"
#include <fstream>

using namespace tinyxml2;

FileOperation::FileOperation()
{
}


FileOperation::~FileOperation()
{
}

void FileOperation::SaveResultsFile(std::string &fileName, std::vector<std::vector<Node>> &nodePerStep, std::vector<Node> &listOfNodes, std::vector<Matrix> &forcePerStep, std::vector<double> &natFreq, Matrix &modeShape) {

	int steps = nodePerStep.size();
	int nodes = nodePerStep[0].size();

	XMLDocument* doc = new XMLDocument();
	XMLDeclaration* dec = doc->NewDeclaration();
	doc->LinkEndChild(dec);

	XMLElement* rootEle = doc->NewElement("result");
	rootEle->SetAttribute("total", steps);
	doc->LinkEndChild(rootEle);

	XMLElement* natFreqEle = doc->NewElement("frequencies");

	for (int i = 0; i < modeShape.GetDimY(); i++) {
		XMLElement* freqEle = doc->NewElement("frequency");
		freqEle->SetAttribute("value", natFreq[i]);
		XMLText* freq = doc->NewText(std::to_string(i + 1).c_str());
		freqEle->LinkEndChild(freq);
		natFreqEle->LinkEndChild(freqEle);
	}

	rootEle->LinkEndChild(natFreqEle);

	XMLElement* modalEle = doc->NewElement("modal");
	modalEle->SetAttribute("total", modeShape.GetDimY());

	for (int i = 0; i < modeShape.GetDimY(); i++) {
		XMLElement* modeEle = doc->NewElement("mode");
		modeEle->SetAttribute("total", modeShape.GetDimX() / 6);
		int count = 0;
		for (int j = 0; j < modeShape.GetDimX(); j += 6) {
			XMLElement* node = doc->NewElement("node");
			node->SetAttribute("x", modeShape.GetMatrixDouble()[j][i]);
			node->SetAttribute("y", modeShape.GetMatrixDouble()[j + 1][i]);
			node->SetAttribute("z", modeShape.GetMatrixDouble()[j + 2][i]);
			XMLText* ID = doc->NewText(std::to_string(count + 1).c_str());
			node->LinkEndChild(ID);
			modeEle->LinkEndChild(node);
			count++;
		}
		modalEle->LinkEndChild(modeEle);
	}

	rootEle->LinkEndChild(modalEle);

	for (int i = 0; i < steps; i++) {
		std::string txtStep = "loadstep" + std::to_string(i + 1);
		XMLElement* stepRoot = doc->NewElement(txtStep.c_str());
		XMLElement* dispRoot = doc->NewElement("displacements");
		dispRoot->SetAttribute("total", nodes);
		for (int j = 0; j < nodePerStep[i].size(); j++) {
			XMLElement* disp = doc->NewElement("displacement");
			std::string txtstr = std::to_string(j + 1);
			XMLText* ID = doc->NewText(txtstr.c_str());
			disp->SetAttribute("x", nodePerStep[i][j].GetX() - listOfNodes[j].GetX());
			disp->SetAttribute("y", nodePerStep[i][j].GetY() - listOfNodes[j].GetY());
			disp->SetAttribute("z", nodePerStep[i][j].GetZ() - listOfNodes[j].GetZ());
			disp->LinkEndChild(ID);
			dispRoot->LinkEndChild(disp);
		}
		stepRoot->LinkEndChild(dispRoot);

		XMLElement* forceRoot = doc->NewElement("forces");
		forceRoot->SetAttribute("total", nodes);

		int j = 0;
		int count = 1;
		while (j < forcePerStep[i].GetDimX()) {
			XMLElement* force = doc->NewElement("force");
			std::string txtstr = std::to_string(count);
			XMLText* ID = doc->NewText(txtstr.c_str());
			force->SetAttribute("x", forcePerStep[i].GetMatrixDouble()[j][0]);
			force->SetAttribute("y", forcePerStep[i].GetMatrixDouble()[j + 1][0]);
			force->SetAttribute("z", forcePerStep[i].GetMatrixDouble()[j + 2][0]);
			force->LinkEndChild(ID);
			forceRoot->LinkEndChild(force);
			j += 6;
			count++;
		}
		stepRoot->LinkEndChild(forceRoot);
		rootEle->LinkEndChild(stepRoot);
	}
	std::string file = fileName + ".res";
	doc->SaveFile(file.c_str());
}

void FileOperation::ReadInputFromXML(std::string fileName, std::vector<MaterialModel*> &listOfMaterials,
	std::vector<OrthotropicElasticMaterial> &listOfShellMaterials,
	std::vector<SpringAxialModel> &listOfSpringAxialMat,
	std::vector<SpringGeneralModel> &listOfSpringGeneralMat,
	std::vector<Load> &listOfLoads,
	std::vector<Support> &listOfSupports,
	std::vector<Node> &listOfNodes,
	std::vector<ShellElement> &listOfShellElements,
	std::vector<Spring3D> &listOfSpringElements,
	std::vector<Mass> &listOfMasses,
	SeismicLoad &sLoad, ImpulseLoad &impLoad,
	StructureManager &structManager) {

	int matTotal, loadTotal, supTotal, eleTotal, nodeTotal, springTotal, massTotal, seismicTotal, impulseTotal, impulseNodeTotal;

	tinyxml2::XMLDocument doc;
	std::string file = fileName + ".xml";
	doc.LoadFile(file.c_str());

	tinyxml2::XMLNode* root = doc.FirstChildElement("structure"); //find the root of the xml file
	tinyxml2::XMLElement* element = root->FirstChildElement("materials"); //finds the next 'child' element, which is the <materials> tag
	element->QueryIntAttribute("total", &matTotal); //finds the attribute of the tag, i.e., the number in between tags eg. <tag> atribute </tag>
	listOfMaterials.reserve(matTotal);
	element = element->FirstChildElement("material"); //finds the next 'child' element, which is the <material> tag
	for (int i = 0; i < matTotal; i++) { //loops through the amount of materials we have
		const char* type;
		char ela[] = "Elastic";
		char springAx[] = "Spring-Axial";
		int ID;
		element->QueryStringAttribute("type", &type); //gets the attribute (the value) associated to the 'elasticity' term.
		if (strcmp("Elastic", type) == 0) {
			double E, v;
			element->QueryDoubleAttribute("elasticity", &E);
			element->QueryDoubleAttribute("poisson", &v);
			element->QueryIntText(&ID); //gets the attribute in between tags.
			ElasticMaterial* mat = new ElasticMaterial(ID, E, v);
			structManager.AddMaterial(mat);
			listOfShellMaterials.emplace_back(ID, E, v); //add to the list of materials
		}
		else if (strcmp("OrthoElastic", type) == 0) {
			double Ex, Ey, vxy, Gxy, Gyz, Gxz;
			element->QueryDoubleAttribute("ex", &Ex);
			element->QueryDoubleAttribute("ey", &Ey);
			element->QueryDoubleAttribute("vxy", &vxy);
			element->QueryDoubleAttribute("gxy", &Gxy);
			element->QueryDoubleAttribute("gyz", &Gyz);
			element->QueryDoubleAttribute("gxz", &Gxz);
			element->QueryIntText(&ID); //gets the attribute in between tags.
			OrthotropicElasticMaterial* mat = new OrthotropicElasticMaterial(ID, Ex, Ey, vxy, Gxy, Gyz, Gxz);
			structManager.AddMaterial(mat);
			listOfShellMaterials.emplace_back(ID, Ex, Ey, vxy, Gxy, Gyz, Gxz); //add to the list of materials
		}
		else if (strcmp("Spring-Axial", type) == 0) {
			double iniStiff, fMax, dMax, degStiff, fRes, dUlt, compStiff, unlStiff, fUnl, conStiff, relStiff;
			element->QueryDoubleAttribute("initialStiffness", &iniStiff);
			element->QueryDoubleAttribute("peakForce", &fMax);
			element->QueryDoubleAttribute("peakDisplacement", &dMax);
			element->QueryDoubleAttribute("degradingStiffness", &degStiff);
			element->QueryDoubleAttribute("residualForce", &fRes);
			element->QueryDoubleAttribute("ultimateDisplacement", &dUlt);
			element->QueryDoubleAttribute("compressiveStiffness", &compStiff);
			element->QueryDoubleAttribute("unloadStiffness", &unlStiff);
			element->QueryDoubleAttribute("unloadForce", &fUnl);
			element->QueryDoubleAttribute("connectStiffness", &conStiff);
			element->QueryDoubleAttribute("reloadStiffness", &relStiff);
			element->QueryIntText(&ID); //gets the attribute in between tags.
			SpringAxialModel* mat = new SpringAxialModel(ID, iniStiff, dMax, fMax, degStiff, fRes, dUlt, compStiff, unlStiff, fUnl, conStiff, relStiff);
			structManager.AddMaterial(mat);
			listOfSpringAxialMat.emplace_back(ID, iniStiff, dMax, fMax, degStiff, fRes, dUlt, compStiff, unlStiff, fUnl, conStiff, relStiff);
			//listOfMaterials.push_back(&listOfSpringAxialMat[listOfSpringAxialMat.size() - 1]); //add to the list of materials
		}
		else {
			double iniStiff, fMax, dMax, degStiff, fRes, dUlt, unlStiff, fUnl, conStiff, relStiff;
			element->QueryDoubleAttribute("initialStiffness", &iniStiff);
			element->QueryDoubleAttribute("peakForce", &fMax);
			element->QueryDoubleAttribute("peakDisplacement", &dMax);
			element->QueryDoubleAttribute("degradingStiffness", &degStiff);
			element->QueryDoubleAttribute("residualForce", &fRes);
			element->QueryDoubleAttribute("ultimateDisplacement", &dUlt);
			element->QueryDoubleAttribute("unloadStiffness", &unlStiff);
			element->QueryDoubleAttribute("unloadForce", &fUnl);
			element->QueryDoubleAttribute("connectStiffness", &conStiff);
			element->QueryDoubleAttribute("reloadStiffness", &relStiff);
			element->QueryIntText(&ID); //gets the attribute in between tags.
			SpringGeneralModel* mat = new SpringGeneralModel(ID, iniStiff, dMax, fMax, degStiff, fRes, dUlt, unlStiff, fUnl, conStiff, relStiff);
			structManager.AddMaterial(mat);
			listOfSpringGeneralMat.emplace_back(ID, iniStiff, dMax, fMax, degStiff, fRes, dUlt, unlStiff, fUnl, conStiff, relStiff);
			//listOfMaterials.push_back(&listOfSpringGeneralMat[listOfSpringGeneralMat.size() - 1]); //add to the list of materials
		}
		element = element->NextSiblingElement("material");
	}

	for (int i = 0; i < listOfSpringAxialMat.size(); i++) {
		listOfMaterials.push_back(&listOfSpringAxialMat[i]); //add to the list of materials
	}
	for (int i = 0; i < listOfSpringGeneralMat.size(); i++) {
		listOfMaterials.push_back(&listOfSpringGeneralMat[i]); //add to the list of materials
	}

	element = root->FirstChildElement("loads");
	element->QueryIntAttribute("total", &loadTotal); //finds the attribute of the tag, i.e., the number in between tags eg. <tag> atribute </tag>
	listOfLoads.reserve(loadTotal);
	element = element->FirstChildElement("load"); //finds the next 'child' element, which is the <material> tag
	for (int i = 0; i < loadTotal; i++) { //loops through the amount of loads we have
		int ID, nodeID, vals;
		const char* status;
		std::vector<std::vector<double>> loadVector;

		element->QueryIntAttribute("ID", &ID);
		element->QueryIntAttribute("nodeID", &nodeID);
		element->QueryStringAttribute("status", &status);
		element->QueryIntAttribute("values", &vals);
		loadVector.reserve(vals);
		element = element->FirstChildElement("value");
		for (int j = 0; j < vals; j++) {
			double val, dir;
			element->QueryDoubleAttribute("direction", &dir);
			element->QueryDoubleText(&val);
			std::vector<double> vec = { dir, val };
			loadVector.emplace_back(vec);
			element = element->NextSiblingElement("value");
		}
		std::string _status(status);
		Load* load = new Load(ID, nodeID, _status);
		load->SetLoadVector(loadVector);
		structManager.AddLoad(load);
		Load l(ID, nodeID, _status);
		l.SetLoadVector(loadVector);
		listOfLoads.emplace_back(l); //add to the list of loads

		element = root->FirstChildElement("loads")->FirstChildElement("load");
		for (int k = 0; k < i + 1; k++) {
			element = element->NextSiblingElement("load");
		}
	}

	if (root->FirstChildElement("masses") != 0) {
		element = root->FirstChildElement("masses");
		element->QueryIntAttribute("total", &massTotal); //finds the attribute of the tag, i.e., the number in between tags eg. <tag> atribute </tag>
		listOfMasses.reserve(massTotal);
		element = element->FirstChildElement("mass"); //finds the next 'child' element, which is the <material> tag
		for (int i = 0; i < massTotal; i++) { //loops through the amount of masses we have
			int ID, nodeID, vals;
			std::vector<std::vector<double>> massVector;

			element->QueryIntAttribute("ID", &ID);
			element->QueryIntAttribute("nodeID", &nodeID);
			element->QueryIntAttribute("values", &vals);
			massVector.reserve(vals);
			element = element->FirstChildElement("value");
			for (int j = 0; j < vals; j++) {
				double val, dir;
				element->QueryDoubleAttribute("direction", &dir);
				element->QueryDoubleText(&val);
				std::vector<double> vec = { dir, val };
				massVector.emplace_back(vec);
				element = element->NextSiblingElement("value");
			}
			Mass* m = new Mass(ID, nodeID);
			m->SetMassVector(massVector);
			structManager.AddMass(m);
			Mass mass(ID, nodeID);
			mass.SetMassVector(massVector);
			listOfMasses.emplace_back(mass); //add to the list of masses

			element = root->FirstChildElement("masses")->FirstChildElement("mass");
			for (int k = 0; k < i + 1; k++) {
				element = element->NextSiblingElement("mass");
			}
		}
	}

	element = root->FirstChildElement("boundaries");
	element->QueryIntAttribute("total", &supTotal); //finds the attribute of the tag, i.e., the number in between tags eg. <tag> atribute </tag>
	listOfSupports.reserve(supTotal);
	element = element->FirstChildElement("boundary"); //finds the next 'child' element, which is the <material> tag
	for (int i = 0; i < supTotal; i++) { //loops through the amount of loads we have
		double val, dir;
		int ID, nodeID, vals;
		std::vector<std::vector<double>> supVector;

		element->QueryIntAttribute("ID", &ID);
		element->QueryIntAttribute("nodeID", &nodeID);
		element->QueryIntAttribute("values", &vals);
		supVector.reserve(vals);
		element = element->FirstChildElement("value");
		for (int j = 0; j < vals; j++) {
			element->QueryDoubleAttribute("direction", &dir);
			element->QueryDoubleText(&val);
			std::vector<double> vec = { dir, val };
			supVector.emplace_back(vec);
			element = element->NextSiblingElement("value");
		}
		Support* s = new Support(ID, nodeID);
		s->SetSupportVector(supVector);
		structManager.AddSupport(s);
		Support sup(ID, nodeID);
		sup.SetSupportVector(supVector);
		listOfSupports.emplace_back(sup); //add to the list of loads

		element = root->FirstChildElement("boundaries")->FirstChildElement("boundary");
		for (int k = 0; k < i + 1; k++) {
			element = element->NextSiblingElement("boundary");
		}
	}

	element = root->FirstChildElement("nodes");
	element->QueryIntAttribute("total", &nodeTotal);
	listOfNodes.reserve(nodeTotal);
	element = element->FirstChildElement("node");
	for (int i = 0; i < nodeTotal; i++) { //loop through every node
		double x, y, z;
		int ID;
		element->QueryDoubleAttribute("X", &x);
		element->QueryDoubleAttribute("Y", &y);
		element->QueryDoubleAttribute("Z", &z);
		element->QueryIntText(&ID);
		Node* n = new Node(ID, x, y, z);
		structManager.AddNode(n);
		listOfNodes.emplace_back(ID, x, y, z); //add to the list of nodes
		element = element->NextSiblingElement("node");
	}

	element = root->FirstChildElement("shell");
	element->QueryIntAttribute("total", &eleTotal);
	listOfShellElements.reserve(eleTotal);
	element = element->FirstChildElement("element");
	for (int i = 0; i < eleTotal; i++) { //loops through every element
		int n1, n2, n3, n4, n5, n6, n7, n8, n9, matID, ID;
		double thick;
		int layers;
		element->QueryDoubleAttribute("thickness", &thick);
		element->QueryIntAttribute("layers", &layers);
		element->QueryIntAttribute("material", &matID);
		element->QueryIntAttribute("N1", &n1);
		element->QueryIntAttribute("N2", &n2);
		element->QueryIntAttribute("N3", &n3);
		element->QueryIntAttribute("N4", &n4);
		element->QueryIntAttribute("N5", &n5);
		element->QueryIntAttribute("N6", &n6);
		element->QueryIntAttribute("N7", &n7);
		element->QueryIntAttribute("N8", &n8);
		element->QueryIntAttribute("N9", &n9);
		element->QueryIntText(&ID);
		std::vector<OrthotropicElasticMaterial> vecMat;
		vecMat.push_back(OrthotropicElasticMaterial::FindElasticMaterialByID(listOfShellMaterials, matID));
		ShellElement* shell = new ShellElement(ID, structManager.Nodes()[n1 - 1], structManager.Nodes()[n2 - 1], structManager.Nodes()[n3 - 1], structManager.Nodes()[n4 - 1], structManager.Nodes()[n5 - 1], structManager.Nodes()[n6 - 1], structManager.Nodes()[n7 - 1], structManager.Nodes()[n8 - 1], structManager.Nodes()[n9 - 1], thick, layers, vecMat);
		structManager.AddShellElement(shell);
		listOfShellElements.emplace_back(ID, &listOfNodes[n1 - 1], &listOfNodes[n2 - 1], &listOfNodes[n3 - 1], &listOfNodes[n4 - 1], &listOfNodes[n5 - 1], &listOfNodes[n6 - 1], &listOfNodes[n7 - 1], &listOfNodes[n8 - 1], &listOfNodes[n9 - 1], thick, layers, vecMat);
		element = element->NextSiblingElement("element");
	}

	element = root->FirstChildElement("spring3D");
	element->QueryIntAttribute("total", &springTotal);
	listOfSpringElements.reserve(springTotal);
	element = element->FirstChildElement("element");
	for (int i = 0; i < springTotal; i++) { //loops through every element
		int n1, n2, matIDX, matIDY, matIDZ, ID;
		const char* xDir;
		const char* yDir;
		element->QueryIntAttribute("N1", &n1);
		element->QueryIntAttribute("N2", &n2);
		element->QueryIntAttribute("materialX", &matIDX);
		element->QueryIntAttribute("materialY", &matIDY);
		element->QueryIntAttribute("materialZ", &matIDZ);
		element->QueryStringAttribute("axial-vec", &xDir);
		element->QueryStringAttribute("shear-vec", &yDir);
		element->QueryIntText(&ID);

		std::vector<SpringMaterialModels*> vecMat;
		vecMat.push_back(SpringMaterialModels::FindSpringMaterialByID(listOfMaterials, matIDX));
		vecMat.push_back(SpringMaterialModels::FindSpringMaterialByID(listOfMaterials, matIDY));
		vecMat.push_back(SpringMaterialModels::FindSpringMaterialByID(listOfMaterials, matIDZ));
		Spring3D* spring = new Spring3D(ID, structManager.Nodes()[n1 - 1], structManager.Nodes()[n2 - 1], vecMat, xDir[0], yDir[0]);
		structManager.AddSpringElement(spring);
		listOfSpringElements.emplace_back(ID, &listOfNodes[n1 - 1], &listOfNodes[n2 - 1], vecMat, xDir[0], yDir[0]);
		element = element->NextSiblingElement("element");
	}


	element = root->FirstChildElement("seismic");
	if (element != 0) {
		element->QueryIntAttribute("total", &seismicTotal); //finds the attribute of the tag, i.e., the number in between tags eg. <tag> atribute </tag>
		element = element->FirstChildElement("data-point"); //finds the next 'child' element, which is the <material> tag

		std::vector<double> recordsX;
		std::vector<double> recordsY;
		std::vector<double> recordsZ;
		std::vector<double> timeVec;

		for (int i = 0; i < seismicTotal; i++) { //loops through the amount of loads we have
			double time, x, y, z;
			element->QueryDoubleAttribute("time", &time);
			timeVec.push_back(time);

			if (element->Attribute("x") != 0) {
				element->QueryDoubleAttribute("x", &x);
				recordsX.push_back(x);
			}
			if (element->Attribute("y") != 0) {
				element->QueryDoubleAttribute("y", &y);
				recordsY.push_back(y);
			}
			if (element->Attribute("z") != 0) {
				element->QueryDoubleAttribute("z", &z);
				recordsZ.push_back(z);
			}
			element = element->NextSiblingElement("data-point");
		}
		if (timeVec.size() != 0) {
			sLoad.SetTime(timeVec);

			if (recordsX.size() != 0) {
				sLoad.SetRecordX(recordsX);
			}
			if (recordsY.size() != 0) {
				sLoad.SetRecordY(recordsY);
			}
			if (recordsZ.size() != 0) {
				sLoad.SetRecordZ(recordsZ);
			}
		}
	}

	element = root->FirstChildElement("impulse");
	if (element != 0) {
		element->QueryIntAttribute("total", &impulseTotal);
		element = element->FirstChildElement("data-point"); //finds the next 'child' element, which is the <material> tag

		std::vector<std::vector<double>> points;

		for (int i = 0; i < impulseTotal; i++) { //loops through the amount of loads we have
			double time, force;
			std::vector<double> vec;
			element->QueryDoubleAttribute("time", &time);
			element->QueryDoubleAttribute("force", &force);
			vec.push_back(time);
			vec.push_back(force);
			points.push_back(vec);

			element = element->NextSiblingElement("data-point");
		}

		impLoad.SetPoints(points);
	}
}

void FileOperation::SaveIterationsResult(std::string fileName, const NodalRecorder<Node*>* disps, const std::map<int, Node*>* nodes) {
	
	std::ofstream myfile;
	std::map<int, std::map<int, Node*>>::const_iterator it = disps->GetRecord().begin();

	while (it != disps->GetRecord().end()) { //number of records that were recorded
		//error handling IO file
		myfile.open(fileName + std::to_string(it->first) + ".txt");

		for (int i = 0; i < nodes->size(); i++) {
			std::map<int, Node*>::const_iterator dispIt = it->second.find(i + 1);
			std::map<int, Node*>::const_iterator nodeIt = nodes->find(i + 1);
			if (nodeIt != nodes->end()) {
				myfile << nodeIt->second->GetX() - nodeIt->second->GetX() << " " << nodeIt->second->GetY() - nodeIt->second->GetY() << " " << nodeIt->second->GetZ() - nodeIt->second->GetZ() << std::endl;
			}
		}
	}
	myfile.close();
}

//I'm pretty sure this function is not used anymore
/*
void FileOperation::SaveIterationsForceResult(std::string fileName, const const NodalRecorder<Load*>* forces) {
	std::ofstream myfile;

	std::map<int, std::map<int, Load*>>::const_iterator it = forces->GetRecord().begin();

	while (it != forces->GetRecord().end()) {

		myfile.open(fileName + std::to_string(it->first + 1) + ".txt");
		double Fx = 0.0, Fy = 0.0, Fz = 0.0;
		for (int j = 0; j < it->second->GetLoadVector().size; j++) {
			if (it->second->GetLoadVector()[j][0] == 1)
			{
				Fx = it->second->GetLoadVector()[j][1];
			}
			else if (it->second->GetLoadVector()[j][0] == 2) {
				Fy = it->second->GetLoadVector()[j][1];
			}
			else if (it->second->GetLoadVector()[j][0] == 3) {
				Fz = it->second->GetLoadVector()[j][1];
			}
		}
		myfile << Fx << " " << Fy << " " << Fz << std::endl;
		myfile.close();
		it++;
	}
}
*/


