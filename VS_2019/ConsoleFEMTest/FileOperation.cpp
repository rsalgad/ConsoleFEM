#include "pch.h"
#include "FileOperation.h"
#include "MonotonicAnalysis.h"
#include "CyclicAnalysis.h"
#include "ReverseCyclic.h"
#include "DynamicAnalysis.h"
#include "tinyxml2.h"
#include <map>
#include <fstream>

using namespace tinyxml2;

FileOperation::FileOperation()
{
}


FileOperation::~FileOperation()
{
}


void FileOperation::SaveResultsFile(std::string& fileName, const StructureManager* structManager, const NodalRecorder<Node>* dispRecorder, const NodalRecorder<Load>* forceRecorder) {
	std::vector<double> emptyVec;
	Matrix emptyMatrix(0);
	SaveResultsFile(fileName, structManager, dispRecorder, forceRecorder, emptyVec, emptyMatrix);

}

void FileOperation::SaveResultsFile(std::string &fileName, const StructureManager* structManager, const NodalRecorder<Node>* dispRecorder, const NodalRecorder<Load>* forceRecorder, std::vector<double> &natFreq, Matrix &modeShape) {

	int steps = dispRecorder->GetRecord()->size();
	int nodes = dispRecorder->GetRecord()->find(1)->second.size();
	int forces = forceRecorder->GetRecord()->find(1)->second.size();

	XMLDocument* doc = new XMLDocument();
	XMLDeclaration* dec = doc->NewDeclaration();
	doc->LinkEndChild(dec);

	XMLElement* rootEle = doc->NewElement("result");
	rootEle->SetAttribute("total", steps);
	doc->LinkEndChild(rootEle);

	if (natFreq.size() != 0) {

		XMLElement* natFreqEle = doc->NewElement("frequencies");

		for (int i = 0; i < modeShape.GetDimY(); i++) {
			XMLElement* freqEle = doc->NewElement("frequency");
			freqEle->SetAttribute("value", natFreq[i]);
			XMLText* freq = doc->NewText(std::to_string(i + 1).c_str());
			freqEle->LinkEndChild(freq);
			natFreqEle->LinkEndChild(freqEle);
		}

		rootEle->LinkEndChild(natFreqEle);

	}

	if (modeShape.GetDimX() != 0) {
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
	}

	for (int i = 0; i < steps; i++) {
		std::string txtStep = "loadstep" + std::to_string(i + 1);
		XMLElement* stepRoot = doc->NewElement(txtStep.c_str());
		XMLElement* dispRoot = doc->NewElement("displacements");
		dispRoot->SetAttribute("total", nodes);
		for (int j = 0; j < nodes; j++) {
			XMLElement* disp = doc->NewElement("displacement");
			std::string txtstr = std::to_string((*dispRecorder->Nodes())[j]);
			XMLText* ID = doc->NewText(txtstr.c_str());
			disp->SetAttribute("x", *dispRecorder->GetRecord()->find(i + 1)->second.find(j + 1)->second.GetX() - *structManager->Nodes()->find(j + 1)->second->GetX());
			disp->SetAttribute("y", *dispRecorder->GetRecord()->find(i + 1)->second.find(j + 1)->second.GetY() - *structManager->Nodes()->find(j + 1)->second->GetY());
			disp->SetAttribute("z", *dispRecorder->GetRecord()->find(i + 1)->second.find(j + 1)->second.GetZ() - *structManager->Nodes()->find(j + 1)->second->GetZ());
			disp->LinkEndChild(ID);
			dispRoot->LinkEndChild(disp);
		}
		stepRoot->LinkEndChild(dispRoot);

		XMLElement* forceRoot = doc->NewElement("forces");
		forceRoot->SetAttribute("total", nodes);

		int j = 0;
		int count = 1;
		for (int j = 0; j < nodes; j++) {
			XMLElement* force = doc->NewElement("force");
			std::string txtstr = std::to_string((*forceRecorder->Nodes())[j]);
			XMLText* ID = doc->NewText(txtstr.c_str());
			force->SetAttribute("x", forceRecorder->GetRecord()->find(i + 1)->second.find(j + 1)->second.GetLoadVector()[0][1]);
			force->SetAttribute("y", forceRecorder->GetRecord()->find(i + 1)->second.find(j + 1)->second.GetLoadVector()[1][1]);
			force->SetAttribute("z", forceRecorder->GetRecord()->find(i + 1)->second.find(j + 1)->second.GetLoadVector()[2][1]);
			force->LinkEndChild(ID);
			forceRoot->LinkEndChild(force);
		}
		stepRoot->LinkEndChild(forceRoot);
		rootEle->LinkEndChild(stepRoot);
	}
	std::string file = fileName + ".res";
	doc->SaveFile(file.c_str());
}


void FileOperation::ReadInputFromXML(std::string fileName, StructureManager& structManager, AnalysisMethod* *analysis) {

	int matTotal, loadTotal, supTotal, eleTotal, nodeTotal, springTotal, massTotal, seismicTotal, impulseTotal, impulseNodeTotal;

	tinyxml2::XMLDocument doc;
	std::string file = fileName + ".xml";
	doc.LoadFile(file.c_str());

	tinyxml2::XMLNode* root = doc.FirstChildElement("structure"); //find the root of the xml file
	tinyxml2::XMLElement* element = root->FirstChildElement("materials"); //finds the next 'child' element, which is the <materials> tag
	element->QueryIntAttribute("total", &matTotal); //finds the attribute of the tag, i.e., the number in between tags eg. <tag> atribute </tag>

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
		}
		element = element->NextSiblingElement("material");
	}

	element = root->FirstChildElement("loads");
	element->QueryIntAttribute("total", &loadTotal); //finds the attribute of the tag, i.e., the number in between tags eg. <tag> atribute </tag>
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
		if (_status != "seismic" || _status != "impulse") {
			*analysis = new MonotonicAnalysis(10, 100, 0.0001);;
		}

		element = root->FirstChildElement("loads")->FirstChildElement("load");
		for (int k = 0; k < i + 1; k++) {
			element = element->NextSiblingElement("load");
		}
	}

	if (root->FirstChildElement("masses") != 0) {
		element = root->FirstChildElement("masses");
		element->QueryIntAttribute("total", &massTotal); //finds the attribute of the tag, i.e., the number in between tags eg. <tag> atribute </tag>
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

			element = root->FirstChildElement("masses")->FirstChildElement("mass");
			for (int k = 0; k < i + 1; k++) {
				element = element->NextSiblingElement("mass");
			}
		}
	}

	element = root->FirstChildElement("boundaries");
	element->QueryIntAttribute("total", &supTotal); //finds the attribute of the tag, i.e., the number in between tags eg. <tag> atribute </tag>
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

		element = root->FirstChildElement("boundaries")->FirstChildElement("boundary");
		for (int k = 0; k < i + 1; k++) {
			element = element->NextSiblingElement("boundary");
		}
	}

	element = root->FirstChildElement("nodes");
	element->QueryIntAttribute("total", &nodeTotal);
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
		element = element->NextSiblingElement("node");
	}

	element = root->FirstChildElement("shell");
	element->QueryIntAttribute("total", &eleTotal);
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
		vecMat.push_back(*OrthotropicElasticMaterial::FindElasticMaterialByID(structManager.Materials(), matID));
		ShellElement* shell = new ShellElement(ID, structManager.Nodes()->find(n1)->second, structManager.Nodes()->find(n2)->second, structManager.Nodes()->find(n3)->second, structManager.Nodes()->find(n4)->second, structManager.Nodes()->find(n5)->second, structManager.Nodes()->find(n6)->second, structManager.Nodes()->find(n7)->second, structManager.Nodes()->find(n8)->second, structManager.Nodes()->find(n9)->second, thick, layers, vecMat);
		structManager.AddShellElement(shell);
		element = element->NextSiblingElement("element");
	}

	element = root->FirstChildElement("spring3D");
	element->QueryIntAttribute("total", &springTotal);
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
		SpringMaterialModels* matX = static_cast<SpringMaterialModels*>(structManager.Materials()->find(matIDX)->second);
		SpringMaterialModels* matY = static_cast<SpringMaterialModels*>(structManager.Materials()->find(matIDY)->second);
		SpringMaterialModels* matZ = static_cast<SpringMaterialModels*>(structManager.Materials()->find(matIDZ)->second);
		vecMat.push_back(matX);
		vecMat.push_back(matY);
		vecMat.push_back(matZ);
		Spring3D* spring = new Spring3D(ID, structManager.Nodes()->find(n1)->second, structManager.Nodes()->find(n2)->second, vecMat, xDir[0], yDir[0]);
		structManager.AddSpringElement(spring);
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

		SeismicLoad* sLoad = new SeismicLoad();

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
			sLoad->SetTime(timeVec);

			if (recordsX.size() != 0) {
				sLoad->SetRecordX(recordsX);
			}
			if (recordsY.size() != 0) {
				sLoad->SetRecordY(recordsY);
			}
			if (recordsZ.size() != 0) {
				sLoad->SetRecordZ(recordsZ);
			}
		}
		*analysis = new DynamicAnalysis(timeVec[timeVec.size() - 1], 0.005, AnalysisTypes::Seismic, IntegrationMethod::AverageNewmark, sLoad, 100, 0.001);
	}

	element = root->FirstChildElement("impulse");
	if (element != 0) {
		element->QueryIntAttribute("total", &impulseTotal);
		element = element->FirstChildElement("data-point"); //finds the next 'child' element, which is the <material> tag

		std::vector<std::vector<double>> points;

		ImpulseLoad* impLoad = new ImpulseLoad();

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

		impLoad->SetPoints(points);
		
		*analysis = new DynamicAnalysis(points[points.size() - 1][0], 0.0005, AnalysisTypes::Impulse, IntegrationMethod::AverageNewmark, impLoad, 100, 0.001);

	}
}

/*
void FileOperation::SaveIterationsResult(std::string fileName, const NodalRecorder<Node>* disps, const std::map<int, Node*>* nodes) {
	
	std::ofstream myfile;
	std::map<int, std::map<int, Node*>>::const_iterator it = disps->GetRecord()->begin();

	while (it != disps->GetRecord()->end()) { //number of records that were recorded
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
*/

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


