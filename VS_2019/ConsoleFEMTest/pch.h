// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

#ifndef PCH_H
#define PCH_H

// TODO: add headers that you want to pre-compile here
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <thread>
#include <mutex>
#include "Matrix.h"
#include "MatrixOperation.h"
#include "VectorOperation.h"
#include "Node.h"
#include "Load.h"
#include "Support.h"
#include "ShellElement.h"
#include "Spring3D.h"
#include "Mass.h"
#include "MaterialModel.h"
#include "ElasticMaterial.h"
#include "OrthotropicElasticMaterial.h"
#include "SpringMaterialModels.h"
#include "SpringAxialModel.h"
#include "SpringGeneralModel.h"
#include "DynamicLoad.h"
#include "SeismicLoad.h"
#include "ImpulseLoad.h"
#include "DynamicAnalysis.h"
#include "StructureManager.h"
#include "NodalRecorder.h"
#include "AnalysisMethod.h"
#include "TimeIntegrationMethod.h"
#include "PreAnalysisSetUp.h"
#include "AnalysisSpringRecorder.h"
#include "Element.h"
#include "Displacement.h"


#endif //PCH_H
