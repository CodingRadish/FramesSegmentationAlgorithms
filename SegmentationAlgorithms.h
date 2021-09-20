#pragma once
#include <ogx/Plugins/EasyPlugin.h>
#include <ogx/Data/Clouds/CloudHelpers.h>
#include <ogx/Data/Clouds/SphericalSearchKernel.h>
#include <ogx/Data/Clouds/KNNSearchKernel.h>
#include <ogx/Data/Clouds/RandomSearchKernel.h>
#include <stack>
#include <chrono>
#include <iostream>

using namespace ogx;
using namespace ogx::Data;


/////////////////////// Custom functions ///////////////////////////////////////////////////////////////

/** \struct
\brief function

Function enables to check if two 3D vectors are pointing in the same direction

\return boolean: true - same direction, false - different directions
*/

bool is_same_direction(Math::Vector3D vector1, Math::Vector3D vector2);

/** \struct
\brief function for std::sort



*/

bool sortbysecdesc(std::pair<int, int>& a, std::pair<int, int>& b);