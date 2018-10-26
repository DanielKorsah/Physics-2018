#pragma once

// GLM
#include <glm/glm.hpp>

// project includes
#include "Application.h"
#include "Mesh.h"
#include "Particle.h"
#include "Body.h"

std::vector < std::vector<glm::vec3>> windMesh;

//vector of vec3 is a triangle
//vector of that is the whole cloth
std::vector<glm::vec3> makeTri(std::vector<std::vector<Particle*>> p_matrix, int size, std::vector < std::vector<glm::vec3>>* windMesh);



