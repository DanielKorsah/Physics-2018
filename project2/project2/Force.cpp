#include <iostream>
#include <cmath>
#include "Force.h"
#include "Body.h"

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include "glm/ext.hpp"

glm::vec3 Gravity::apply(float mass, const glm::vec3 & pos, const glm::vec3 & vel)
{

	return m_gravity * mass;
}

glm::vec3 Drag::apply(float mass, const glm::vec3 & pos, const glm::vec3 & vel)
{
	float velSquared = glm::length(vel) * glm::length(vel);
	glm::vec3 velNorm = glm::normalize(vel);
	glm::vec3 dragForce = 0.5f * density * velSquared * coEff * 1 * velNorm;
	return dragForce;
}

glm::vec3 Force::apply(float mass, const glm::vec3 & pos, const glm::vec3 & vel)
{
	return glm::vec3(0);
}

