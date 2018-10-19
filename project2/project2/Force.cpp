#include <iostream>
#include <cmath>
#include "Force.h"
#include "Body.h"
#include <iostream>

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

glm::vec3 Hooke::apply(float mass, const glm::vec3 & pos, const glm::vec3 & vel)
{
	glm::vec3 direction;
	direction = m_b2->getPos() - m_b1->getPos();

	float length = glm::length(direction);	//length of spring
	direction = direction / length;	//normalise direction

	//velocities used
	float v1 = glm::dot(m_b1->getVel(), direction);
	float v2 = glm::dot(m_b2->getVel(), direction);

	//spring and damper components of formula
	glm::vec3 spring(-m_ks * (m_rest - length));
	glm::vec3 damper(-m_kd * (v1 - v2));

	glm::vec3 hookeForce(spring + damper);
	if (hookeForce * direction != hookeForce * direction)
		std::cout << std::endl;
	return hookeForce * direction;
}
