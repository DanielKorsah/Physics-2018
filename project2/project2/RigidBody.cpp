// Math constants
#define _USE_MATH_DEFINES
#include <cmath>  
#include "RigidBody.h"


// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include "glm/ext.hpp"

void RigidBody::setInvInertia(const glm::mat3 &invInertia)
{
	float twelfthMass = 1.0f / 12.0f * Body::getMass();
	//get x scale value for width
	float w = Body::getScale()[0][0];
	//get the y scale value for height
	float h = Body::getScale()[1][1];
	//get the z scale value for depth
	float d = Body::getScale()[2][2];

	//set invInertia
	glm::mat3 inertia;
	inertia = {
		twelfthMass*(h*h + d*d), 0, 0,
		0, twelfthMass*(w*w + d*d), 0,
		0, 0, twelfthMass*(w*w + h*h)
	};

	m_invInertia = getRotate() * glm::inverse(inertia) * glm::transpose(getRotate());
}

glm::mat3 RigidBody::getinvInertia() 
{ 
	setInvInertia(m_invInertia);
	return m_invInertia; 
}


