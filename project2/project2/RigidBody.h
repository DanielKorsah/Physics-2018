#pragma once
#include "Body.h"

class RigidBody : public Body
{
public:
	RigidBody();
	~RigidBody();

	//setters and getters
	void setAngVel(const glm::vec3 &omega) { m_angVel = omega; }
	void setAngAcc(const glm::vec3 &alpha) { m_angAcc = alpha; }
	void setInvIntertia(const glm::mat3 &invInertia) { m_invInertia; }

	glm::vec3 getAngVel() { return m_angVel; }
	glm::vec3 getAngAcc() { return m_angAcc; }
	glm::mat3 getinvInertia() { return m_invInertia; }
	void scale(glm::vec3 vect);

private:
	float m_density;
	glm::mat3 m_invInertia;	//inverse inertia
	glm::vec3 m_angVel;	//angular velocity
	glm::vec3 m_angAcc;	//angular acceleration
};