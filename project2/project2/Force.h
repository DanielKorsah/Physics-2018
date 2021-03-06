#pragma once
#include <glm/glm.hpp>
#include <iostream>

class Body; //forward declaration to avoid circular dependencies

class Force
{
public:
	Force() {}
	~Force() {}

	virtual glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);
};


//Gravity class
class Gravity : public Force 
{
public:
	//constructors
	Gravity() {}
	Gravity(const glm::vec3 &gravity) { m_gravity = gravity; }

	//get and set
	glm::vec3 getGravity() const { return m_gravity; }
	void setGravity(glm::vec3 gravity) { m_gravity = gravity; }

	//physics
	////////////////fix???????????????????????????????
	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);

private:
	glm::vec3 m_gravity = glm::vec3(0.0f, -9.8f, 0.0);
};


//Drag Class
class Drag : public Force
{
public:
	Drag() {}

	//physics
	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);

private:
	float density = 0.001225f; //water is 0.846, air is 0.001225
	float coEff = 0.47f;
};

//Hooke force
class Hooke : public Force 
{
public:
	Hooke() {}
	Hooke(Body* b1, Body* b2, float ks, float kd, float rest) 
	{
		m_ks = ks;
		m_kd = kd;
		m_rest = rest;
		m_b1 = b1;
		m_b2 = b2;
	}

	//physics
	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel);

private:
	float m_ks; //spring stiffness
	float m_kd; //damping coefficient
	float m_rest; //spring rest length

	Body* m_b1; //pointer to body connected to one extremity of the spring
	Body* m_b2; //pointer to body connected to the other extremity of the spring
};


