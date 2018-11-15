#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "Mesh.h"
#include "Force.h"

class Body
{
public:
	Body();
	~Body();

	/*
	** GET METHODS
	*/
	// mesh
	Mesh& getMesh() { return m_mesh; }

	// transform matrices
	glm::mat4 getTranslate() const { return m_mesh.getTranslate(); }
	glm::mat3 getRotate() const { return m_mesh.getRotate(); }
	glm::mat3 getScale() const { return m_mesh.getScale(); }
	glm::mat3 getprevRotate() { return m_prevRotate; };

	//set rotation matrix
	void setRotate(const glm::mat4 &mat) { m_mesh.setRotate(mat); }

	// dynamic variables
	glm::vec3& getAcc() { return m_acc; }
	glm::vec3& getVel() { return m_vel; }
	glm::vec3& getPos() { return m_pos; }


	// physical properties
	float getMass() const { return m_mass; }
	float getCor() { return m_cor; }

	/*
	** SET METHODS
	*/
	// mesh
	void setMesh(Mesh m) { m_mesh = m; }

	// dynamic variables
	void setAcc(const glm::vec3 &vect) { m_acc = vect; }
	void setVel(const glm::vec3 &vect) { m_vel = vect; }
	void setVel(int i, float v) { m_vel[i] = v; } //set the ith coordinate of the velocity vector
	void setPos(const glm::vec3 &vect) { m_pos = vect; m_mesh.setPos(vect); }
	void setPos(int i, float p) { m_pos[i] = p; m_mesh.setPos(i, p); } //set the ith coordinate of the position vector
	void setPrevRotate(glm::mat4 r) { m_prevRotate = r; }

	// physical properties
	void setCor(float cor) { m_cor = cor; }
	void setMass(float mass) { m_mass = mass; }

	/*
	** OTHER METHODS
	*/

	// transformation methods
	void translate(const glm::vec3 &vect);
	void rotate(float angle, const glm::vec3 &vect);
	void scale(const glm::vec3 &vect);


	//apply forces
	glm::vec3 applyForces(glm::vec3 x, glm::vec3 v, float t, float dt)
	{
		glm::vec3 fAccumulator = glm::vec3(0.0f);
		for (auto &f : m_forces)
		{
			fAccumulator += f->apply(getMass(), x, v);
		}
		//guard against nans
		if (getMass() == 0) { setMass(1.0f); }
		return fAccumulator / getMass();
	}
	void addForce(Force* f) { m_forces.push_back(f); }

	glm::vec3 sumOfForces(glm::vec3 x, glm::vec3 v)
	{
		glm::vec3 fAccumulator = glm::vec3(0.0f);
		for (auto &f : m_forces)
		{
			fAccumulator += f->apply(getMass(), x, v);
		}
		return fAccumulator;
	}
	
	
private:

	Mesh m_mesh; // mesh used to represent the body

	float m_mass; // mass
	float m_cor; // coefficient of restitution

	glm::mat3 m_prevRotate;
	glm::vec3 m_acc; // acceleration
	glm::vec3 m_vel; // velocity
	glm::vec3 m_pos; // position

	std::vector<Force*> m_forces;
	

};

