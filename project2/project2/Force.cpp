#include "Force.h"
#include "Body.h"


Force::Force()
{
}

S
Force::~Force()
{
}


class Gravity : public Force
{
public:
	//constructors
	Gravity() {}
	Gravity(const glm::vec3 &gravity) { m_gravity = gravity; }

	//physics
	glm::vec3 apply(float mass, const glm::vec3 &pos, const glm::vec3 &vel)
	{

	}

private:
	glm::vec3 m_gravity = glm::vec3(0.0f, -9.8f, 0.0);
};

