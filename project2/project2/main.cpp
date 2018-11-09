#pragma once
// Math constants
#define _USE_MATH_DEFINES
#include <cmath>  
#include <random>

// Std. Includes
#include <string>
#include <time.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/matrix_operation.hpp>
#include "glm/ext.hpp"


// Other Libs
#include "SOIL2/SOIL2.h"
#include "windows.h"

// project includes
#include "Application.h"
#include "Shader.h"
#include "Mesh.h"
#include "Particle.h"
#include "Body.h"
#include "RigidBody.h"

Shader redShader;
Shader blueShader;
Shader transparent;

//forawrd declare demos
void Integration(Application app);
void BlowDryer(Application app);

// time
GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

//forward declare mode switching method
void CheckMode(Application app);

struct Cube
{
	glm::vec3 origin = glm::vec3(-5, 0.0f, -5);
	glm::vec3 bound = glm::vec3(5, 10.0f, 5);
};

struct Cone
{
	glm::vec3 origin = glm::vec3(0, -1.0f, 0);
	glm::vec3 base = glm::vec3(0, 3.0f, 0);
	float radius = 4.0f;
};

//initialise cube
Cube cube;
//init cone
Cone cone;

//toggle variable for moving cone
bool toggle = false;





bool inCone(Particle p, float* cm)
{
	//length of cone and related stuff
	glm::vec3 originToBase = cone.base - cone.origin;
	float coneHeight = glm::length(originToBase);
	glm::vec3 originToBaseNormal = glm::normalize(originToBase);
	
	//distance along axis of cone
	float distanceAlongCone = glm::dot(p.getPos() - cone.origin, originToBaseNormal);

	//if less than zero or farther along than the height of the cone, not in cone
	if (distanceAlongCone <= 0 || distanceAlongCone > coneHeight)
		return false;

	//getthe radius of the cone at the height of the point
	float coneRadiusAtPoint = (distanceAlongCone / coneHeight) * cone.radius;
	float distanceFromCentre = glm::length((p.getPos() - cone.origin) - distanceAlongCone * originToBaseNormal);

	//if if closer to the centre axis than the radius of the cone, must be inside cone
	if (distanceFromCentre < coneRadiusAtPoint)
	{
		*cm = ((1 - (distanceAlongCone/coneHeight)) * 30);
		*cm += (1 - distanceFromCentre/coneRadiusAtPoint) * 30;
		//std::cout << *cm << std::endl;
		return true;
	}
	else 
	{
		return false;
	}
		
}

bool isCorner(int row, int column, int squareSize)
{
	if (row == 0 && column == 0)
		return true;
	else if (row == 0 && column == squareSize - 1)
		return true;
	else if (row == squareSize - 1 && column == 0)
		return true;
	else if (row == squareSize - 1 && column == squareSize - 1)
		return true;
	else
		return false;
}





void upHooke(std::vector<std::vector<Particle*>> p_matrix, int i, int j, float spring, float damper, float rest)
{
	Hooke* hookeUp = new Hooke(p_matrix[i][j], p_matrix[i-1][j], spring, damper, rest);
	p_matrix[i][j]->addForce(hookeUp);
}

void downHooke(std::vector<std::vector<Particle*>> p_matrix, int i, int j, float spring, float damper, float rest)
{
	Hooke* hookeDown = new Hooke(p_matrix[i][j], p_matrix[i + 1][j], spring, damper, rest);
	p_matrix[i][j]->addForce(hookeDown);
}

void rightHooke(std::vector<std::vector<Particle*>> p_matrix, int i, int j, float spring, float damper, float rest)
{
	Hooke* hookeRight = new Hooke(p_matrix[i][j], p_matrix[i][j + 1], spring, damper, rest);
	p_matrix[i][j]->addForce(hookeRight);
}

void upRightHooke(std::vector<std::vector<Particle*>> p_matrix, int i, int j, float spring, float damper, float rest)
{
	Hooke* hookeUpRight = new Hooke(p_matrix[i][j], p_matrix[i - 1][j + 1], spring, damper, rest);
	p_matrix[i][j]->addForce(hookeUpRight);
}

void downRightHooke(std::vector<std::vector<Particle*>> p_matrix, int i, int j, float spring, float damper, float rest)
{
	Hooke* hookeDownRight = new Hooke(p_matrix[i][j], p_matrix[i + 1][j + 1], spring, damper, rest);
	p_matrix[i][j]->addForce(hookeDownRight);
}

void leftHooke(std::vector<std::vector<Particle*>> p_matrix, int i, int j, float spring, float damper, float rest)
{
	Hooke* hookeLeft = new Hooke(p_matrix[i][j], p_matrix[i][j - 1], spring, damper, rest);
	p_matrix[i][j]->addForce(hookeLeft);
}

void upLeftHooke(std::vector<std::vector<Particle*>> p_matrix, int i, int j, float spring, float damper, float rest)
{
	Hooke* hookeUpLeft = new Hooke(p_matrix[i][j], p_matrix[i - 1][j - 1], spring, damper, rest);
	p_matrix[i][j]->addForce(hookeUpLeft);
}

void downLeftHooke(std::vector<std::vector<Particle*>> p_matrix, int i, int j, float spring, float damper, float rest)
{
	Hooke* hookeDownLeft = new Hooke(p_matrix[i][j], p_matrix[i + 1][j - 1], spring, damper, rest);
	p_matrix[i][j]->addForce(hookeDownLeft);
}



void Integration(Application app)
{
	deltaTime = 0.0f;
	lastFrame = 0.0f;

	// create ground plane
	Mesh plane = Mesh::Mesh(Mesh::QUAD);
	// scale it up x5
	plane.scale(glm::vec3(5.0f, 5.0f, 5.0f));
	plane.setShader(Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag"));

	glm::vec3 gravity = glm::vec3(0.0f, -9.8f, 0.0f);

	// create particle
	std::vector<Particle> particles;
	int particleNum = 2;
	for (int i = 0; i < particleNum; i++)
	{
		Particle p = Particle::Particle();
		particles.push_back(p);
		//scale it down (x.1), translate it up by 2.5 and rotate it by 90 degrees around the x axis
		//particles[i].setPos(glm::vec3(0.0f, 4.0f, 0.0f));
		particles[i].scale(glm::vec3(0.5f, 0.5f, 0.5f));
		//particles[i].rotate((GLfloat) M_PI_2, glm::vec3(0.0f, 1.0f, 0.0f));
		particles[i].getMesh().setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));


		//initial velocty
		//particles[i].setVel(glm::vec3(sin(i), 0.0f, cos(i)));

		//make ring
		particles[i].setPos(glm::vec3(i + 1, 3, 0));
		//particles[i].setVel(glm::vec3(sin(i)*1.5f, .0f, cos(i)*1.5f));

		//set start acceleration to gravity
		particles[i].setAcc(gravity);
	}

	//height marker particle
	Particle m = Particle::Particle();
	m.setPos(glm::vec3(0, 3, 0));
	m.scale(glm::vec3(0.5f, 0.5f, 0.5f));
	m.getMesh().setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));


	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();

	//fixed timestep
	double physicsTime = 1.0f;
	const double fixedDeltaTime = 0.01f;
	double currentTime = (GLfloat)glfwGetTime();
	double accumulator = 0.0f;

	while (!glfwWindowShouldClose(app.getWindow()))
	{

		//fixed timstep
		double newTime = (GLfloat)glfwGetTime();
		double frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;

		while (accumulator >= fixedDeltaTime)
		{
			particles[0].setAcc((gravity));
			particles[1].setAcc(gravity);

			//Semi-Implicit Euler integration
			particles[0].getVel() += particles[0].getAcc() * fixedDeltaTime;
			particles[0].setPos(particles[0].getPos() + particles[0].getVel() * fixedDeltaTime);


			//forward Euler integration
			particles[1].setPos(particles[1].getPos() + particles[1].getVel() * fixedDeltaTime);
			particles[1].getVel() += particles[1].getAcc() * fixedDeltaTime;

			for (int i = 0; i < particleNum; i++)
			{


				//collisions to bound within the box
				for (int j = 0; j < 3; j++)
				{
					if (particles[i].getTranslate()[3][j] < cube.origin[j])
					{
						glm::vec3 diff = glm::vec3(0.0f);
						diff[j] = cube.origin[j] - particles[i].getPos()[j];
						particles[i].setPos(j, cube.origin[j] + diff[j]);
						particles[i].getVel()[j] *= -1.0f;
					}

					if (particles[i].getTranslate()[3][j] > cube.bound[j])
					{
						glm::vec3 diff = glm::vec3(0.0f);
						diff[j] = cube.bound[j] - particles[i].getPos()[j];
						particles[i].setPos(j, cube.bound[j] + diff[j]);
						particles[i].getVel()[j] *= -1.0f;
					}
				}

			}

			accumulator -= fixedDeltaTime;
			physicsTime += fixedDeltaTime;

		}

		// Set frame time
		GLfloat currentFrame = (GLfloat)glfwGetTime() - firstFrame;
		// the animation can be sped up or slowed down by multiplying currentFrame by a factor.
		currentFrame *= 1.5f;
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(deltaTime);

		//stwitch mode
		CheckMode(app);


		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
		// draw groud plane
		app.draw(plane);
		// draw particles
		for (int i = 0; i < particleNum; i++)
		{
			app.draw(particles[i].getMesh());
		}

		//render height marker
		app.draw(m.getMesh());

		app.display();
	}
	app.terminate();
}

void BlowDryer(Application app) 
{
	deltaTime = 0.0f;
	lastFrame = 0.0f;

	// create ground plane
	Mesh plane = Mesh::Mesh(Mesh::QUAD);
	// scale it up x5
	plane.scale(glm::vec3(5.0f, 5.0f, 5.0f));
	plane.setShader(Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag"));

	glm::vec3 gravity = glm::vec3(0.0f, -9.8f, 0.0f);


	// create particle
	std::vector<Particle> particles;
	int particleNum = 400;
	for (int i = 0; i < particleNum; i++)
	{
		Particle p = Particle::Particle();
		particles.push_back(p);
		//scale it down (x.1), translate it up by 2.5 and rotate it by 90 degrees around the x axis
		//particles[i].setPos(glm::vec3(0.0f, 4.0f, 0.0f));
		particles[i].scale(glm::vec3(0.5f, 0.5f, 0.5f));
		//particles[i].rotate((GLfloat) M_PI_2, glm::vec3(0.0f, 1.0f, 0.0f));
		particles[i].getMesh().setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));


		//initial velocty
		//particles[i].setVel(glm::vec3(sin(i), 0.0f, cos(i)));

		//make ring
		//particles[i].setPos(glm::vec3(sin(i), 3.0f, cos(i)));
		particles[i].setPos(glm::vec3(0, 5.0f, 0));
		particles[i].setVel(glm::vec3(sin(i)*1.5f, .0f, cos(i)*1.5f));

		//set start acceleration to gravity
		particles[i].setAcc(gravity);
		particles[i].setMass(1.0f);
	}

	//visualise the box
	Mesh cubeMesh = Mesh::Mesh("resources/models/cube.obj");
	cubeMesh.translate(glm::vec3(0.0f, 5.0f, 0.0f));
	cubeMesh.scale(glm::vec3(10.0f, 10.0f, 10.0f));
	cubeMesh.setShader(transparent);


	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();

	//fixed timestep
	double physicsTime = 1.0f;
	const double fixedDeltaTime = 0.01f;
	double currentTime = (GLfloat)glfwGetTime();
	double accumulator = 0.0f;

	

	// Game loop
	while (!glfwWindowShouldClose(app.getWindow()))
	{

		//fixed timstep
		double newTime = (GLfloat)glfwGetTime();
		double frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;

		while (accumulator >= fixedDeltaTime)
		{
			//move origin of cone
			if (toggle == true)
			{
				cone.origin = glm::vec3(sin(physicsTime - firstFrame) * 4, -0.5f, cos(physicsTime - firstFrame) * 4);
				cone.base = glm::vec3(sin(physicsTime - firstFrame) * 4, 3.0f, cos(physicsTime - firstFrame) * 4);
				cone.radius = 3.0f;
			}
			else 
			{ 
				cone.origin = glm::vec3(0, -1.0f, 0);
				cone.base = glm::vec3(0, 4.0f, 0);
				cone.radius = 4.0f;
			}

			for (int i = 0; i < particleNum; i++)
			{
				/*
				**	SIMULATION
				*/
				glm::vec3 gravityForce;
				gravityForce = gravity * particles[i].getMass();

				//check if in cone
				glm::vec3 coneForce = glm::vec3(0, 0, 0);
				float coneMultiplier = 1;
				

				if (inCone(particles[i], &coneMultiplier))
				{

					coneForce = (glm::vec3(glm::normalize(particles[i].getPos()-cone.origin) * coneMultiplier))*particles[i].getMass();
					particles[i].getMesh().setShader(redShader);
				}
				else 
				{
					coneForce = glm::vec3(0, 0, 0);
					particles[i].getMesh().setShader(blueShader);
				}

				
				//aero force
				float density = 0.001225f; //water is 0.846, air is 0.001225
				float velSquared = glm::length(particles[i].getVel()) * glm::length(particles[i].getVel());
				float coEff = 0.47f;
				glm::vec3 velNorm = glm::normalize(particles[i].getVel());
				glm::vec3 dragForce = 0.5f * density * velSquared * coEff * 1 * velNorm;


				//apply forces
				glm::vec3 acc = (gravityForce + coneForce - dragForce) / particles[i].getMass();
				particles[i].setAcc(acc);

				//Semi-Implicit Euler integration
				particles[i].translate(particles[i].getVel() * fixedDeltaTime);

				//collisions to bound within the box
				for (int j = 0; j < 3; j++)
				{
					if (particles[i].getTranslate()[3][j] < cube.origin[j])
					{
						glm::vec3 diff = glm::vec3(0.0f);
						diff[j] = cube.origin[j] - particles[i].getPos()[j];
						particles[i].setPos(j, cube.origin[j] + diff[j]);
						particles[i].getVel()[j] *= -0.8f;
					}

					if (particles[i].getTranslate()[3][j] > cube.bound[j])
					{
						glm::vec3 diff = glm::vec3(0.0f);
						diff[j] = cube.bound[j] - particles[i].getPos()[j];
						particles[i].setPos(j, cube.bound[j] + diff[j]);
						particles[i].getVel()[j] *= -0.8f;
					}
				}

				//set next frame velocity
				particles[i].setVel(particles[i].getVel() + particles[i].getAcc() * fixedDeltaTime);

			}

			accumulator -= fixedDeltaTime;
			physicsTime += fixedDeltaTime;
		}

		// Set frame time
		GLfloat currentFrame = (GLfloat)glfwGetTime() - firstFrame;
		// the animation can be sped up or slowed down by multiplying currentFrame by a factor.
		currentFrame *= 1.5f;
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(deltaTime);

		//switch mode
		CheckMode(app);

		//toggle cone movement
		if (app.keys[GLFW_KEY_8])
		{
			toggle = true;
		}

		if (app.keys[GLFW_KEY_9])
		{
			toggle = false;
		}


		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
		// draw groud plane
		app.draw(plane);
		
		// draw particles
		for (int i = 0; i < particleNum; i++)
		{
			app.draw(particles[i].getMesh());
		}

		//draw the cube
		app.draw(cubeMesh);

		app.display();
	}

	app.terminate();
}

void Rope(Application app) 
{
	
	float t = 0.0f;

	deltaTime = 0.0f;
	lastFrame = 0.0f;

	// create ground plane
	Mesh plane = Mesh::Mesh(Mesh::QUAD);
	// scale it up x5
	plane.scale(glm::vec3(5.0f, 5.0f, 5.0f));
	plane.setShader(Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag"));

	
	//create forces
	Gravity* g = new Gravity(glm::vec3(0.0f, -9.8f, 0.0f));
	Drag* d = new Drag();

	//Hooke parameter controls
	float spring = 50.0f;
	float damper = 30.0f;
	float rest = 0.5f;

	// create particles
	std::vector<Particle*> particles;
	int particleNum = 10;
	for (int i = 0; i < particleNum; i++)
	{
		Particle *p = new Particle();
		particles.push_back(p);
		//scale it down (x.1), translate it up by 2.5 and rotate it by 90 degrees around the x axis
		//particles[i].setPos(glm::vec3(0.0f, 4.0f, 0.0f));
		particles[i]->scale(glm::vec3(0.5f, 0.5f, 0.5f));
		//particles[i]->rotate((GLfloat) M_PI_2, glm::vec3(0.0f, 1.0f, 0.0f));
		particles[i]->getMesh().setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));


		//initial velocty
		//particles[i]->setVel(glm::vec3(sin(i), 0.0f, cos(i)));

		//make ring
		particles[i]->setPos(glm::vec3(sin(i), 6, cos(i)));
		//particles[i]->setVel(glm::vec3(sin(i)*1.5f, .0f, cos(i)*1.5f));
	}

	for (int i = 0; i < particleNum; i++)
	{
		//point 0 must get no forces
		if ((i > 0))
		{
			//add gravity, drag and hooke forces
			particles[i]->addForce(g);

			/*Drag* d = new Drag();
			particles[i]->addForce(d);*/
			Hooke* hooke1 = new Hooke(particles[i], particles[i - 1], spring, damper, rest);
			particles[i]->addForce(hooke1);
			if (i > 1)
			{
				Hooke* hooke2 = new Hooke(particles[i - 1], particles[i], spring, damper, rest);
				particles[i - 1]->addForce(hooke2);
			}
			if (i == particleNum - 2)
			{
				Hooke* hooke3 = new Hooke(particles[i], particles[i + 1], spring, damper, rest);
				particles[i]->addForce(hooke3);
			}
		}
	}



	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();

	//fixed timestep
	double physicsTime = 0.0f;
	const double fixedDeltaTime = 0.001f;
	double currentTime = (GLfloat)glfwGetTime();
	double accumulator = 0.0f;

	while (!glfwWindowShouldClose(app.getWindow()))
	{
		

		//fixed timstep
		double newTime = (GLfloat)glfwGetTime();
		double frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;

		while (accumulator >= fixedDeltaTime)
		{
			
			for (int i = 0; i < particleNum; i++)
			{
				particles[i]->setAcc(particles[i]->applyForces(particles[i]->getPos(), particles[i]->getVel(), physicsTime, fixedDeltaTime));
				if (particles[i]->getAcc() != particles[i]->getAcc())
				{
					std::cout << i << " is broke :  " << glm::to_string(particles[i]->getAcc()) << std::endl;
				}
			}

			for (int i = 0; i < particleNum; i++)
			{
				
				
				//Semi-Implicit Euler integration
				particles[i]->getVel() += particles[i]->getAcc() * fixedDeltaTime;
				particles[i]->translate(particles[i]->getVel() * fixedDeltaTime);

				//diagnosing
				//std::cout<< i << " "<< glm::to_string(particles[i]->getVel())<<std::endl;

				//collisions to bound within the box
				for (int j = 0; j < 3; j++)
				{
					if (particles[i]->getTranslate()[3][j] < cube.origin[j])
					{
						glm::vec3 diff = glm::vec3(0.0f);
						diff[j] = cube.origin[j] - particles[i]->getPos()[j];
						particles[i]->setPos(j, cube.origin[j] + diff[j]);
						particles[i]->getVel()[j] *= -0.8f;
					}

					if (particles[i]->getTranslate()[3][j] > cube.bound[j])
					{
						glm::vec3 diff = glm::vec3(0.0f);
						diff[j] = cube.bound[j] - particles[i]->getPos()[j];
						particles[i]->setPos(j, cube.bound[j] + diff[j]);
						particles[i]->getVel()[j] *= -0.8f;
					}
				}

			}

			accumulator -= fixedDeltaTime;
			physicsTime += fixedDeltaTime;

		}

		// Set frame time
		GLfloat currentFrame = (GLfloat)glfwGetTime() - firstFrame;
		// the animation can be sped up or slowed down by multiplying currentFrame by a factor.
		currentFrame *= 1.5f;
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(deltaTime);

		//toggle modes
		CheckMode(app);


		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
		// draw groud plane
		app.draw(plane);
		// draw particles
		for (int i = 0; i < particleNum; i++)
		{
			app.draw(particles[i]->getMesh());
		}


		app.display();

	}
	app.terminate();
}

void Chain(Application app)
{

	float t = 0.0f;

	deltaTime = 0.0f;
	lastFrame = 0.0f;

	// create ground plane
	Mesh plane = Mesh::Mesh(Mesh::QUAD);
	// scale it up x5
	plane.scale(glm::vec3(5.0f, 5.0f, 5.0f));
	plane.setShader(Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag"));


	//create forces
	Gravity* g = new Gravity(glm::vec3(0.0f, -9.8f, 0.0f));
	Drag* d = new Drag();

	//Hooke parameter controls
	float spring = 50.0f;
	float damper = 30.0f;
	float rest = 0.5f;

	// create particles
	std::vector<Particle*> particles;
	int particleNum = 10;
	for (int i = 0; i < particleNum; i++)
	{
		Particle *p = new Particle();
		particles.push_back(p);
		//scale it down (x.1), translate it up by 2.5 and rotate it by 90 degrees around the x axis
		//particles[i].setPos(glm::vec3(0.0f, 4.0f, 0.0f));
		particles[i]->scale(glm::vec3(0.5f, 0.5f, 0.5f));
		//particles[i]->rotate((GLfloat) M_PI_2, glm::vec3(0.0f, 1.0f, 0.0f));
		particles[i]->getMesh().setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));


		//initial velocty
		//particles[i]->setVel(glm::vec3(sin(i), 0.0f, cos(i)));

		//make ring
		particles[i]->setPos(glm::vec3(sin(i), 2, cos(i)));
		//particles[i]->setVel(glm::vec3(sin(i)*1.5f, .0f, cos(i)*1.5f));



	}

	for (int i = 0; i < particleNum; i++)
	{
		//point 0 must get no forces
		if ((i > 0) && (i != particleNum - 1))
		{
			//add gravity, drag and hooke forces
			particles[i]->addForce(g);

			/*Drag* d = new Drag();
			particles[i]->addForce(d);*/
			Hooke* hooke1 = new Hooke(particles[i], particles[i - 1], spring, damper, rest);
			particles[i]->addForce(hooke1);
			if (i > 1)
			{
				Hooke* hooke2 = new Hooke(particles[i - 1], particles[i], spring, damper, rest);
				particles[i - 1]->addForce(hooke2);
			}
			if (i == particleNum - 2)
			{
				Hooke* hooke3 = new Hooke(particles[i], particles[i + 1], spring, damper, rest);
				particles[i]->addForce(hooke3);
			}
		}
	}



	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();

	//fixed timestep
	double physicsTime = 0.0f;
	const double fixedDeltaTime = 0.001f;
	double currentTime = (GLfloat)glfwGetTime();
	double accumulator = 0.0f;

	while (!glfwWindowShouldClose(app.getWindow()))
	{


		//fixed timstep
		double newTime = (GLfloat)glfwGetTime();
		double frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;

		while (accumulator >= fixedDeltaTime)
		{

			for (int i = 0; i < particleNum; i++)
			{
				particles[i]->setAcc(particles[i]->applyForces(particles[i]->getPos(), particles[i]->getVel(), physicsTime, fixedDeltaTime));
				if (particles[i]->getAcc() != particles[i]->getAcc())
				{
					std::cout << i << " is broke :  " << glm::to_string(particles[i]->getAcc()) << std::endl;
				}
			}

			for (int i = 0; i < particleNum; i++)
			{


				//Semi-Implicit Euler integration
				particles[i]->getVel() += particles[i]->getAcc() * fixedDeltaTime;
				particles[i]->translate(particles[i]->getVel() * fixedDeltaTime);

				//diagnosing
				//std::cout<< i << " "<< glm::to_string(particles[i]->getVel())<<std::endl;

				//collisions to bound within the box
				for (int j = 0; j < 3; j++)
				{
					if (particles[i]->getTranslate()[3][j] < cube.origin[j])
					{
						glm::vec3 diff = glm::vec3(0.0f);
						diff[j] = cube.origin[j] - particles[i]->getPos()[j];
						particles[i]->setPos(j, cube.origin[j] + diff[j]);
						particles[i]->getVel()[j] *= -0.8f;
					}

					if (particles[i]->getTranslate()[3][j] > cube.bound[j])
					{
						glm::vec3 diff = glm::vec3(0.0f);
						diff[j] = cube.bound[j] - particles[i]->getPos()[j];
						particles[i]->setPos(j, cube.bound[j] + diff[j]);
						particles[i]->getVel()[j] *= -0.8f;
					}
				}

			}

			accumulator -= fixedDeltaTime;
			physicsTime += fixedDeltaTime;

		}

		// Set frame time
		GLfloat currentFrame = (GLfloat)glfwGetTime() - firstFrame;
		// the animation can be sped up or slowed down by multiplying currentFrame by a factor.
		currentFrame *= 1.5f;
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(deltaTime);

		//toggle modes
		CheckMode(app);


		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
		// draw groud plane
		app.draw(plane);
		// draw particles
		for (int i = 0; i < particleNum; i++)
		{
			app.draw(particles[i]->getMesh());
		}


		app.display();

	}
	app.terminate();
}

void Trampoline(Application app)
{
	float t = 0.0f;

	deltaTime = 0.0f;
	lastFrame = 0.0f;

	// create ground plane
	Mesh plane = Mesh::Mesh(Mesh::QUAD);
	// scale it up x5
	plane.scale(glm::vec3(5.0f, 5.0f, 5.0f));
	plane.setShader(Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag"));


	//create forces
	Gravity* g = new Gravity(glm::vec3(0.0f, -9.8f, 0.0f));
	Drag* d = new Drag();

	//Hooke parameter controls
	float spring = 70.0f;
	float damper = 65.0f;
	float rest = 0.5f;

	
	int clothSize = 10;

	//declare particle matrix as vector of vectors
	std::vector<std::vector<Particle*>> p_matrix;

	//**NOTE**: I'm indexing first by row then by column

	//Create the matrix of particles
	for (int i = 0; i < clothSize; i++)
	{
		//declare new row
		std::vector<Particle*> row;
		for (int j = 0; j < clothSize;  j++)
		{
			Particle *particle = new Particle();
			row.push_back(particle);
			row[j]->scale(glm::vec3(0.5f, 0.5f, 0.5f));
			row[j]->getMesh().setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));
			row[j]->setPos(glm::vec3(cube.origin.x + 0.5f + j , 2.5f, cube.origin.z + 0.5f + i ));
		}

		//add row to matrix
		p_matrix.push_back(row);
	}

	//set forces for each particle
	for (int i = 0; i < clothSize; i++)
	{
		for (int j = 0; j < clothSize; j++)
		{
			//give all but the corners forces
			if (!isCorner(i, j, clothSize))
			{
				//gravity
				p_matrix[i][j]->addForce(g);

				//hookes for top edge
				if (i == 0)
				{	
					//none with up
					downHooke(p_matrix, i, j, spring, damper, rest);
					leftHooke(p_matrix, i, j, spring, damper, rest);
					rightHooke(p_matrix, i, j, spring, damper, rest);

					downLeftHooke(p_matrix, i, j, spring, damper, rest);
					downRightHooke(p_matrix, i, j, spring, damper, rest);
				}
				else if (i == clothSize - 1)
				{
					//hookes for bottom edge

					//none with down component
					upHooke(p_matrix, i, j, spring, damper, rest);
					leftHooke(p_matrix, i, j, spring, damper, rest);
					rightHooke(p_matrix, i, j, spring, damper, rest);

					upLeftHooke(p_matrix, i, j, spring, damper, rest);
					upRightHooke(p_matrix, i, j, spring, damper, rest);
				}
				else if (j == 0)
				{
					//hookes for left edge

					//none that have an left component
					upHooke(p_matrix, i, j, spring, damper, rest);
					downHooke(p_matrix, i, j, spring, damper, rest);
					rightHooke(p_matrix, i, j, spring, damper, rest);

					upRightHooke(p_matrix, i, j, spring, damper, rest);
					downRightHooke(p_matrix, i, j, spring, damper, rest);
				}
				else if (j == clothSize - 1)
				{
					//hookes for right edge

					//none with right compnent
					upHooke(p_matrix, i, j, spring, damper, rest);
					downHooke(p_matrix, i, j, spring, damper, rest);
					leftHooke(p_matrix, i, j, spring, damper, rest);

					upLeftHooke(p_matrix, i, j, spring, damper, rest);
					downLeftHooke(p_matrix, i, j, spring, damper, rest);
				}
				else
				{
					//all other hookes

					//all types
					upHooke(p_matrix, i, j, spring, damper, rest);
					downHooke(p_matrix, i, j, spring, damper, rest);
					leftHooke(p_matrix, i, j, spring, damper, rest);
					rightHooke(p_matrix, i, j, spring, damper, rest);

					upLeftHooke(p_matrix, i, j, spring, damper, rest);
					upRightHooke(p_matrix, i, j, spring, damper, rest);
					downLeftHooke(p_matrix, i, j, spring, damper, rest);
					downRightHooke(p_matrix, i, j, spring, damper, rest);

				}
			}
		}

	}

	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();

	//fixed timestep
	double physicsTime = 0.0f;
	const double fixedDeltaTime = 0.001f;
	double currentTime = (GLfloat)glfwGetTime();
	double accumulator = 0.0f;

	while (!glfwWindowShouldClose(app.getWindow()))
	{


		//fixed timstep
		double newTime = (GLfloat)glfwGetTime();
		double frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;

		while (accumulator >= fixedDeltaTime)
		{

			for (int i = 0; i < clothSize; i++)
			{
				for (int j = 0; j < clothSize; j++)
				{
					p_matrix[i][j]->setAcc(p_matrix[i][j]->applyForces(p_matrix[i][j]->getPos(), p_matrix[i][j]->getVel(), physicsTime, fixedDeltaTime));
				}
				
			}

			for (int i = 0; i < clothSize; i++)
			{
				for (int j = 0; j < clothSize; j++)
				{
					//Semi-Implicit Euler integration
					p_matrix[i][j]->getVel() += p_matrix[i][j]->getAcc() * fixedDeltaTime;
					p_matrix[i][j]->translate(p_matrix[i][j]->getVel() * fixedDeltaTime);

				}
			}

			for (int i = 0; i < clothSize; i++)
			{
				for (int j = 0; j < clothSize; j++)
				{


					//Semi-Implicit Euler integration
					p_matrix[i][j]->getVel() += p_matrix[i][j]->getAcc() * fixedDeltaTime;
					p_matrix[i][j]->translate(p_matrix[i][j]->getVel() * fixedDeltaTime);



					//collisions to bound within the box
					for (int k = 0; k < 3; k++)
					{
						if (p_matrix[i][j]->getTranslate()[3][k] < cube.origin[k])
						{
							glm::vec3 diff = glm::vec3(0.0f);
							diff[k] = cube.origin[k] - p_matrix[i][j]->getPos()[k];
							p_matrix[i][j]->setPos(k, cube.origin[k] + diff[k]);
							p_matrix[i][j]->getVel()[k] *= -0.8f;
						}

						if (p_matrix[i][j]->getTranslate()[3][k] > cube.bound[k])
						{
							glm::vec3 diff = glm::vec3(0.0f);
							diff[j] = cube.bound[j] - p_matrix[i][k]->getPos()[k];
							p_matrix[i][j]->setPos(k, cube.bound[k] + diff[k]);
							p_matrix[i][j]->getVel()[k] *= -0.8f;
						}
					}
				}
			}

			accumulator -= fixedDeltaTime;
			physicsTime += fixedDeltaTime;

		}

		// Set frame time
		GLfloat currentFrame = (GLfloat)glfwGetTime() - firstFrame;
		// the animation can be sped up or slowed down by multiplying currentFrame by a factor.
		currentFrame *= 1.5f;
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(deltaTime);

		//toggle modes
		CheckMode(app);


		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
		// draw groud plane
		app.draw(plane);
		// draw particles
		for (int i = 0; i < clothSize; i++)
		{
			for (int j = 0; j < clothSize; j++)
			{
				app.draw(p_matrix[i][j]->getMesh());
			}
		}


		app.display();

	}
	app.terminate();
}

void Flag(Application app)
{
	float t = 0.0f;

	deltaTime = 0.0f;
	lastFrame = 0.0f;

	// create ground plane
	Mesh plane = Mesh::Mesh(Mesh::QUAD);
	// scale it up x5
	plane.scale(glm::vec3(5.0f, 5.0f, 5.0f));
	plane.setShader(Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag"));


	//create forces
	Gravity* g = new Gravity(glm::vec3(0.0f, -9.8f, 0.0f));
	Drag* d = new Drag();

	//Hooke parameter controls
	float spring = 200.0f;
	float damper = 10.0f;
	float rest = 1.0f;

	float diagonal_spring = 00.0f;
	float diagonal_damper = 10.0f;
	float diagonal_rest = 2.8f;


	int clothSize = 10;

	//declare particle matrix as vector of vectors
	std::vector<std::vector<Particle*>> p_matrix;

	//**NOTE**: I'm indexing first by row then by column

	//Create the matrix of particles
	for (int i = 0; i < clothSize; i++)
	{
		//declare new row
		std::vector<Particle*> row;
		for (int j = 0; j < clothSize; j++)
		{
			Particle *particle = new Particle();
			row.push_back(particle);
			row[j]->scale(glm::vec3(0.5f, 0.5f, 0.5f));
			row[j]->getMesh().setShader(Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag"));
			row[j]->setPos(glm::vec3(cube.origin.x + 0.5f + j, 8.0f, cube.origin.z + 5.0f + i));
		}

		//add row to matrix
		p_matrix.push_back(row);
	}

	//set forces for each particle
	for (int i = 0; i < clothSize; i++)
	{
		for (int j = 0; j < clothSize; j++)
		{
			//give all but the corners forces
			if (!(i == 0 && j == 0) && !(i == 0 && j == clothSize - 1))
			{
				//gravity
				p_matrix[i][j]->addForce(g);

				//hookes for top edge
				if (i == 0)
				{
					//none with up
					downHooke(p_matrix, i, j, spring, damper, rest);
					leftHooke(p_matrix, i, j, spring, damper, rest);
					rightHooke(p_matrix, i, j, spring, damper, rest);

					downLeftHooke(p_matrix, i, j, diagonal_spring, damper, rest);
					downRightHooke(p_matrix, i, j, diagonal_spring, damper, rest);
				}
				else if (i == clothSize - 1)
				{
					//hookes for bottom edge

					//none with down component
					upHooke(p_matrix, i, j, spring, damper, rest);
					

					if (j != clothSize - 1)
					{
						rightHooke(p_matrix, i, j, spring, damper, rest);
						upRightHooke(p_matrix, i, j, diagonal_spring, damper, rest);
					}

					if (j != 0)
					{
						leftHooke(p_matrix, i, j, spring, damper, rest);
						upLeftHooke(p_matrix, i, j, diagonal_spring, damper, rest);
					}
					
				}
				else if (j == 0)
				{
					//hookes for left edge

					//none that have an left component
					upHooke(p_matrix, i, j, spring, damper, rest);
					downHooke(p_matrix, i, j, spring, damper, rest);
					rightHooke(p_matrix, i, j, spring, damper, rest);

					upRightHooke(p_matrix, i, j, diagonal_spring, damper, rest);
					downRightHooke(p_matrix, i, j, diagonal_spring, damper, rest);
				}
				else if (j == clothSize - 1)
				{
					//hookes for right edge

					//none with right compnent
					upHooke(p_matrix, i, j, spring, damper, rest);
					downHooke(p_matrix, i, j, spring, damper, rest);
					leftHooke(p_matrix, i, j, spring, damper, rest);

					upLeftHooke(p_matrix, i, j, diagonal_spring, damper, rest);
					downLeftHooke(p_matrix, i, j, diagonal_spring, damper, rest);
				}
				else
				{
					//all other hookes

					//all types
					upHooke(p_matrix, i, j, spring, damper, rest);
					downHooke(p_matrix, i, j, spring, damper, rest);
					leftHooke(p_matrix, i, j, spring, damper, rest);
					rightHooke(p_matrix, i, j, spring, damper, rest);

					upLeftHooke(p_matrix, i, j, diagonal_spring, damper, rest);
					upRightHooke(p_matrix, i, j, diagonal_spring, damper, rest);
					downLeftHooke(p_matrix, i, j, diagonal_spring, damper, rest);
					downRightHooke(p_matrix, i, j, diagonal_spring, damper, rest);
				}


			}
		}

	}

	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();

	//fixed timestep
	double physicsTime = 0.0f;
	const double fixedDeltaTime = 0.001f;
	double currentTime = (GLfloat)glfwGetTime();
	double accumulator = 0.0f;

	while (!glfwWindowShouldClose(app.getWindow()))
	{


		//fixed timstep
		double newTime = (GLfloat)glfwGetTime();
		double frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;

		while (accumulator >= fixedDeltaTime)
		{

			for (int i = 0; i < clothSize; i++)
			{
				for (int j = 0; j < clothSize; j++)
				{
					p_matrix[i][j]->setAcc(p_matrix[i][j]->applyForces(p_matrix[i][j]->getPos(), p_matrix[i][j]->getVel(), physicsTime, fixedDeltaTime));
				}
			}

			
			for (int i = 0; i < clothSize; i++)
			{
				for (int j = 0; j < clothSize; j++)
				{
					//Semi-Implicit Euler integration
					p_matrix[i][j]->getVel() += p_matrix[i][j]->getAcc() * fixedDeltaTime;
					p_matrix[i][j]->translate(p_matrix[i][j]->getVel() * fixedDeltaTime);

				}
			}

			for (int i = 0; i < clothSize; i++)
			{
				for (int j = 0; j < clothSize; j++)
				{


					//Semi-Implicit Euler integration
					p_matrix[i][j]->getVel() += p_matrix[i][j]->getAcc() * fixedDeltaTime;
					p_matrix[i][j]->translate(p_matrix[i][j]->getVel() * fixedDeltaTime);



					//collisions to bound to above plane
					if (p_matrix[i][j]->getPos().y < 0)
					{
							
						p_matrix[i][j]->setPos(glm::vec3(p_matrix[i][j]->getPos().x , 0, p_matrix[i][j]->getPos().z));
						p_matrix[i][j]->getVel().y *= -0.2f;
					}
						
				}
			}

			accumulator -= fixedDeltaTime;
			physicsTime += fixedDeltaTime;

		}

		// Set frame time
		GLfloat currentFrame = (GLfloat)glfwGetTime() - firstFrame;
		// the animation can be sped up or slowed down by multiplying currentFrame by a factor.
		currentFrame *= 1.5f;
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(deltaTime);

		//toggle modes
		CheckMode(app);

		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
		// draw groud plane
		app.draw(plane);
		// draw particles
		for (int i = 0; i < clothSize; i++)
		{
			for (int j = 0; j < clothSize; j++)
			{
				app.draw(p_matrix[i][j]->getMesh());
			}
		}


		app.display();

	}
	app.terminate();
}

void RigidBody1(Application app)
{
	deltaTime = 0.0f;
	lastFrame = 0.0f;

	// create ground plane
	Mesh plane = Mesh::Mesh(Mesh::QUAD);
	// scale it up x5
	plane.scale(glm::vec3(5.0f, 5.0f, 5.0f));
	plane.setShader(Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag"));

	glm::vec3 gravity = glm::vec3(0.0f, -9.8f, 0.0f);

	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();

	//fixed timestep
	double physicsTime = 1.0f;
	const double fixedDeltaTime = 0.01f;
	double currentTime = (GLfloat)glfwGetTime();
	double accumulator = 0.0f;

	//set up cubic rigidbody
	RigidBody rb = RigidBody();
	Mesh m = Mesh::Mesh(Mesh::CUBE);
	rb.setMesh(m);
	Shader rbShader = Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag");
	rb.getMesh().setShader(rbShader);
	rb.scale(glm::vec3(1.0f, 1.0f, 1.0f));

	//rigid body motion values
	rb.translate(glm::vec3(0.0f, 4.0f, 0.0f));
	rb.setVel(glm::vec3(0.0f, 0.0f, 0.0f));
	rb.setAngVel(glm::vec3(0.0f, 2.0f, 0.0f));
	Gravity* g = new Gravity(glm::vec3(0.0f, -9.8f, 0.0f));

	rb.addForce(g);

	// Game loop
	while (!glfwWindowShouldClose(app.getWindow()))
	{

		//fixed timstep
		double newTime = (GLfloat)glfwGetTime();
		double frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;


		while (accumulator >= fixedDeltaTime)
		{

			rb.setAcc(rb.applyForces(rb.Body::getPos(), rb.Body::getVel(), physicsTime, fixedDeltaTime));

			//Semi - Implicit Euler integration
			rb.Body::getVel() += rb.Body::getAcc() * fixedDeltaTime;
			rb.translate(rb.Body::getVel() * fixedDeltaTime);

			//rotation integration
			rb.setAngVel(rb.getAngVel() + fixedDeltaTime * rb.getAngAcc());
			//create skew symetric matrix for w
			glm::mat3 angVelSkew = glm::matrixCross3(rb.getAngVel());
			//create 3x3 rotation matrix from rb rotation matrix
			glm::mat3 R = glm::mat3(rb.getRotate());
			//update rotation matrix
			R += fixedDeltaTime * angVelSkew * R;
			R = glm::orthonormalize(R);
			rb.setRotate(glm::mat4(R));


			//collision
			//scale[1][1] is y scale
			if (rb.getPos().y <= plane.getPos().y + rb.getScale()[1][1])
			{
				rb.Body::setPos(1, plane.getPos().y + rb.getScale()[1][1]);
			}

			accumulator -= fixedDeltaTime;
			physicsTime += fixedDeltaTime;
		}

		// Set frame time
		GLfloat currentFrame = (GLfloat)glfwGetTime() - firstFrame;
		// the animation can be sped up or slowed down by multiplying currentFrame by a factor.
		currentFrame *= 1.5f;
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(deltaTime);

		//switch mode
		CheckMode(app);

		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
		// draw groud plane
		app.draw(plane);
		app.draw(rb.getMesh());

		app.display();
	}

	app.terminate();
}


void RigidBody2(Application app)
{
	deltaTime = 0.0f;
	lastFrame = 0.0f;

	// create ground plane
	Mesh plane = Mesh::Mesh(Mesh::QUAD);
	// scale it up x5
	plane.scale(glm::vec3(5.0f, 5.0f, 5.0f));
	plane.setShader(Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag"));

	glm::vec3 gravity = glm::vec3(0.0f, -9.8f, 0.0f);

	// time
	GLfloat firstFrame = (GLfloat)glfwGetTime();

	//fixed timestep
	double physicsTime = 1.0f;
	const double fixedDeltaTime = 0.01f;
	double currentTime = (GLfloat)glfwGetTime();
	double accumulator = 0.0f;

	//set up cubic rigidbody
	RigidBody rb = RigidBody();
	Mesh m = Mesh::Mesh(Mesh::CUBE);
	rb.setMesh(m);
	Shader rbShader = Shader("resources/shaders/physics.vert", "resources/shaders/physics.frag");
	rb.getMesh().setShader(rbShader);
	rb.scale(glm::vec3(2.0f, 1.0f, 1.0f));
	rb.setMass(1.0f);

	//rigid body motion values
	rb.translate(glm::vec3(0.0f, 4.0f, 0.0f));
	rb.setVel(glm::vec3(0.0f, 0.0f, 0.0f));
	rb.setAngVel(glm::vec3(0.0f, 0.0f, 2.0f));
	Gravity* g = new Gravity(glm::vec3(0.0f, -9.8f, 0.0f));

	rb.addForce(g);

	// Game loop
	while (!glfwWindowShouldClose(app.getWindow()))
	{

		//fixed timstep
		double newTime = (GLfloat)glfwGetTime();
		double frameTime = newTime - currentTime;
		currentTime = newTime;

		accumulator += frameTime;


		while (accumulator >= fixedDeltaTime)
		{

			rb.setAcc(rb.applyForces(rb.Body::getPos(), rb.Body::getVel(), physicsTime, fixedDeltaTime));

			//Semi - Implicit Euler integration
			rb.Body::getVel() += rb.Body::getAcc() * fixedDeltaTime;
			rb.translate(rb.Body::getVel() * fixedDeltaTime);

			//rotation integration

			//set w
			rb.setAngVel(rb.getAngVel() + fixedDeltaTime * rb.getAngAcc());
			//create skew symetric matrix for w
			glm::mat3 angVelSkew = glm::matrixCross3(rb.getAngVel());
			//create 3x3 rotation matrix from rb rotation matrix
			glm::mat3 R = glm::mat3(rb.getRotate());
			//update rotation matrix
			R += fixedDeltaTime * angVelSkew * R;
			R = glm::orthonormalize(R);
			rb.setRotate(glm::mat4(R));

			

			//collision

			std::vector<glm::vec3> collisionPoints = {};
			std::vector<Vertex> vertices = rb.getMesh().Mesh::getVertices();
			glm::vec3 vertexSum;
			//for each vertex of the rigidbody, if it's below the plane add it to a vector of collision points
			//at the same time calculate centre of mass
			for (Vertex v : vertices)
			{
				vertexSum += v.getCoord();
				glm::mat3 m = glm::mat3(rb.getMesh().getModel());
				glm::vec3 worldSpaceVertex = m * v.getCoord() + rb.getPos();
				//glm::mat4 worldPoint = v.getCoord() * modelMat;
				if (worldSpaceVertex.y < plane.getPos().y)
				{
					collisionPoints.push_back(worldSpaceVertex);
				}
			}
			//set centre of mass
			glm::vec3 centreOfMass = vertexSum / vertices.size();

			if (!collisionPoints.empty())
			{
				//get average collision point
				glm::vec3 colPoint;
				if (collisionPoints.size() > 1)
				{
					glm::vec3 sumVec;
					for (glm::vec3 p : collisionPoints)
					{
						sumVec += p;
					}

					glm::vec3 avgVec = sumVec / collisionPoints.size();

					colPoint = avgVec;
				}
				else
				{
					colPoint = collisionPoints[0];
				}

				//get collision normal (straight down for plane collision only)
				glm::vec3 colNormal;
				colNormal = glm::vec3(0, -1, 0);
				colNormal = glm::normalize(colNormal);

				//resolve overlap 
				glm::vec3 colOverlap = glm::vec3(colPoint.x, 0, colPoint.z) - colPoint;
				rb.setPos(rb.getPos() + colOverlap);

				

				//generate the impulse
				//e = coefficient of restitution
				float e = rb.getRestitution();
				//vr = v for plane collision only
				glm::vec3 vr = rb.getVel();
				//get distance r from CoM to point: c to r = r-c
				glm::vec3 r = colPoint - centreOfMass;
				glm::vec3 numerator = -(1 + e) * vr - colNormal;
				//break it down a bit to rduce nested glm functions and brackets
				glm::vec3 rCrossN = glm::cross(r, colNormal);
				float denominator = glm::pow(rb.getMass(), -1.0f) + glm::dot(colNormal, (rb.getinvInertia() * glm::cross(rCrossN, rb.getVel())));

				glm::vec3 impulse = numerator / denominator;

				//set new velocity and rotation using impulse
				rb.Body::setVel(rb.getVel() - glm::dot(impulse/rb.getMass(), colNormal));

				//set new angular velocity using impulse
				rb.setAngVel(rb.getAngVel() - impulse * rb.getinvInertia() * (glm::cross(r, colNormal)));

			}

			accumulator -= fixedDeltaTime;
			physicsTime += fixedDeltaTime;
		}

		// Set frame time
		GLfloat currentFrame = (GLfloat)glfwGetTime() - firstFrame;
		// the animation can be sped up or slowed down by multiplying currentFrame by a factor.
		currentFrame *= 1.5f;
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		/*
		**	INTERACTION
		*/
		// Manage interaction
		app.doMovement(deltaTime);

		//switch mode
		CheckMode(app);

		/*
		**	RENDER
		*/
		// clear buffer
		app.clear();
		// draw groud plane
		app.draw(plane);
		app.draw(rb.getMesh());

		app.display();
	}

	app.terminate();
}

//demo switching method
void CheckMode(Application app)
{
	glfwPollEvents();
	if (app.keys[GLFW_KEY_1])
	{
		app.clear();
		Integration(app);
	}

	if (app.keys[GLFW_KEY_2])
	{
		app.clear();;
		BlowDryer(app);
	}

	if (app.keys[GLFW_KEY_3])
	{
		app.clear();
		Rope(app);
	}

	if (app.keys[GLFW_KEY_4])
	{
		app.clear();
		Chain(app);
	}

	if (app.keys[GLFW_KEY_5])
	{
		app.clear();
		Trampoline(app);
	}

	if (app.keys[GLFW_KEY_6])
	{
		app.clear();
		Flag(app);
	}
}

// main function
int main()
{
	// create application
	Application app = Application::Application();
	app.initRender();
	Application::camera.setCameraPosition(glm::vec3(0.0f, 4.0f, 15.0f));

	redShader = Shader("resources/shaders/solid.vert", "resources/shaders/solid_red.frag");
	blueShader = Shader("resources/shaders/solid.vert", "resources/shaders/solid_blue.frag");
	transparent = Shader("resources/shaders/solid.vert", "resources/shaders/solid_transparent.frag");

	//start from the beginning
	//Integration(app);

	//mass spring shortcut
	//Rope(app);

	//rigidbody shortcut
	RigidBody2(app);

	
	return EXIT_SUCCESS;
}

