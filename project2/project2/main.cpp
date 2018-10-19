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

// project includes
#include "Application.h"
#include "Shader.h"
#include "Mesh.h"
#include "Particle.h"
#include "Body.h"

Shader redShader;
Shader blueShader;
Shader transparent;

//forawrd declare demos
void Integration(Application app);
void BlowDryer(Application app);

// time
GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

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

bool inCone(Particle p, float *cm)
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
		glfwPollEvents();
		if (app.keys[GLFW_KEY_2])
		{
			app.clear();
			Integration(app);
		}

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
		particles[i].setPos(glm::vec3(i+1, 3, 0));
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
		glfwPollEvents();
		if (app.keys[GLFW_KEY_1])
		{
			app.clear();
			BlowDryer(app);
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

		//render height marker
		app.draw(m.getMesh());

		app.display();
	}
	app.terminate();
}

void MassSpring(Application app) 
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
	float spring = 8.0f;
	float damper = 2.0f;
	float rest = 0.5f;

	// create particles
	std::vector<Particle> particles;
	int particleNum = 3;
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
		particles[i].setPos(glm::vec3(i + 1, 6, 0));
		//particles[i].setVel(glm::vec3(sin(i)*1.5f, .0f, cos(i)*1.5f));

		
		//point 0 must get no forces
		if (i > 0) 
		{
			//add gravity, drag and hooke forcesparticles[i].addForce(g);
			particles[i].addForce(g);
			Hooke* hooke = new Hooke(&particles[i], &particles[i-1], spring, damper, rest);
			particles[i].addForce(hooke);
			if (i > 1)
			{
				particles[i-1].addForce(hooke);
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
				
				particles[i].setAcc(particles[i].applyForces(particles[i].getPos(), particles[i].getVel(), physicsTime, fixedDeltaTime));
				/*if (i == 0)
					particles[i].setAcc(glm::vec3(0));*/
				//Semi-Implicit Euler integration
				particles[i].getVel() += particles[i].getAcc() * fixedDeltaTime;
				particles[i].translate(particles[i].getVel() * fixedDeltaTime);

				//collisions to bound within the box
				/*for (int j = 0; j < 3; j++)
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
				}*/

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


		app.display();
	}
	app.terminate();
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

	//coursework 1, particles
	//Integration(app);

	//use force class
	MassSpring(app);

	
	return EXIT_SUCCESS;
}

