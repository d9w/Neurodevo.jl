#include <iostream>
#include <exception>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "Viewer.h"
using namespace glm;

void Viewer::display() {

  glBegin(GL_LINES);
  glColor4f(1.0, 0.0, 0.0, 1.0);
  glVertex3f(-0.6, -0.6, -0.6);
  glVertex3f(0.6, -0.6, -0.6);
  glEnd();

  glBegin(GL_LINES);
  glColor4f(0.0, 1.0, 0.0, 1.0);
  glVertex3f(-0.6, -0.6, -0.6);
  glVertex3f(-0.6, 0.6, -0.6);
  glEnd();

  glBegin(GL_LINES);
  glColor4f(0.0, 0.0, 1.0, 1.0);
  glVertex3f(-0.6, -0.6, -0.6);
  glVertex3f(-0.6, -0.6, 0.6);
  glEnd();

}

void Viewer::KeyInput(){

	if (glfwGetKey( window, GLFW_KEY_UP ) == GLFW_PRESS)
    rotate_x += 5;
	if (glfwGetKey( window, GLFW_KEY_DOWN ) == GLFW_PRESS)
    rotate_x -= 5;
	if (glfwGetKey( window, GLFW_KEY_RIGHT ) == GLFW_PRESS)
    rotate_y += 5;
	if (glfwGetKey( window, GLFW_KEY_LEFT ) == GLFW_PRESS)
    rotate_y -= 5;
	if(glfwGetKey(window, GLFW_KEY_ESCAPE ) == GLFW_PRESS)
    glfwSetWindowShouldClose(window, 1);

  glRotatef( rotate_x, 1.0, 0.0, 0.0 );
  glRotatef( rotate_y, 0.0, 1.0, 0.0 );
}

Viewer::Viewer() {
  if( !glfwInit() ) {
    throw std::logic_error("Failed to initialize GLFW.");
  }

  window = glfwCreateWindow( 1024, 768, "nge", NULL, NULL);
  if( window == NULL ){
    glfwTerminate();
    throw std::logic_error("Failed to open GLFW window.");
  }
  glfwMakeContextCurrent(window);

	glewExperimental = true;
	if (glewInit() != GLEW_OK) {
		glfwTerminate();
		throw std::logic_error("Failed to initialize GLEW.");
	}

  glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
  glfwPollEvents();

  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glEnable(GL_CULL_FACE);
}

void Viewer::run() {

  do {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
		KeyInput();
    display();
    glFlush();
    glfwSwapBuffers(window);
    glfwPollEvents();
  }
	while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
         glfwWindowShouldClose(window) == 0 );

  glfwTerminate();
}
