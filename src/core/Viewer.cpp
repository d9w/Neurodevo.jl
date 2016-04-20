#include <iostream>
#include <exception>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "Viewer.h"
#include "ANN.h"
using namespace glm;

double Viewer::ce(int position, int dim) {
  double scale;
  if (dim >= 3) {
    scale = 50.0;
  } else {
    scale = (double)lengths[dim];
  }
  return position/scale-0.5;
}

// taken straight from GLUT
void Viewer::DrawBox(GLfloat size, GLenum type) {
  static GLfloat n[6][3] =
  {
    {-1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {1.0, 0.0, 0.0},
    {0.0, -1.0, 0.0},
    {0.0, 0.0, 1.0},
    {0.0, 0.0, -1.0}
  };
  static GLint faces[6][4] =
  {
    {0, 1, 2, 3},
    {3, 2, 6, 7},
    {7, 6, 5, 4},
    {4, 5, 1, 0},
    {5, 6, 2, 1},
    {7, 4, 0, 3}
  };
  GLfloat v[8][3];
  GLint i;

  v[0][0] = v[1][0] = v[2][0] = v[3][0] = -size / 2;
  v[4][0] = v[5][0] = v[6][0] = v[7][0] = size / 2;
  v[0][1] = v[1][1] = v[4][1] = v[5][1] = -size / 2;
  v[2][1] = v[3][1] = v[6][1] = v[7][1] = size / 2;
  v[0][2] = v[3][2] = v[4][2] = v[7][2] = -size / 2;
  v[1][2] = v[2][2] = v[5][2] = v[6][2] = size / 2;

  for (i = 5; i >= 0; i--) {
    glBegin(type);
    glNormal3fv(&n[i][0]);
    glVertex3fv(&v[faces[i][0]][0]);
    glVertex3fv(&v[faces[i][1]][0]);
    glVertex3fv(&v[faces[i][2]][0]);
    glVertex3fv(&v[faces[i][3]][0]);
    glEnd();
  }
}

void Viewer::DrawAxes() {
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

void Viewer::DrawANN(ANN ann) {
  // draw soma
  for (auto& soma : ann.somas) {
    int s_x = soma.position[0];
    int s_y = soma.position[1];
    int s_z = soma.position[2];
    glColor4f(ann.morphogens[s_x][s_y][s_z][0], ann.morphogens[s_x][s_y][s_z][1], ann.morphogens[s_x][s_y][s_z][2],
              0.4);
    glPushMatrix ();
    glTranslatef(ce(s_x,0), ce(s_y,1), ce(s_z,2));
    DrawBox(0.08, GL_LINE_LOOP);
    glPopMatrix ();
  }

  // draw axons
  for (auto& soma : ann.somas) {
    for (auto& axon : soma.axons) {
      vector<double> slope;
      for (unsigned int d=0; d<Config::N_D; d++) {
        slope.push_back((axon.position[d]-soma.position[d])/11.0);
      }
      for (int i=0; i<11; i++) {
        double x = soma.position[0]+i*slope[0];
        double y = soma.position[1]+i*slope[1];
        double z = soma.position[2]+i*slope[2];
        auto morphs = ann.morphogens[std::round(x)][std::round(y)][std::round(z)];
        glBegin(GL_LINES);
        glColor4f(1.0, 1.0, 1.0, 1.0);
        glVertex3f(ce(x,0), ce(y, 1), ce(z, 2));
        glVertex3f(ce(soma.position[0]+(i+1)*slope[0],0), ce(soma.position[1]+(i+1)*slope[1],1),
                   ce(soma.position[2]+(i+1)*slope[2],2));
        glEnd();
      }
    }
  }
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

void Viewer::run(void (*draw)()) {

  do {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
		KeyInput();
    draw();
    glFlush();
    glfwSwapBuffers(window);
    glfwPollEvents();
  }
	while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
         glfwWindowShouldClose(window) == 0 );

  glfwTerminate();
}
