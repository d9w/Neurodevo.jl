#include <stdio.h>
#include <iostream>
#include <stdarg.h>
#include <math.h>
#include <chrono>
#include <thread>
#define GL_GLEXT_PROTOTYPES
#include <GL/glut.h>
#include <GL/glu.h>
#include "easylogging++.h"
#include "color.h"
#include "Constants.h"
#include "Environment.h"
#include "Spikes.h"

INITIALIZE_EASYLOGGINGPP

void display();
void specialKeys();

double rotate_y=0;
double rotate_x=0;
vector<int> lengths = {X_SIZE, Y_SIZE, Z_SIZE};
vector<vector<double> > inputs;

Environment env;
GRN soma_grn;
GRN axon_grn;
Spikes spikes;
int t_action = 0;
int seed = 0;
bool play=true;
int dev_count = 0;
vector<double> max_morphs = {0.01, 0.01, 0.01};
int num_axons = 0;
double fit = 0.0;
int window_1, window_2;

double ce(int position, int dim) {
  double scale;
  if (dim >= 3) {
    scale = 50.0;
  } else {
    scale = (double)lengths[dim];
  }
  return position/scale-0.5;
}

void draw_morphogens() {
   for (auto& soma : env.somas) {
    int s_x = soma.position[0];
    int s_y = soma.position[1];
    int s_z = soma.position[2];
    //glColor4f(env.morphogens[s_x][s_y][s_z][0]/max_morphs[0], env.morphogens[s_x][s_y][s_z][1]/max_morphs[1], env.morphogens[s_x][s_y][s_z][2]/max_morphs[2], 0.2);
    glColor4f(1.0, 1.0, 1.0, 0.4);
    glPushMatrix ();
    glTranslatef(ce(s_x,0), ce(s_y,1), ce(s_z,2));
    glutWireCube(0.09);
    glPopMatrix ();
  }
}

void draw_axes() {
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

void draw_axons() {
  for (auto& soma : env.somas) {
    for (auto& axon : soma.axons) {
      int s_x = axon.position[0];
      int s_y = axon.position[1];
      int s_z = axon.position[2];
      glColor4f(env.morphogens[s_x][s_y][s_z][0]/max_morphs[0], env.morphogens[s_x][s_y][s_z][1]/max_morphs[1], env.morphogens[s_x][s_y][s_z][2]/max_morphs[2], 0.5);
      glPushMatrix ();
      glTranslatef(ce(s_x,0), ce(s_y,1), ce(s_z,2));
      glutSolidCube(0.009);
      glPopMatrix ();
      vector<double> slope;
      for (int d=0; d<N_D; d++) {
        slope.push_back((axon.position[d]-soma.position[d])/11.0);
      }
      for (int i=0; i<11; i++) {
        double x = soma.position[0]+i*slope[0];
        double y = soma.position[1]+i*slope[1];
        double z = soma.position[2]+i*slope[2];
        auto morphs = env.morphogens[std::round(x)][std::round(y)][std::round(z)];
        glBegin(GL_LINES);
        glColor4f(morphs[0]/max_morphs[0], morphs[1]/max_morphs[1], morphs[2]/max_morphs[2], 0.5);
        glVertex3f(ce(x,0), ce(y, 1), ce(z, 2));
        glVertex3f(ce(soma.position[0]+(i+1)*slope[0],0), ce(soma.position[1]+(i+1)*slope[1],1), ce(soma.position[2]+(i+1)*slope[2],2));
        glEnd();
        //LOG(INFO) << "Going from (" << soma.position[0] << " " << soma.position[1] << " " << soma.position[2] << ") and ending up (" << x << " " << y << " " << z << " on the way to (" << axon.position[0] << " " << axon.position[1] << " " << axon.position[2] << ")";
      }
      /*
      glBegin(GL_LINES);
      glColor4f(1.0, 1.0, 1.0, 0.2);
      glVertex3f(ce(soma.position[0],0), ce(soma.position[1],1), ce(soma.position[2],2));
      glVertex3f(ce(axon.position[0],0), ce(axon.position[1],1), ce(axon.position[2],2));
      glEnd();
      */
    }
  }
}

void draw_fired_somas() {
  for (auto& soma : env.somas) {
    if (soma.fired) {
      int s_x = soma.position[0];
      int s_y = soma.position[1];
      int s_z = soma.position[2];
      //double f_t = soma.grn.proteins[SOMA_GRN_OUTPUT_FIRING_THRESH].concentration;
      glColor4f(1.0, 1.0, 1.0, 0.8);
      glPushMatrix ();
      glTranslatef(ce(s_x,0), ce(s_y,1), ce(s_z,2));
      glutSolidCube(0.045);
      glPopMatrix ();
      for (auto& axon : soma.axons) {
        auto rec_soma = env.soma_at(axon.position);
        if (rec_soma != NULL) {
          glBegin(GL_LINES);
          glColor4f(1.0, 1.0, 1.0, 0.8);
          glVertex3f(ce(soma.position[0],0), ce(soma.position[1],1), ce(soma.position[2],2));
          glVertex3f(ce(axon.position[0],0), ce(axon.position[1],1), ce(axon.position[2],2));
          glEnd();
        }
      }
    }
  }
}

void draw_spikes() {
  auto r = spikes.robot;

  glColor4f(0.0, 0.0, 1.0, 1.0);
  for (auto food : spikes.food) {
    glPushMatrix();
    glTranslatef(ce(food[0],3), ce(food[1],3), ce(0.0,3));
    glutSolidCube(1.0/50.0);
    glPopMatrix();
  }

  for (auto food : spikes.nearby_food()) {
    glBegin(GL_LINES);
    glColor4f(0.0, 1.0, 0.0, 1.0);
    glVertex3f(ce(r.x,3),ce(r.y,3),ce(0.0,3));
    glVertex3f(ce(r.x+food[0]*std::cos(r.theta+food[1]),3),
               ce(r.y+food[0]*std::sin(r.theta+food[1]),3),
               ce(0.0,3));
    glEnd();
  }

  glColor4f(1.0, 1.0, 1.0, 1.0);
  glPushMatrix();
  glTranslatef(ce(r.x,3), ce(r.y,3), ce(0.0,3));
  glutSolidCube(r.size/50.0);
  glPopMatrix();

  // draw the robot sight
  for (int i=0; i<=8; i++) {
    double coeff = (double)i/8.0 * M_PI;
    glBegin(GL_LINES);
    glColor4f(1.0, 0.0, 0.0, 1.0);
    glVertex3f(ce(r.x,3),ce(r.y,3),ce(0.0,3));
    glVertex3f(ce(r.x+r.sight*std::cos(r.theta+coeff),3),
               ce(r.y+r.sight*std::sin(r.theta+coeff),3),
               ce(0.0,3));
    glEnd();
  }
}

void draw_counter() {
  glMatrixMode( GL_PROJECTION ) ;
  glPushMatrix() ; // save
  glLoadIdentity();// and clear
  glMatrixMode( GL_MODELVIEW ) ;
  glPushMatrix() ;
  glLoadIdentity() ;

  glDisable( GL_DEPTH_TEST ) ; // also disable the depth test so renders on top

  glRasterPos2f( -0.8,0.8 ) ; // center of screen. (-1,0) is center left.
  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  char buf[300];
  sprintf( buf, "%d", dev_count) ;
  const char * p = buf ;
  do glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, *p ); while( *(++p) ) ;

  glEnable( GL_DEPTH_TEST ) ; // Turn depth testing back on

  glMatrixMode( GL_PROJECTION ) ;
  glPopMatrix() ; // revert back to the matrix I had before.
  glMatrixMode( GL_MODELVIEW ) ;
  glPopMatrix() ;
}

void display(){
  //LOG(DEBUG) << "displaying";
  //  Clear screen and Z-buffer
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  // Reset transformations
  glLoadIdentity();

  glRotatef( rotate_x, 1.0, 0.0, 0.0 );
  glRotatef( rotate_y, 0.0, 1.0, 0.0 );

  draw_morphogens();

  draw_axes();

  //draw_axons();

  draw_fired_somas();

  draw_counter();

  // Rotate when user changes rotate_x and rotate_y
  glFlush();
  glutSwapBuffers();
}

void display_spikes(){
  //  Clear screen and Z-buffer
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();

  draw_axes();

  draw_spikes();

  draw_counter();

  glFlush();
  glutSwapBuffers();
}

void envDev(int iter) {
  int step = spikes.tsteps;
  dev_count = step;
  inputs.clear();
  inputs.push_back({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
  for (int i=1; i<X_SIZE; i++) {
    inputs.push_back({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  }
  auto nearby_food = spikes.nearby_food();
  double reward = 0.0;
  for (auto food : nearby_food) {
    double dist = 1.0 - (food[0] / (double)spikes.robot.sight);
    double angle = food[1];
    int sensor = std::floor(angle/(M_PI/8.0));
    std::cout << " dist " << dist << " angle " << angle << " sensor " << sensor << std::endl;

    reward = std::max(reward, 1.0-dist);
    for (int i=1; i<X_SIZE; i++) {
      if ((double)i/(X_SIZE-1) < dist) {
        inputs[i][sensor] = 1.0;
      }
    }
  }

  /*
  for (auto& input : inputs) {
    for (auto i : input) {
      std::cout << i << " ";
    }
    std::cout << std::endl;
  }
  */

  env.develop_grns(reward);

  env.set_nt_concentration(inputs);

  env.fire_ann();

  std::cout << " OUTPUT " << std::endl;

  double right = 0;
  double left = 0;
  for (auto& soma : env.somas) {
    if (soma.position[2] == lengths[2]-1) {
      if (soma.fired) {
        if (soma.position[0] < 5) {
          left += 1.0;
        }
        if (soma.position[0] > 2) {
          right += 1.0;
        }
      }
    }
  }

  for (int i=0; i<Y_SIZE; i++) {
    for (int j=0; j<X_SIZE; j++) {
      auto soma = env.soma_at({j, i, Z_SIZE-1});
      if (soma->fired) {
        std::cout << "1 ";
      } else {
        std::cout << "0 ";
      }
    }
    std::cout << std::endl;
  }
  std::cout << "right " << right << " left " << left << std::endl;

  right /= (8.0*20.0);
  left /= (8.0*20.0);

  spikes.robot.move(right, left);
  spikes.step();

  glutSetWindow(window_1);
  glutPostRedisplay();

  glutSetWindow(window_2);
  glutPostRedisplay();

  LOG(INFO) << spikes.robot;

  if (spikes.robot.life < 0) {
    seed += 1;
    if (seed < 5) {
      env = Environment(lengths, soma_grn, axon_grn);
      env.set_random_connectivity(seed);
      spikes = Spikes(50, 50, seed);
    }
  }
  if (play && seed < 5) {
    glutTimerFunc(25, envDev, 0);
  }
}

void specialKeys( int key, int x, int y ) {
  if (key == GLUT_KEY_RIGHT)
    rotate_y += 5;
  else if (key == GLUT_KEY_LEFT)
    rotate_y -= 5;
  else if (key == GLUT_KEY_UP)
    rotate_x += 5;
  else if (key == GLUT_KEY_DOWN)
    rotate_x -= 5;
  else if (key == GLUT_KEY_F1) {
    if (play) {
      play = false;
    } else {
      play = true;
      glutTimerFunc(25, envDev, 0);
    }
  }
  glutSetWindow(window_1);
  glutPostRedisplay();
}

int main(int argc, char* argv[]){
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  el::Configurations conf("log.conf");
  el::Loggers::reconfigureAllLoggers(conf);

  srand(time(NULL));
  soma_grn = GRN(argv[1]);
  axon_grn = GRN(argv[2]);
  std::cout << "Soma grn: " << std::endl;
  std::cout << soma_grn.toString();
  std::cout << std::endl << "Axon grn: " << std::endl;
  std::cout << axon_grn.toString();
  env = Environment(lengths, soma_grn, axon_grn);
  env.set_random_connectivity(seed);
  spikes = Spikes(50, 50, seed);
  t_action = (int)std::round((soma_grn.t_action+axon_grn.t_action)/2.0);

  //glutPostRedisplay();

  window_1 = glutCreateWindow("brainz");
  glutSetWindow(window_1);
  glutReshapeWindow(960, 960);
  //glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH)-960)/4,540);
  glutPositionWindow(0, 0);

  //  Enable Z-buffer depth test
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // Callback functions
  glutDisplayFunc(display);
  glutSpecialFunc(specialKeys);
  glutTimerFunc(25, envDev, 0);

  window_2 = glutCreateWindow("brainz");
  glutSetWindow(window_2);
  glutReshapeWindow(960, 960);
  //glutInitWindowPosition(3*(glutGet(GLUT_SCREEN_WIDTH)-960)/4, 540);
  glutPositionWindow(960, 0);

  glutDisplayFunc(display_spikes);

  //  Pass control to GLUT for events
  glutMainLoop();
  //glutPostRedisplay();
  return 0;
}
