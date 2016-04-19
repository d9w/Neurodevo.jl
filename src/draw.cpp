#include "external/cxxopts/src/cxxopts.hpp"
#include "core/Viewer.h"
#include "core/config.hpp"
#include "core/DNA.h"
#include "core/ANN.h"
#include "problems/Forage.h"

using std::vector;
using std::cout;
using std::endl;
using std::map;

DNA dna;
Forage forage{0};
ANN ann{dna, 0};
vector<vector<double>> outputs;
vector<double> max_morphs = {0.01, 0.01, 0.01};
vector<int> lengths = {Config::X_SIZE, Config::Y_SIZE, Config::Z_SIZE};

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
  for (auto& soma : ann.somas) {
    int s_x = soma.position[0];
    int s_y = soma.position[1];
    int s_z = soma.position[2];
    glColor4f(ann.morphogens[s_x][s_y][s_z][0]/max_morphs[0], ann.morphogens[s_x][s_y][s_z][1]/max_morphs[1],
              ann.morphogens[s_x][s_y][s_z][2]/max_morphs[2], 0.4);
    //glColor4f(1.0, 1.0, 1.0, 0.4);
    glPushMatrix ();
    glTranslatef(ce(s_x,0), ce(s_y,1), ce(s_z,2));
    //glutWireCube(0.09);
    glPopMatrix ();
  }
}

void draw_axons() {
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
        //glColor4f(morphs[0]/max_morphs[0], morphs[1]/max_morphs[1], morphs[2]/max_morphs[2], 0.5);
        glColor4f(1.0, 1.0, 1.0, 1.0);
        glVertex3f(ce(x,0), ce(y, 1), ce(z, 2));
        glVertex3f(ce(soma.position[0]+(i+1)*slope[0],0), ce(soma.position[1]+(i+1)*slope[1],1), ce(soma.position[2]+(i+1)*slope[2],2));
        glEnd();
      }
    }
  }
}

void display() {
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

  draw_axons();
}

void step() {
  if (!forage.stop()) {
    auto ins = forage.getInputs();
    ann.step(ins, forage.getReward());
    ann.set_outputs(&outputs);
    forage.step(outputs);
  }

  display();
}

int main(int, char**) {
  Viewer v;
  v.run(step);
}
