#include <random>
#include "external/cxxopts/src/cxxopts.hpp"
#include "core/config.hpp"
#include "problems/forage.hpp"
#include "core/Evaluator.h"
#include "core/types.hpp"
#include "core/Viewer.h"

using std::vector;
using std::cout;
using std::endl;
using std::map;


void step() {
  std::cout << "evaluating " << std::endl;
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
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glVertex3f(-0.6, -0.6, -0.6);
  glVertex3f(-0.6, -0.6, 0.6);
  glEnd();
}

int main(int argc, char** argv) {
  using dna_t = Types::DNAType;
  using problem_t = Types::ProblemType;
  using controller_t = Types::ControllerType;
  Types t;

  dna_t dna = t.random_dna();

  Evaluator<problem_t, controller_t> e(dna, 0);

  Viewer v;
  v.run(display);
}
