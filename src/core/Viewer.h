#ifndef VIEWER_HPP
#define VIEWER_HPP
#include <vector>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include "config.hpp"
#include "ANN.h"

class Viewer {
private:
  GLFWwindow* window;
  int rotate_y = 0;
  int rotate_x = 0;
  bool play = true;
  std::vector<int> lengths = {Config::X_SIZE, Config::Y_SIZE, Config::Z_SIZE};

  void KeyInput();

 public:
  double ce(int position, int dim);
  void DrawBox(GLfloat size, GLenum type);
  void DrawAxes();
  void WireBox(double x, double y, double z, double radius);
  void SampleDisplay();
  void DrawANN(ANN ann);

  Viewer();
  void run(void (*draw)());

};
#endif
