#ifndef VIEWER_HPP
#define VIEWER_HPP
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

class Viewer {
private:
  GLFWwindow* window;
  int rotate_y = 0;
  int rotate_x = 0;
  bool play = true;

  void SampleDisplay();
  void KeyInput();

public:

  Viewer();
  void run(void (*draw)());

};
#endif
