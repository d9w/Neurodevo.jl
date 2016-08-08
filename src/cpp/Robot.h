#ifndef ROBOT_H
#define ROBOT_H

#include <ostream>

class Robot {
 public:
  double x;
  double y;
  double theta;
  double size;
  int life;
  int max_life;
  int sight;

  Robot();
  Robot(double s, double q, double t, double b, int ml, int si);

  void move(double left, double right);

  friend std::ostream& operator<<(std::ostream& out, const Robot& r);
};

#endif
