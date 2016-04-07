#ifndef ROBOT_H
#define ROBOT_H

#include <vector>
#include <cmath>
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

  Robot() {}

  Robot(double s, double q, double t, double b, int ml, int si) {
    x = s;
    y = q;
    theta = t;
    size = b;
    max_life = ml;
    life = max_life;
    sight = si;
  }

  void move(double left, double right) {
    double phi = theta + M_PI/2.0;
    if (left == right) {
      x = x + right * size * std::cos(phi);
      y = y + right * size * std::sin(phi);
    } else {
      x = x + size * (right + left)/(2.0 * (right - left)) *
        (std::sin((right-left)/size + phi) - std::sin(phi));
      y = y - size * (right + left)/(2.0 * (right - left)) *
        (std::cos((right-left)/size + phi) - std::cos(phi));
    }

    theta = theta + (left - right) / size;
  }

  friend std::ostream& operator<<(std::ostream& out, const Robot& r) {
    return out << "Robot at (" << r.x << "," << r.y << ";" << r.theta << ") life: " << r.life;
  }
};
#endif
