#ifndef MNIST_H
#define MNIST_H

#include <vector>
#include <string>

class MNIST {
 public:
  std::vector<std::vector<double> > images;
  std::vector<int> labels;
  std::vector<std::vector<double> > centers;
  int number_of_images;
  int n_rows;
  int n_cols;

  MNIST();

  void normalize_images();
  void get_images(std::string filename);
  void get_labels(std::string filename);
  void bit_images(int bits);
  void remove_zeros();
  void random_circles(int x, int y, int imgs);
  void irises();
};

#endif
