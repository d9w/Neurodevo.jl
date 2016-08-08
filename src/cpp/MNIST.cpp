#include "MNIST.h"
#include <fstream>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include "easylogging++.h"

int reverse_int (int i) {
  unsigned char c1, c2, c3, c4;

  c1 = i & 255;
  c2 = (i >> 8) & 255;
  c3 = (i >> 16) & 255;
  c4 = (i >> 24) & 255;

  return ((int)c1 << 24) + ((int)c2 << 16) + ((int)c3 << 8) + c4;
}

MNIST::MNIST() {};

void MNIST::normalize_images() {
  int imsize = images[0].size();
  std::vector<double> maxs (imsize, 0.0);

  for (auto& image : images) {
    for (int i=0; i<imsize; i++) {
      if (image[i] > maxs[i]) {
        maxs[i] = image[i];
      }
    }
  }

  for (auto& image : images) {
    for (int i=0; i<imsize; i++) {
      if (maxs[i] > 0.0) {
        image[i] = image[i]/maxs[i];
      }
    }
  }
}

void MNIST::get_images(std::string filename) {
  int magic_number=0;
  number_of_images=0;
  n_rows=0;
  n_cols=0;

  std::ifstream file (filename, std::ios::binary);
  if (file.is_open()) {
    file.read((char*)&magic_number,sizeof(magic_number));
    magic_number= reverse_int(magic_number);
    file.read((char*)&number_of_images,sizeof(number_of_images));
    number_of_images= reverse_int(number_of_images);
    file.read((char*)&n_rows,sizeof(n_rows));
    n_rows= reverse_int(n_rows);
    file.read((char*)&n_cols,sizeof(n_cols));
    n_cols= reverse_int(n_cols);

    images.resize(number_of_images, std::vector<double>(n_rows*n_cols));

    for(int i=0; i<number_of_images; i++) {
      for(int r=0; r<n_rows; r++) {
        for(int c=0; c<n_cols; c++) {
          unsigned char temp=0;
          file.read((char*)&temp,sizeof(temp));
          images[i][(n_rows*r)+c]= (double)temp;
        }
      }
    }
  }
}

void MNIST::get_labels(std::string filename) {
  int magic_number=0;
  number_of_images=0;

  std::ifstream file (filename, std::ios::binary);
  if (file.is_open()) {
    file.read((char*)&magic_number,sizeof(magic_number));
    magic_number= reverse_int(magic_number);
    file.read((char*)&number_of_images,sizeof(number_of_images));
    number_of_images= reverse_int(number_of_images);

    labels.resize(number_of_images);

    for(int i=0;i<number_of_images;++i) {
      unsigned char temp=0;
      file.read((char*)&temp,sizeof(temp));
      labels[i] = (int)temp;
    }
  }
}

void MNIST::bit_images(int bits) {
  for (int j=0; j<pow(2,bits); j++) {
    std::vector<double> image;
    int num=j;
    int r;
    int label=0;
    for (int z=0; z<bits; z++) {
      r = num%2;
      image.insert(image.begin(), r);
      num /= 2;
      label += r;
    }
    label = label%2;
    labels.push_back(label);
    images.push_back(image);
  }
}

void MNIST::remove_zeros() {
  for (int i=0; i<labels.size(); i++) {
    if (labels[i]==0) {
      labels.erase(labels.begin()+i);
      images.erase(images.begin()+i);
      i = i-1;
    }
  }
}

bool overlaps(double ic_x, double ic_y, double ir, double jc_x, double jc_y, double jr) {
  return (ic_x - jc_x)*(ic_x - jc_x) + (ic_y - jc_y)*(ic_y - jc_y) <= (ir + jr)*(ir + jr);
}

void MNIST::random_circles(int x, int y, int imgs) {
  for (int i=0; i<imgs; i++) {
    double ic_x, ic_y, ir;
    double jc_x, jc_y, jr;
    do {
      ic_x = (double)rand()/RAND_MAX * x;
      ic_y = (double)rand()/RAND_MAX * y;
      ir = 1.5 + 0.5 * (double)rand()/RAND_MAX;
      jc_x = 1.0 * (rand() % (x+1));
      jc_y = 1.0 * (rand() % (y+1));
      jr = 1;
    } while (overlaps(ic_x, ic_y, ir, jc_x, jc_y, jr));
    std::vector<double> image;
    for (double ix=0; ix<x; ix++) {
      for (double iy=0; iy<y; iy++) {
        if (overlaps(ic_x, ic_y, ir, ix, iy, 1.0) || overlaps(jc_x, jc_y, jr, ix, iy, 1.0)) {
          image.push_back(1.0);
        } else {
          image.push_back(0.0);
        }
      }
    }
    LOG(DEBUG) << i << " map: " << ic_x << " " << ic_y << " " << ir << " | " << jc_x << " " << jc_y << " " << jr;
    images.push_back(image);
    centers.push_back({ic_x, ic_y});
  }
}

void MNIST::irises() {
  images.push_back({5.1,3.5,1.4,0.2});
  labels.push_back(0);
  images.push_back({4.9,3.0,1.4,0.2});
  labels.push_back(0);
  images.push_back({4.7,3.2,1.3,0.2});
  labels.push_back(0);
  images.push_back({4.6,3.1,1.5,0.2});
  labels.push_back(0);
  images.push_back({5.0,3.6,1.4,0.2});
  labels.push_back(0);
  images.push_back({5.4,3.9,1.7,0.4});
  labels.push_back(0);
  images.push_back({4.6,3.4,1.4,0.3});
  labels.push_back(0);
  images.push_back({5.0,3.4,1.5,0.2});
  labels.push_back(0);
  images.push_back({4.4,2.9,1.4,0.2});
  labels.push_back(0);
  images.push_back({4.9,3.1,1.5,0.1});
  labels.push_back(0);
  images.push_back({5.4,3.7,1.5,0.2});
  labels.push_back(0);
  images.push_back({4.8,3.4,1.6,0.2});
  labels.push_back(0);
  images.push_back({4.8,3.0,1.4,0.1});
  labels.push_back(0);
  images.push_back({4.3,3.0,1.1,0.1});
  labels.push_back(0);
  images.push_back({5.8,4.0,1.2,0.2});
  labels.push_back(0);
  images.push_back({5.7,4.4,1.5,0.4});
  labels.push_back(0);
  images.push_back({5.4,3.9,1.3,0.4});
  labels.push_back(0);
  images.push_back({5.1,3.5,1.4,0.3});
  labels.push_back(0);
  images.push_back({5.7,3.8,1.7,0.3});
  labels.push_back(0);
  images.push_back({5.1,3.8,1.5,0.3});
  labels.push_back(0);
  images.push_back({5.4,3.4,1.7,0.2});
  labels.push_back(0);
  images.push_back({5.1,3.7,1.5,0.4});
  labels.push_back(0);
  images.push_back({4.6,3.6,1.0,0.2});
  labels.push_back(0);
  images.push_back({5.1,3.3,1.7,0.5});
  labels.push_back(0);
  images.push_back({4.8,3.4,1.9,0.2});
  labels.push_back(0);
  images.push_back({5.0,3.0,1.6,0.2});
  labels.push_back(0);
  images.push_back({5.0,3.4,1.6,0.4});
  labels.push_back(0);
  images.push_back({5.2,3.5,1.5,0.2});
  labels.push_back(0);
  images.push_back({5.2,3.4,1.4,0.2});
  labels.push_back(0);
  images.push_back({4.7,3.2,1.6,0.2});
  labels.push_back(0);
  images.push_back({4.8,3.1,1.6,0.2});
  labels.push_back(0);
  images.push_back({5.4,3.4,1.5,0.4});
  labels.push_back(0);
  images.push_back({5.2,4.1,1.5,0.1});
  labels.push_back(0);
  images.push_back({5.5,4.2,1.4,0.2});
  labels.push_back(0);
  images.push_back({4.9,3.1,1.5,0.1});
  labels.push_back(0);
  images.push_back({5.0,3.2,1.2,0.2});
  labels.push_back(0);
  images.push_back({5.5,3.5,1.3,0.2});
  labels.push_back(0);
  images.push_back({4.9,3.1,1.5,0.1});
  labels.push_back(0);
  images.push_back({4.4,3.0,1.3,0.2});
  labels.push_back(0);
  images.push_back({5.1,3.4,1.5,0.2});
  labels.push_back(0);
  images.push_back({5.0,3.5,1.3,0.3});
  labels.push_back(0);
  images.push_back({4.5,2.3,1.3,0.3});
  labels.push_back(0);
  images.push_back({4.4,3.2,1.3,0.2});
  labels.push_back(0);
  images.push_back({5.0,3.5,1.6,0.6});
  labels.push_back(0);
  images.push_back({5.1,3.8,1.9,0.4});
  labels.push_back(0);
  images.push_back({4.8,3.0,1.4,0.3});
  labels.push_back(0);
  images.push_back({5.1,3.8,1.6,0.2});
  labels.push_back(0);
  images.push_back({4.6,3.2,1.4,0.2});
  labels.push_back(0);
  images.push_back({5.3,3.7,1.5,0.2});
  labels.push_back(0);
  images.push_back({5.0,3.3,1.4,0.2});
  labels.push_back(0);
  images.push_back({7.0,3.2,4.7,1.4});
  labels.push_back(1);
  images.push_back({6.4,3.2,4.5,1.5});
  labels.push_back(1);
  images.push_back({6.9,3.1,4.9,1.5});
  labels.push_back(1);
  images.push_back({5.5,2.3,4.0,1.3});
  labels.push_back(1);
  images.push_back({6.5,2.8,4.6,1.5});
  labels.push_back(1);
  images.push_back({5.7,2.8,4.5,1.3});
  labels.push_back(1);
  images.push_back({6.3,3.3,4.7,1.6});
  labels.push_back(1);
  images.push_back({4.9,2.4,3.3,1.0});
  labels.push_back(1);
  images.push_back({6.6,2.9,4.6,1.3});
  labels.push_back(1);
  images.push_back({5.2,2.7,3.9,1.4});
  labels.push_back(1);
  images.push_back({5.0,2.0,3.5,1.0});
  labels.push_back(1);
  images.push_back({5.9,3.0,4.2,1.5});
  labels.push_back(1);
  images.push_back({6.0,2.2,4.0,1.0});
  labels.push_back(1);
  images.push_back({6.1,2.9,4.7,1.4});
  labels.push_back(1);
  images.push_back({5.6,2.9,3.6,1.3});
  labels.push_back(1);
  images.push_back({6.7,3.1,4.4,1.4});
  labels.push_back(1);
  images.push_back({5.6,3.0,4.5,1.5});
  labels.push_back(1);
  images.push_back({5.8,2.7,4.1,1.0});
  labels.push_back(1);
  images.push_back({6.2,2.2,4.5,1.5});
  labels.push_back(1);
  images.push_back({5.6,2.5,3.9,1.1});
  labels.push_back(1);
  images.push_back({5.9,3.2,4.8,1.8});
  labels.push_back(1);
  images.push_back({6.1,2.8,4.0,1.3});
  labels.push_back(1);
  images.push_back({6.3,2.5,4.9,1.5});
  labels.push_back(1);
  images.push_back({6.1,2.8,4.7,1.2});
  labels.push_back(1);
  images.push_back({6.4,2.9,4.3,1.3});
  labels.push_back(1);
  images.push_back({6.6,3.0,4.4,1.4});
  labels.push_back(1);
  images.push_back({6.8,2.8,4.8,1.4});
  labels.push_back(1);
  images.push_back({6.7,3.0,5.0,1.7});
  labels.push_back(1);
  images.push_back({6.0,2.9,4.5,1.5});
  labels.push_back(1);
  images.push_back({5.7,2.6,3.5,1.0});
  labels.push_back(1);
  images.push_back({5.5,2.4,3.8,1.1});
  labels.push_back(1);
  images.push_back({5.5,2.4,3.7,1.0});
  labels.push_back(1);
  images.push_back({5.8,2.7,3.9,1.2});
  labels.push_back(1);
  images.push_back({6.0,2.7,5.1,1.6});
  labels.push_back(1);
  images.push_back({5.4,3.0,4.5,1.5});
  labels.push_back(1);
  images.push_back({6.0,3.4,4.5,1.6});
  labels.push_back(1);
  images.push_back({6.7,3.1,4.7,1.5});
  labels.push_back(1);
  images.push_back({6.3,2.3,4.4,1.3});
  labels.push_back(1);
  images.push_back({5.6,3.0,4.1,1.3});
  labels.push_back(1);
  images.push_back({5.5,2.5,4.0,1.3});
  labels.push_back(1);
  images.push_back({5.5,2.6,4.4,1.2});
  labels.push_back(1);
  images.push_back({6.1,3.0,4.6,1.4});
  labels.push_back(1);
  images.push_back({5.8,2.6,4.0,1.2});
  labels.push_back(1);
  images.push_back({5.0,2.3,3.3,1.0});
  labels.push_back(1);
  images.push_back({5.6,2.7,4.2,1.3});
  labels.push_back(1);
  images.push_back({5.7,3.0,4.2,1.2});
  labels.push_back(1);
  images.push_back({5.7,2.9,4.2,1.3});
  labels.push_back(1);
  images.push_back({6.2,2.9,4.3,1.3});
  labels.push_back(1);
  images.push_back({5.1,2.5,3.0,1.1});
  labels.push_back(1);
  images.push_back({5.7,2.8,4.1,1.3});
  labels.push_back(1);
  images.push_back({6.3,3.3,6.0,2.5});
  labels.push_back(2);
  images.push_back({5.8,2.7,5.1,1.9});
  labels.push_back(2);
  images.push_back({7.1,3.0,5.9,2.1});
  labels.push_back(2);
  images.push_back({6.3,2.9,5.6,1.8});
  labels.push_back(2);
  images.push_back({6.5,3.0,5.8,2.2});
  labels.push_back(2);
  images.push_back({7.6,3.0,6.6,2.1});
  labels.push_back(2);
  images.push_back({4.9,2.5,4.5,1.7});
  labels.push_back(2);
  images.push_back({7.3,2.9,6.3,1.8});
  labels.push_back(2);
  images.push_back({6.7,2.5,5.8,1.8});
  labels.push_back(2);
  images.push_back({7.2,3.6,6.1,2.5});
  labels.push_back(2);
  images.push_back({6.5,3.2,5.1,2.0});
  labels.push_back(2);
  images.push_back({6.4,2.7,5.3,1.9});
  labels.push_back(2);
  images.push_back({6.8,3.0,5.5,2.1});
  labels.push_back(2);
  images.push_back({5.7,2.5,5.0,2.0});
  labels.push_back(2);
  images.push_back({5.8,2.8,5.1,2.4});
  labels.push_back(2);
  images.push_back({6.4,3.2,5.3,2.3});
  labels.push_back(2);
  images.push_back({6.5,3.0,5.5,1.8});
  labels.push_back(2);
  images.push_back({7.7,3.8,6.7,2.2});
  labels.push_back(2);
  images.push_back({7.7,2.6,6.9,2.3});
  labels.push_back(2);
  images.push_back({6.0,2.2,5.0,1.5});
  labels.push_back(2);
  images.push_back({6.9,3.2,5.7,2.3});
  labels.push_back(2);
  images.push_back({5.6,2.8,4.9,2.0});
  labels.push_back(2);
  images.push_back({7.7,2.8,6.7,2.0});
  labels.push_back(2);
  images.push_back({6.3,2.7,4.9,1.8});
  labels.push_back(2);
  images.push_back({6.7,3.3,5.7,2.1});
  labels.push_back(2);
  images.push_back({7.2,3.2,6.0,1.8});
  labels.push_back(2);
  images.push_back({6.2,2.8,4.8,1.8});
  labels.push_back(2);
  images.push_back({6.1,3.0,4.9,1.8});
  labels.push_back(2);
  images.push_back({6.4,2.8,5.6,2.1});
  labels.push_back(2);
  images.push_back({7.2,3.0,5.8,1.6});
  labels.push_back(2);
  images.push_back({7.4,2.8,6.1,1.9});
  labels.push_back(2);
  images.push_back({7.9,3.8,6.4,2.0});
  labels.push_back(2);
  images.push_back({6.4,2.8,5.6,2.2});
  labels.push_back(2);
  images.push_back({6.3,2.8,5.1,1.5});
  labels.push_back(2);
  images.push_back({6.1,2.6,5.6,1.4});
  labels.push_back(2);
  images.push_back({7.7,3.0,6.1,2.3});
  labels.push_back(2);
  images.push_back({6.3,3.4,5.6,2.4});
  labels.push_back(2);
  images.push_back({6.4,3.1,5.5,1.8});
  labels.push_back(2);
  images.push_back({6.0,3.0,4.8,1.8});
  labels.push_back(2);
  images.push_back({6.9,3.1,5.4,2.1});
  labels.push_back(2);
  images.push_back({6.7,3.1,5.6,2.4});
  labels.push_back(2);
  images.push_back({6.9,3.1,5.1,2.3});
  labels.push_back(2);
  images.push_back({5.8,2.7,5.1,1.9});
  labels.push_back(2);
  images.push_back({6.8,3.2,5.9,2.3});
  labels.push_back(2);
  images.push_back({6.7,3.3,5.7,2.5});
  labels.push_back(2);
  images.push_back({6.7,3.0,5.2,2.3});
  labels.push_back(2);
  images.push_back({6.3,2.5,5.0,1.9});
  labels.push_back(2);
  images.push_back({6.5,3.0,5.2,2.0});
  labels.push_back(2);
  images.push_back({6.2,3.4,5.4,2.3});
  labels.push_back(2);
  images.push_back({5.9,3.0,5.1,1.8});
  labels.push_back(2);

  for (int i=0; i<images.size(); i++) {
    int index = rand() % (i+1);
    auto temp_img = images[index];
    auto temp_label = labels[index];
    images[index] = images[i];
    labels[index] = labels[i];
    images[i] = temp_img;
    labels[i] = temp_label;
  }
}
