// Project 3 Yu Deng and Shaowei Qiu
#include "CImg.h"
#include "p2.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#define TWO_PI 6.28318530
#define PI 3.14159265359
#define SIGN(a) ((a)>0 ? 1:-1)
using namespace cimg_library;
using namespace std;

typedef std::pair<std::pair<int, int>, double> tuple;

struct clr { 
  int r;
  int g;
  int b;
};

// declaration of functions
// 2
// --2.1
CImg<double> edge_detector(CImg<double> input);
//2,2
CImg<double> circle_detector(CImg<double> input, CImg<double> sx, CImg<double> sy);
CImg<double> edgeusing2ndder(CImg<double> input,CImg<double> sobel);
CImg<double> removesmallervalues(CImg<double>input, int x, int y);
CImg<double> simpleBW(CImg<double> input, int thresh);
CImg<double> filter(CImg<double> input, CImg<double> filter);
CImg<double> Gaussian_Generator(double sigma);
CImg<double> averageGrayscale(CImg<double> input);
CImg<double> gaussianClrFilterGen(double sigma);
void fft(const CImg<double> &input, CImg<double> &fft_real, CImg<double> &fft_imag);
void ifft(const CImg<double> &input_real, const CImg<double> &input_imag, CImg<double> &output_real);
// 3
// --3.1
std::vector<int> prod_segmentation(CImg<double> input);
std::map<std::pair<int, int>, double> graph_generator(CImg<double> input);
std::vector<tuple> sort_graph(std::map<std::pair<int, int>, double> graph);
double calc_Threshold(std::vector<tuple> input);
//--3.2
std::vector<int> Real_Segmentation(CImg<double> input, double k);
int calculate_size(std::vector<int> data, int length, int grp);

CImg<double> simpleBW(CImg<double> input, int thresh){
        int x, y;
        int input_w = input.width();
        int input_h = input.height();
        CImg<double> output(input, "xyzc", 0);

        for(x = 0; x < input_w; x++) {
            for(y = 0; y < input_h; y++) {
              if(input(x, y, 0, 0) > thresh) {
                output(x,y, 0, 0) = 255;
              }
            }
        }
        return output;
}

// 2.1
CImg<double> edge_detector(CImg<double> input) {
  int x,y;
  CImg<double> output(input, "xyzc", 0);
  CImg<double> y_deriv(input);
  CImg<double> x_deriv(input);
  CImg<double> sx(3, 3, 1, 1);
  CImg<double> sy(3, 3, 1, 1);

 sx(0,0,0,0) = -1;  sx(1,0,0,0) = -2;  sx(2,0,0,0) = -1;
  sx(0,1,0,0) = 0;   sx(1,1,0,0) = 0;  sx(2,1,0,0) = 0;
  sx(0,2,0,0) = 1;  sx(1,2,0,0) = 2;  sx(2,2,0,0) = 1;

  sy(0,0,0,0) = 1;   sy(1,0,0,0) = 0;  sy(2,0,0,0) = -1;
  sy(0,1,0,0) = 2;   sy(1,1,0,0) = 0;  sy(2,1,0,0) = -2;
  sy(0,2,0,0) = 1;  sy(1,2,0,0) = 0;  sy(2,2,0,0) = -1;

  x_deriv.convolve(sx);
  y_deriv.convolve(sy);

  x_deriv.normalize(0,255);
  y_deriv.normalize(0,255);
  
  x_deriv.save("sx.jpg");
  y_deriv.save("sy.jpg");

  for(x = 0; x < input.width(); x++) {
    for(y = 0; y < input.height(); y++) {
      output(x,y,0,0) = sqrt( pow(x_deriv(x,y,0,0), 2) + pow(y_deriv(x,y,0,0), 2));
    }
  }
  return simpleBW(output, 200);
}

double calculate_Gaussian(int x, int y, int sigma)
{
    return ((1/(TWO_PI*(sigma*sigma)))*exp(-(((double)(x*x)+(y*y))/(2*(sigma*sigma)))));

}

CImg<double> Gaussian_Generator(double sigma){
  int filterSize = (2*3*sigma)+1;
  CImg<double> H(filterSize, filterSize, 1, 1);
  double value=0,value2=0.0;
  for(int y = 0; y < filterSize; y++) {
    for(int x = 0; x < filterSize; x++)
      {
H(x,y,0,0) = calculate_Gaussian(x-(filterSize/2),y-(filterSize/2),sigma);
value+=H(x,y,0,0);
      }
  }

  for(int y = 0; y < filterSize; y++) {
    for(int x = 0; x < filterSize; x++)
      {
H(x,y,0,0) = H(x,y,0,0)/value;
      }
  }
  return H;
}



CImg<double> gaussianClrFilterGen(double sigma){
        int filterSize = 3*sigma;
        CImg<double> H(filterSize, filterSize, 1, 3);
        double value=0,value2=0.0;
        double valuea=0,value2a=0.0;
        double valueb=0,value2b=0.0;
        for(int y = 0; y < filterSize; y++) {
            for(int x = 0; x < filterSize; x++)
            {
              H(x,y,0,0) = calculate_Gaussian(x-(filterSize/2),y-(filterSize/2),sigma);
              value+=H(x,y,0,0);
              H(x,y,0,2) = calculate_Gaussian(x-(filterSize/2),y-(filterSize/2),sigma);
              valuea+=H(x,y,0,1);
              H(x,y,0,1) = calculate_Gaussian(x-(filterSize/2),y-(filterSize/2),sigma);
              valueb+=H(x,y,0,2);
            }
        }

        for(int y = 0; y < filterSize; y++) {
            for(int x = 0; x < filterSize; x++)
            {
                H(x,y,0,0) = H(x,y,0,0)/value;
                H(x,y,0,1) = H(x,y,0,1)/valuea;
                H(x,y,0,2) = H(x,y,0,2)/valueb;
            }
        }
        return H;

}

//2.2 reference https://files.nyu.edu/jb4457/public/files/research/bristol/hough-report.pdf reference
CImg<double> circle_detector(CImg<double> input, CImg<double> sx, CImg<double> sy) {
//define hough parameters a,b,r
int a,b;
double Min_r=10.0;
double Max_r=60.0;
double dx;
double dy;
double dx1;
double dy1;
double x1,y1;
int x2,y2,x3,y3;
//set threshold value of the vote to determine whether a circle is detected
int threshold=2;
//Create Hough image space
CImg<double> Hough_space(input,"xyzc",0);
CImg<double> edge_direction(input,"xyzc",0);
for(int x=0;x<input.width();x++){
for(int y=0;y<input.height();y++){
// calculating the direction map
if(input(x,y,0,0)==0){
edge_direction(x,y,0,0)=0;
}
else{
edge_direction(x,y,0,0)=atan2(sy(x,y,0,0),sx(x,y,0,0));
if(edge_direction(x,y,0,0) > PI/2 && edge_direction(x,y,0,0)<PI){
edge_direction(x,y,0,0) -= PI;}
if(edge_direction(x,y,0,0) < -PI/2 && edge_direction(x,y,0,0)>-PI){
edge_direction(x,y,0,0) += PI;}
}
}
}
for(int x=10;x<input.width()-10;x++){
for(int y=10;y<input.height()-10;y++){
if(input(x,y,0,0)>0){
x1 = Min_r * cos(edge_direction(x,y,0,0));
y1 = Min_r * sin(edge_direction(x,y,0,0));
if((edge_direction(x,y,0,0) > -PI/4) && (edge_direction(x,y,0,0)< PI/4)){
dx = SIGN(x1);
dy=round(dx*tan(edge_direction(x,y,0,0)));
}
else{
dy = SIGN(y1);
dx=round(dy/tan(edge_direction(x,y,0,0)));
}
while(sqrt(x1*x1+y1*y1)<Max_r){
x2 = (int)(x +x1); y2 = (int)(y-y1);
x3 = (int)(x -x1); y3 = (int)(y+y1);
if((x2<input.width())&&(x2>=0)&&(y2< input.height())&&(y2>=0)){
Hough_space(x2,y2,0,0)+=input(x,y,0,0)/sqrt(x1*x1+y1*y1);}
if((x3<input.width())&&(x3>=0)&&(y3< input.height())&&(y3>=0)){
Hough_space(x3,y3,0,0)+=input(x,y,0,0)/sqrt(x1*x1+y1*y1);}
x1=x1+dx;
y1=y1+dy;
}
}
}
}
Hough_space = filter(Hough_space, Gaussian_Generator(2));//apply the enhanceabspace
/*for (int x = 0; x < input.width(); x++){ //remove individual dots
for (int y = 0; y < input.height(); y++)
{
int maxval = 0;
for (int k = -3; k < 4; k++)
for (int l = -3; l < 4; l++)
if ((x + k >= 0) && (y + l >= 0) && (x + k < input.height()) && (y + l < input.width()) &&
((k != 0) || (l != 0)))
if (Hough_space(x+k,y+l,0,0) > 0) maxval++;
if (maxval < 3) Hough_space(x,y,0,0) = 0;
}
}*/
return Hough_space.normalize(0,255);
}

CImg<double> edgeusing2ndder(CImg<double> input, CImg<double> sobel){
int x,y;
CImg<double> output(input, "xyzc", 0);

CImg<double> F(5, 5, 1, 1);
F(0,0,0,0) = 0; F(1,0,0,0) = 0; F(2,0,0,0) = -1; F(3,0,0,0)=0;  F(4,0,0,0)=0;
F(0,1,0,0) = 0; F(1,1,0,0) = -1; F(2,1,0,0) = -2; F(3,1,0,0)=-1; F(4,1,0,0)=0;
F(0,2,0,0) = -1; F(1,2,0,0) = -2; F(2,2,0,0) = 16; F(3,2,0,0)=-2;F(4,2,0,0)=-1;
F(0,3,0,0) = 0; F(1,3,0,0) = -1; F(2,3,0,0) = -2;  F(3,3,0,0)=-1;F(4,3,0,0)=0;
F(0,4,0,0) = 0; F(1,4,0,0) = 0; F(2,4,0,0) = -1;  F(3,4,0,0)=0;F(4,4,0,0)=0;

input.convolve(F);
for(int x=0;x<input.width();x++){
for(int y=0;y<input.height();y++){
if((x>1)&&(y>1)&&(x<input.width()-2)&&(y<input.height()-2)){
if(abs(input(x,y,0,0) +input(x+1,y,0,0) + input(x,y+1,0,0)) == abs(input(x,y,0,0))+ abs(input(x+1,y,0,0))+ abs(input(x,y+1,0,0))){
sobel(x,y,0,0) = 0;
}
else
{
  sobel(x,y,0,0)= sobel(x,y,0,0);
  }
  }
  }
  }
  return sobel.normalize(0,255);
  }
  
  
CImg<double> averageGrayscale(CImg<double> input){
  int x, y, col;
  double sum;
  int input_w = input.width();
  int input_h = input.height();
  CImg<double> output(input_w, input_h, 1, 1);
  for(x = 0; x < input_w; x++) {
    for(y = 0; y < input_h; y++) {
      sum = 0.0;
      sum += input(x, y, 0, 0) * 0.21;
      sum += input(x, y, 0, 1) * 0.71;
      sum += input(x, y, 0, 2) * 0.07;
      output(x,y,0,0) = sum;
    }
  }
  return output;
  }

double calc_mvng_avg(double val, double n, double prev_avg) {
  return ((val + n * prev_avg) / (n + 1.0));
}
//3.1
CImg<double> segmentation(CImg<double> input, int threshold){

  CImg<double> graph(input.width(),input.height(),1,4);
  // int number_block = input.width() * input.height();
  //0-left 1-right 2-up 3-down
  //CImg<double> connect_set(number_block,number_block,1,1);
  double MAX_WEIGHT = -1;
  for(int x= 0; x< input.width(); x++){
    for(int y= 0; y< input.height(); y++){
      double r_diff, g_diff, b_diff, color_diff;
      if(x-1 >= 0){
  r_diff = pow(input(x,y,0,0)-input(x-1,y,0,0),2);
        g_diff = pow(input(x,y,0,1)-input(x-1,y,0,1,0),2);
        b_diff = pow(input(x,y,0,2)-input(x-1,y,0,2,0),2);
        color_diff = sqrt(r_diff + g_diff + b_diff);
        if(color_diff < threshold) graph(x,y,0,0) = 0;
        else { graph(x,y,0,0) = color_diff;}
      }
      else{ graph(x,y,0,0) = -1;}
      if(x + 1 <= input.width()){
        r_diff = pow(input(x,y,0,0)-input(x+1,y,0,0),2);
        g_diff = pow(input(x,y,0,1)-input(x+1,y,0,1),2);
        b_diff = pow(input(x,y,0,2)-input(x+1,y,0,2),2);
        color_diff = sqrt(r_diff + g_diff + b_diff);
  if(color_diff < threshold) graph(x,y,0,1) = 0;
        else{graph(x,y,0,1) = color_diff;}
      }
      else{ graph(x,y,0,1) = -1;}
      if(y - 1 >= 0){
        r_diff = pow(input(x,y,0,0)-input(x,y-1,0,0),2);
        g_diff = pow(input(x,y,0,1)-input(x,y-1,0,1),2);
        b_diff = pow(input(x,y,0,2)-input(x,y-1,0,2),2);
        color_diff = sqrt(r_diff + g_diff + b_diff);
        if(color_diff < threshold) graph(x,y,0,2) = 0;
        else{ graph(x,y,0,2) = color_diff;}
      }
      else{ graph(x,y,0,2) = -1;}
      if(y + 1 <= input.height()){
        r_diff = pow(input(x,y,0,0)-input(x-1,y,0,0),2);
        g_diff = pow(input(x,y,0,1)-input(x-1,y,0,1),2);
        b_diff = pow(input(x,y,0,2)-input(x-1,y,0,2),2);
        color_diff = sqrt(r_diff + g_diff + b_diff);
  if(color_diff < threshold) graph(x,y,0,3) = 0;
        else{ graph(x,y,0,3) = color_diff;}
      }
      else{ graph(x,y,0,3) = -1;}
    }
  }
  return graph;
}

// 3.1
// using kruskal's algorithm to sort the input via segmentations produced
std::vector<int> prod_segmentation(CImg<double> input) {
  int i, x, pid1, pid2, tot_pix, prev_grp_num;
  double sx;
  tuple cur_tup;
  std::map<std::pair<int, int>, double> graph;
  std::vector<tuple> out_vec;
  std::pair<int,int> cur_pair;
  tot_pix = input.width() * input.height();
  graph = graph_generator(input);
  std::vector<tuple> sort_vec = sort_graph(graph);  
  sx = calc_Threshold(sort_vec);
  std::vector<int> grouping;
  grouping.resize(tot_pix);
  // put each node in its own group
  for(i = 0; i < tot_pix; i++) {
    grouping[i] = i;
  }
  for(i = 0; i < sort_vec.size(); i++) {
    cur_tup = sort_vec[i];
    cur_pair = cur_tup.first;
    pid1 = cur_pair.first;
    pid2 = cur_pair.second;
    // never merge the groups if they are above the average difference between pixels
    if(out_vec.size() > tot_pix || cur_tup.second > sx * 1.5) {
      break;
    }

    // manage merging the groups
    if(grouping[pid1] != grouping[pid2]) {
      out_vec.push_back(cur_tup);
      prev_grp_num = grouping[pid2];
      for(x = 0; x < tot_pix; x++) {
  if(grouping[x] == prev_grp_num) {
    grouping[x] = grouping[pid1];
  }
  //printf("%i ",grouping[x]);
      }
      //printf("\n");
    }   
  }
  /*for(x = 0; x < tot_pix; x++) {
    if(x % input.width() == 0) {
      printf("\n");
    }
    printf("%i ", grouping[x]);
  }
  printf("avg: %f \n", sx);*/
  return grouping;
}

// Thresholding is done by not grouping if the difference is greater than
// the average difference
double calc_Threshold(std::vector<tuple> input) {
  int i;
  double sum = 0;
  double squared_diff = 0;
  double avg;

  for(i = 0; i < input.size(); i++) {
    sum += input[i].second;
  }
  avg = sum / double(i);

  return avg;
}

// generates the initial fully direct graph for an image
std::map<std::pair<int, int>, double> graph_generator(CImg<double> input) {
  int x,y,gx,gy;
  int pid = 0;
  double val;
  std::map<std::pair<int, int>, double> graph;
  std::pair<int, int> cur_pair;
  for(x = 0; x < input.width(); x++) {
    for(y = 0; y < input.height(); y++) {            
      /* add the coordinates below and to the right of the current coordinate */
      cur_pair.first = pid;
      gx = x;
      gy = y + 1;
      // spot below current pixel
      if(gy < input.height()) {
  if(input.spectrum() != 1) {
    val = sqrt(pow(input(x,y,0,0) - input(gx, gy, 0, 0), 2) +
         pow(input(x,y,0,1) - input(gx, gy, 0, 1), 2) +
         pow(input(x,y,0,2) - input(gx, gy, 0, 2), 2));
  } else {
    val = sqrt(pow(input(x,y,0,0) - input(gx, gy, 0, 0), 2));
  }
  cur_pair.second = pid + 1;
  graph[cur_pair] = val;
      }
      // spot right of current pixel
      gx = x + 1;
      gy = y;
      if(gx < input.width()) {
  if(input.spectrum() != 1) {
    val = sqrt(pow(input(x,y,0,0) - input(gx, gy, 0, 0), 2) +
         pow(input(x,y,0,1) - input(gx, gy, 0, 1), 2) +
         pow(input(x,y,0,2) - input(gx, gy, 0, 2), 2));
  } else {
    val = sqrt(pow(input(x,y,0,0) - input(gx, gy, 0, 0), 2));
  }
  cur_pair.second = pid + input.height();
  graph[cur_pair] = val;
      }
      pid++;
    }
  }  
  return graph;
}
// Using datastructure to copy a map to a vector
// http://stackoverflow.com/questions/684475/c-how-to-copy-a-map-to-a-vector

template<class tup>
struct compare: std::binary_function<tup,tup,bool> {
  bool operator()(const tup& t1, const tup& t2) {
    return t1.second < t2.second;
  }
};

// sort the input map into a vector in increasing order
std::vector<tuple> sort_graph(std::map<std::pair<int, int>, double> graph) {
  vector<tuple> output;
  copy(graph.begin(), graph.end(), back_inserter(output));
  sort(output.begin(), output.end(), compare<tuple>());
  return output;
}

// 3.2
std::vector<int> Real_Segmentation(CImg<double> input, double k) {
  int i, x, pid1, pid2, tot_pix, prev_grp_num;
  tuple cur_tup;
  std::map<std::pair<int, int>, double> graph;
  std::map<int, double> max_grp_weight;
  std::vector<tuple> out_vec;
  std::pair<int,int> cur_pair;
  double n_of_cs[2];
  int mag_of_cs[2];
  tot_pix = input.width() * input.height();
  graph = graph_generator(input);
  std::vector<tuple> sort_vec = sort_graph(graph);  
  std::vector<int> grouping;
  grouping.resize(tot_pix);
  // put each node in its own group
  for(i = 0; i < tot_pix; i++) {
    grouping[i] = i;
  }

  for(i = 0; i < sort_vec.size(); i++) {
    cur_tup = sort_vec[i];
    cur_pair = cur_tup.first;
    pid1 = cur_pair.first;
    pid2 = cur_pair.second;

    // don't merge the groups if they are above the average difference between pixels
    if(out_vec.size() > tot_pix) {
      break;
    }
    
    mag_of_cs[0] = calculate_size(grouping, tot_pix, grouping[pid1]);
    mag_of_cs[1] = calculate_size(grouping, tot_pix, grouping[pid2]);
    if(max_grp_weight.find(pid1) == max_grp_weight.end()) {
      n_of_cs[0] = cur_tup.second;
    } else {
      n_of_cs[0] = max_grp_weight[pid1];
    }
    if(max_grp_weight.find(pid2) == max_grp_weight.end()) {
      n_of_cs[1] = cur_tup.second;
    } else {
      n_of_cs[1] = max_grp_weight[pid2];
    }
    // manage merging the groups
    if(grouping[pid1] != grouping[pid2] && 
       cur_tup.second <= min((n_of_cs[0] + k / double(mag_of_cs[0])),
           (n_of_cs[1] + k / double(mag_of_cs[1])))) {
      out_vec.push_back(cur_tup);
      prev_grp_num = grouping[pid2];
      max_grp_weight[pid1] = cur_tup.second;
      max_grp_weight[pid2] = cur_tup.second;
      for(x = 0; x < tot_pix; x++) {
  if(grouping[x] == prev_grp_num) {
    grouping[x] = grouping[pid1];
  }
  //printf("%i ",grouping[x]);
      }
      //printf("\n");
    }   
  }
  /*for(x = 0; x < tot_pix; x++) {
    if(x % input.width() == 0) {
      printf("\n");
    }
    printf("%i ", grouping[x]);
    }*/
  return grouping;
}

// 3.2 a
// calculate the size of component c
int calculate_size(std::vector<int> data, int length, int grp) {
  int i;
  int size = 0;

  for(i = 0; i < length; i++) {
    if(data[i] == grp) {
      size++;
    }
  }

  return size;
}

CImg<int> return_image(std::vector<int> input, int width, int height) {
  CImg<int> output(width, height, 1, 1);
  int x,y,i;
  x = 0;
  y = 0;
  for(i = 0; i < input.size(); i++) {
    if(i != 0 && i % height == 0) {
      x = 0;
      y++;
    }
    output(y,x,0,0) = input[i];
    x++;
  }
  return output;
}

std::set<int> distinct_values(const std::vector<int> input) {
  std::set<int> dist_vals;
  int i;
  for(i = 0; i < input.size(); i++) {
    dist_vals.insert(input[i]);
  }
  
  return dist_vals;
}


CImg<double> Final_segmentation(std::vector<int> input, int width, int height) {
  int x, y;
  double clrs_left;
  int i = 0;
  int tot_clrs = pow(255, 3);
  int capped = 0;
  clr curr_clr;
  CImg<double> output(width, height, 1, 3);
  std::set<int> distinct = distinct_values(input);
  std::map<int, clr> clr_mapping;
  set<int>::iterator it;
  printf("The number of segments: %i \n", distinct.size());
  clrs_left = (double(tot_clrs) / double(distinct.size()));
  for (it = distinct.begin(); it != distinct.end(); it++) {
    curr_clr.r = rand() % 255;
    curr_clr.g = rand() % 255;
    curr_clr.b = rand() % 255;
    i++;
    clr_mapping[*it] = curr_clr;
  }
  // map the vector into a cimg
  x = 0;
  y = 0;
  for(i = 0; i < input.size(); i++) {
    if(i != 0 && i % height == 0) {
      x = 0;
      y++;
    }
    curr_clr = clr_mapping[input[i]];
    output(y,x,0,0) = curr_clr.r;
    output(y,x,0,1) = curr_clr.g;
    output(y,x,0,2) = curr_clr.b;
    x++;
  }
  return output;
}

// proportionally resize an image to the given dimension
CImg<double> Final_resize(CImg<double> input, int dim) {
  CImg<double> output(input);
  double ratio, w, h;
  if(input.width() > dim || input.height() > dim) {
    if(input.width() > input.height()) {
      ratio = dim / double(input.width());
    } else {
      ratio = dim / double(input.height());
    }
    w = input.width() * ratio;
    h = input.height() * ratio;
    output.resize((int)w, (int)h);
  }
  return output.normalize(0,255);
}

int main(int argc, char **argv)
{
  try {
    if(argc < 4)
      {
        cout << "wrong number of arguments; correct usage:" << endl;
        cout << " p3 problemID inputfile outputfile" << endl;
        return -1;
      }
    string part = argv[1];
    string inputFile = argv[2];
    string outputFile = argv[3];
    string inputFile2;
    double ratio;
    cout << "In: " << inputFile <<" Out: " << outputFile << endl;
    
    CImg<double> input_image(inputFile.c_str());
    CImg<double> output_image;
    std::vector<int> output_vec;
    if(part == "2.1") {
      /* Edge detector:
* Estimates partial derivatives in the x and y directions using the Sobel operator
* at each pixel, and a threshhold to produce a binary image.
*/
      if(input_image.spectrum() != 1) {
input_image = averageGrayscale(input_image);
      }
      output_image = edge_detector(input_image);    
      output_image.save(outputFile.c_str());
    }
    else if(part == "2.2") {
    if(input_image.spectrum() != 1) {
input_image = averageGrayscale(input_image);
      }
      CImg<double> sobel = edge_detector(input_image);
      input_image = edgeusing2ndder(input_image,sobel);
     // input_image = edgeusing2ndder(input_image,sobel);
      input_image.save("input1.jpg");
      CImg<double> x_deriv("sx.jpg");
      CImg<double> y_deriv("sy.jpg");
     output_image = circle_detector(input_image, x_deriv, y_deriv);    
      output_image.save(outputFile.c_str());      
     } else if(part == "3.1") {
      // Simple segmentation algorithm
      input_image = Final_resize(input_image, 300);
      output_vec = prod_segmentation(input_image);
      output_image = Final_segmentation(output_vec, input_image.width(), input_image.height());
      output_image.save(outputFile.c_str());
    } else if(part == "3.2") {
      /* Improved segmentation algorithm */
      // apply gaussian filter to all color layers
      input_image = Final_resize(input_image, 300);
      output_image = input_image;
      CImg<double> gaus;
      double k_param = 500;
      if(input_image.spectrum() != 1) {
  gaus = gaussianClrFilterGen(5);
      } else {
  gaus = Gaussian_Generator(5);
      }
      if(argc > 4) {
  k_param = atof(argv[4]);
      }
      output_image = input_image.convolve(gaus);
      output_image.normalize(0,255);

      // use segmentation with the specified formula
      output_vec = Real_Segmentation(output_image, k_param);
      output_image = Final_segmentation(output_vec, input_image.width(), input_image.height());
      output_image.save(outputFile.c_str());
    }
      else if(part == "1.0") {
      output_image = averageGrayscale(input_image);
      output_image.save(outputFile.c_str());
    }
  } catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
  return 0;
}