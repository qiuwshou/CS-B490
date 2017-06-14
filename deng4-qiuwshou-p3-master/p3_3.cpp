#include "CImg.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>

using namespace cimg_library;
using namespace std;

CImg<double> gray(CImg<double>image){
  
  CImg<double>gray(image.width(),image.height(),1,1);
  for(int x=0;x<image.width();x++){
    for(int y=0;y<image.height();y++){
      if(image.spectrum() ==3){
	gray(x,y,0,0) =(image(x,y,0,1)+image(x,y,0,0)+image(x,y,0,2))/3;}
      else{
	gray(x,y,0,0) = image(x,y,0,0);}
    }
  }

  return gray;
}


CImg<double> filter(CImg<double> input, CImg<double> filter) {
  CImg<double> output(input, "wyzc", 0);

  int fw = filter.width()/2;
  int fh = filter.height()/2;
   
  for (int i = fw; i < input.width() - fw; i++) {
    for (int j = fh; j < input.height() - fh; j++) {
      double sum = 0.0;
      for (int u = 0 - fw; u <= fw; u++) {
	for (int v = 0 - fh; v <= fh; v++) {
	  sum += filter(u+fw,v+fh,0,0) * input(i-u,j-v,0,0);
	  //printf("[%d][%d]\t[%d][%d]\n", u+fw, v+fh, i-u, j-v);
	}
            
      }
         
      if (sum < 0) sum = 0;
      if (sum > 255) sum = 255;
         
      output(i,j,0,0) = sum;
    }
  } 
  return output;
}

//part2.1
CImg<double> edge_detector(CImg<double> input, const double &threshold){

  CImg<double> y_deriv(input); 
  CImg<double> x_deriv(input);
  CImg<double> sx(3, 3, 1, 1);
  CImg<double> sy(3, 3, 1, 1);

  sx(0,0,0,0) = -.125;  sx(1,0,0,0) = 0.00;  sx(2,0,0,0) = .125;
  sx(0,1,0,0) = -.25;   sx(1,1,0,0) = 0.00;  sx(2,1,0,0) = .25;
  sx(0,2,0,0) = -.125;  sx(1,2,0,0) = 0.00;  sx(2,2,0,0) = .125;

  sy(0,0,0,0) = .125;   sy(1,0,0,0) = 0.25;  sy(2,0,0,0) = .125;
  sy(0,1,0,0) = 0.00;   sy(1,1,0,0) = 0.00;  sy(2,1,0,0) = 0.00;
  sy(0,2,0,0) = -.125;  sy(1,2,0,0) = -.25;  sy(2,2,0,0) = -.125;

  x_deriv = filter(input, sx);
  y_deriv = filter(input, sy);

  x_deriv.normalize(0,255);
  y_deriv.normalize(0,255);

  CImg<double> output(input,"wyzc",0);

  for(int i = 0; i < input.width(); i++){
    for(int j = 0; j < input.height(); j++){
      output(i,j,0,0) = sqrt( pow(x_deriv(i,j,0,0), 2) + pow(y_deriv(i,j,0,0), 2));
      if(output(i,j,0,0) > threshold) output(i,j,0,0) = 255;
      else output(i,j,0,0) = 0;
    }
  }
  return output;
}

//part2.2
//increment all circles that pass through (x,y)
void increment(int x, int y,  CImg<double> &houghSpace){
  for(int i = 0; i < houghSpace.width(); i++){
    for(int j = 0; j < houghSpace.height(); j++){
	int local_R = floor(sqrt(pow(i-x,2)+pow(j-y,2)));
	//to check if the circle is out of the boundary
        int x_up = i + local_R - houghSpace.width(); int x_down = i - local_R;
        int y_up = j + local_R - houghSpace.height();int y_down = j - local_R;
        int mult = x_up * x_down * y_up * y_down;
        if( mult > 0) houghSpace(i,j,local_R,0,0) += 1;       

    }
  }
}
 
void drawCircle(int x, int y, int z, CImg<double> &image){
  for(int i = 0; i < image.width(); i++){
    for(int j = 0; j < image.height(); j++){
      int dist = floor(sqrt(pow(x-i,2)+pow(y-j,2)));
      if(dist == z) image(i,j,0,0) = 255;
    }
  }
}
//part3.1
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
//part3.2
CImg<double> gaussianFilter(CImg<double> input, double sigma) {
  //FIXME Determine filter size, see Part 4 2b for how
  int filterSize = sigma * 3;
  if (!(filterSize % 2))
    filterSize -= 1;
  //Creates a new grayscale image (just a matrix) to be our filter 
  CImg<double> H(filterSize, filterSize, 1, 1); 
  //FIXME Fill filter values
  double sum = 0.0;
  for (int x = 0; x < filterSize; x++) {
    for (int y = 0; y < filterSize; y++) {
      double mean = filterSize/2;
      double g = exp(-0.5 * (pow((x-mean)/sigma, 2.0) + pow((y-mean)/sigma, 2.0))) 
	/ (2 * M_PI * sigma * sigma);
         
      H(x,y,0,0) = g;
      sum += H(x,y,0,0);
    }
  }
  // Normalize kernel values
  // printf("GAUSSIAN KERNEL VALUES:\n");
  for (int i = 0; i < filterSize; i++) {
    for (int j = 0; j < filterSize; j++) {
      H(i,j,0,0) = H(i,j,0,0) / sum;
      //  printf("\t%f", H(i,j,0,0));
    }
    //  printf("\n");
  }
  // Convolve with filter and return
  return filter(input, H);
}





int main(int argc, char **argv){ 
  try {
  if(argc < 4)
    {
      cout << "Insufficent number of arguments; correct usage:" << endl;
      cout << "    p3 problemID inputfile outputfile" << endl;
      return -1;
    }

  string part = argv[1];
  string inputFile = argv[2];
  string outputFile = argv[3];
  cout << "In: " << inputFile <<"  Out: " << outputFile << endl;
  CImg<double>input_image(inputFile.c_str());
  if(input_image.spectrum() != 1) input_image = gray(input_image);

  if(part == "2.1")
    {
      CImg<double>result(input_image);
      result = edge_detector(input_image, 50);
      result.save(outputFile.c_str());

    }
  if(part == "2.2")
    {
      //sth wrong for resize!!!!!!!
      input_image.resize(-30,-30,1,1,0);
      input_image.save("111.jpg");
      int R_BIN_SIZE = 1;
      double image_diagonal = sqrt(pow(input_image.width(),2)+pow(input_image.height(),2));
      int R_BINS = image_diagonal/R_BIN_SIZE + 1;
      CImg<double> image(input_image, "xyzc", 0);
      CImg<double> houghSpace(input_image.width(),input_image.height(),R_BINS,1,1);
      
      for(int i = 0; i < input_image.width(); i++){
	for(int j = 0; j < input_image.height(); j++){
	  if(input_image(i,j,0,0)>0)     
	    increment(i, j, houghSpace);
	}
      }     

      int threshold = 50;
      for(int x = 0; x < image.width(); x++){
	for(int y = 0; y < image.height(); y++){
	  for(int z = 1; z <= R_BINS; z++){
    //if a spot get more than 50 votes from other points, draw circle around it
	    if(houghSpace(x,y,z,0,0) > threshold)
	      drawCircle(x, y, z, image);
	  }
	}
      }
      image.normalize(0,255).save(outputFile.c_str());

    }
  if(part == "3.1")
    {
      input_image = segmentation(input_image, 50);
      input_image.save(outputFile.c_str());
    }
  if(part == "3.2")
    { //input is a graph
      input_image = gaussianFilter(input_image, 10);
      double MAX_WEIGHT = -1;
      for(int x = 0; x < input_image.width(); x++){
	for(int y = 0; y < input_image.height(); y++){
	  for(int i = 0; i <= 3; i++){
            if(input_image(x,y,0,i) > MAX_WEIGHT){
	      MAX_WEIGHT = input(x,y,0,i);
	    }
	  }  
	}

    }
  if(part == "3.3")
    {



    }
  if(part == "4")
    {



    }

  } 
catch(const string &err) {
  cerr << "Error: " << err << endl;
 }
}
 
