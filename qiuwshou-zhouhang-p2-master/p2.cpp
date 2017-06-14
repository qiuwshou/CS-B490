/*
	B490/B659 Project 2 Skeleton Code    (2/2015)
	
	Be sure to read over the project document and this code (should you choose to use it) before 
	starting to program. 
	

	Compiling:
		A simple console command 'make' will excute the compilation instructions
		found in the Makefile bundled with this code. It will result in an executable
		named p2.

	Running:
		The executable p2 should take commands of the form:
			./p2 problem_ID input_File ouput_File additional_Arguments
	
*/


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fft.h>
#include <math.h>
#include <cmath>
//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;




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


  //printf("\n");
  return output;
}



CImg<double> gaussianFilter(double sigma) {


  //FIXME Determine filter size, see Part 4 2b for how
  int filterSize = 3*sigma;

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
  //printf("GAUSSIAN KERNEL VALUES:\n");
   for (int i = 0; i < filterSize; i++) {
    for (int j = 0; j < filterSize; j++) {
      H(i,j,0,0) = H(i,j,0,0) / sum;
      //  printf("\t%f", H(i,j,0,0));
    }
    //   printf("\n");
   }


  // Convolve with filter and return
  //return filter(input, H);
   return H;
}



// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//


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
CImg<double> resize(CImg<double>image){
 
  CImg<double>output(image.width(),image.height(),1,1);
  for(int x=0; x < image.width(); x++){
    for(int y = 0; y< image.height(); y++){
      output(x,y,0,0) = image(x,y,0,0);
    }
  }  


  int size;  
  if(image.width()>image.height()){
    size = image.width();
  }
  else{
    size = image.height();
  }

  int intpart;
  intpart = ceil(log(size)/log(2));
  size =pow(2,intpart);
  output.resize(size,size,1,1,0);

  return output;
};



void fft(const CImg<double> &input, CImg<double> &fft_real, CImg<double> &fft_imag)
{
  fft_real = input;
  fft_imag = input;
  fft_imag = 0.0;

  FFT_2D(1, fft_real, fft_imag);
};

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const CImg<double> &input_real, const CImg<double> &input_imag, CImg<double> &output_real)
{
  output_real = input_real;
  CImg<double> output_imag = input_imag;

  FFT_2D(0, output_real, output_imag);
};



// Write this in Part 3.1
CImg<double> fft_magnitude(const CImg<double> &fft_real, const CImg<double> &fft_imag){
  CImg<double>magnitude(fft_real, "xyzc", 0);
  for(int x = 0; x < fft_real.width(); x++){
    for(int y = 0; y <fft_real.height(); y++){
      magnitude(x,y,0,0) = log(sqrt(fft_real(x,y,0,0)*fft_real(x,y,0,0)+fft_imag(x,y,0,0)*fft_imag(x,y,0,0)));
    }
  }
  return magnitude;
};

// Write this in Part 3.2
CImg<double> remove_interference(const CImg<double> &input){

  CImg<double>real(input, "xyzc", 0);
  CImg<double>imag(input, "xyzc", 0);
  fft(input, real, imag);
  
  CImg<double>output(input, "xyzc",0);

  int x_center = input.width()/2;
  int y_center = input.height()/2;

  for(int x = 0; x < input.width(); x++){
    for(int y = 0; y < input.height(); y++){
      if(x< x_center-90 or x> x_center+90){
        real(x,y,0,0) = 0;
	imag(x,y,0,0) = 0;}
      if(y< y_center-90 or y> y_center+90){
	real(x,y,0,0) = 0;
	imag(x,y,0,0) = 0;}
     
     }
  }
  ifft(real,imag,output);
  output.normalize(0,255);
  return output;
};

// Write this in Part 3.3
CImg<double> fft_filter(const CImg<double> &input, const CImg<double> &filter){
  CImg<double>output(input, "xyzc", 0);
  CImg<double>output_real(input, "xyzc", 0);
  CImg<double>output_imag(input, "xyzc", 0);

  CImg<double>filter_real(input, "xyzc", 0);
  CImg<double>filter_imag(input, "xyzc", 0);
  CImg<double>input_real(input, "xyzc", 0);
  CImg<double>input_imag(input, "xyzc", 0);  

  fft(input, input_real, input_imag);
  fft(filter, filter_real, filter_imag);

  for(int x=0; x<input.width(); x++){
    for(int y= 0; y<input.height(); y++){
      double a = input_real(x,y,0,0);
      double b = input_imag(x,y,0,0);
      double c = filter_real(x,y,0,0);
      double d = filter_imag(x,y,0,0);
      //complex multip
      output_real(x,y,0,0)= a*c-b*d;
      output_imag(x,y,0,0)= c*b+a*d;
      }
  }
  
  
  ifft(output_real, output_imag, output);
  output.normalize(0,255);
  return output;
}
;

// Write this in Part 3.4
CImg<double> fft_defilter(const CImg<double> &input, const CImg<double> &filter){
  CImg<double>output(input, "xyzc", 0);
  CImg<double>output_real(input, "xyzc", 0);
  CImg<double>output_imag(input, "xyzc", 0);

  CImg<double>filter_real(input, "xyzc", 0);
  CImg<double>filter_imag(input, "xyzc", 0);
  CImg<double>input_real(input, "xyzc", 0);
  CImg<double>input_imag(input, "xyzc", 0);

  fft(input, input_real, input_imag);
  fft(filter, filter_real, filter_imag);

  for(int x=0; x<input.width(); x++){
    for(int y= 0; y<input.height(); y++){
      double a = input_real(x,y,0,0);
      double b = input_imag(x,y,0,0);
      double c = filter_real(x,y,0,0);
      double d = filter_imag(x,y,0,0);
      // complex division
      if(c*c+d*d != 0){
        output_real(x,y,0,0)= (a*c+b*d)/(c*c+d*d);
        output_imag(x,y,0,0)= (c*b-a*d)/(c*c+d*d);}
      else{
	output_real(x,y,0,0)= 0;
	output_imag(x,y,0,0)= 0;}
    }
  }

  ifft(output_real, output_imag, output);
  output.normalize(0,255);
  return output;

};

// Write this in Part 4 -- add watermark N to image
CImg<double> mark_image(const CImg<double> &input, int N){

  CImg<double> input_real(input, "xyzc", 0);
  CImg<double> input_imag(input, "xyzc", 0);
  fft(input,input_real, input_imag);
  
  vector<int> binary;

 
  while(N > 0){
    int remainder = N%2;
    //   cout<< remainder<<endl;
    binary.push_back(remainder);   
    N= N/2;
  }

  double alpha = 10;
  double r = input.width()/4;
  int l = binary.size();

  int x_center = input.width() / 2;
  int y_center = input.height() / 2;

  for(int i = 0;i < l;i++){
    double theta =2* (i+1)/l;
    int x1 = x_center + floor(r*cos(theta * M_PI));
    int y1 = y_center + floor(r*sin(theta * M_PI));
    int x2 = x_center - floor(r*cos(theta * M_PI));
    int y2 = y_center - floor(r*sin(theta * M_PI));
    input_real(x1,y1,0,0) = input_real(x1,y1,0,0) + alpha*(abs(input_real(x1,y1,0,0)))*binary[i];
    input_real(x2,y2,0,0) = input_real(x2,y2,0,0) + alpha*(abs(input_real(x2,y2,0,0)))*binary[i];

  }

  CImg<double> output(input, "xyzc",0);
  ifft(input_real, input_imag, output);
  output.normalize(0,255);
  return output;

};
// Write this in Part 4 -- check if watermark N is in image
CImg<double> check_image(const CImg<double> &input, int N) {

  vector<int> v;
  while(N > 0){
    int remainder = N%2;
    //   cout<< remainder<<endl;
    v.push_back(remainder);
    N= N/2;
  }

  CImg<double> input_real(input, "xyzc", 0);
  CImg<double> input_imag(input, "xyzc", 0);
  fft(input,input_real, input_imag);

  vector<double> c;

  double alpha = 20;
  double r = input.width()/4;
  int l = v.size();

  int x_center = input.width() / 2;
  int y_center = input.height() / 2;

  //since Vi = 0 or 1, so it is easy to set Ci between 0 to 1
  input_real.normalize(0,1); 
  for(int i = 0;i < l;i++){
    double theta = 2* (i+1)/l;
    int x1 = x_center + floor(r*cos(theta * M_PI));
    int y1 = y_center + floor(r*sin(theta * M_PI));
    double real = input_real(x1,y1,0,0);
    c.push_back(real);
  }

  //if V and C have correlation, Vi=1 gives larger Ci, Vi=0 gives smaller Ci
  // since we normalize Ci so when Vi=0, Ci will be close to 0
  // so we add up value of Ci/Vi when Vi=1 
  //and threshold t depends on the size of Vi 
  double sum = 0;
  for(int i = 0; i < l ; i++){
    if(v[i] == 1){
      sum = sum + c[i];
    }
  }
  cout << sum << endl;
  if(sum > 0.9){
    cout<< "there is a watermark" << endl;
  }
  else{
    cout<< "there is no watermark" << endl; 
  
  }

};


int main(int argc, char **argv)
{ 
  srand(time(NULL));
  try {

    if(argc < 4)
      {
	cout << "Insufficent number of arguments; correct usage:" << endl;
	cout << "    p2 problemID inputfile outputfile" << endl;
	return -1;
      }
    
    string part = argv[1];
    string inputFile = argv[2];
    string outputFile = argv[3];
    cout << "In: " << inputFile <<"  Out: " << outputFile << endl;
    CImg<double>input_image(inputFile.c_str());
    input_image  = gray(input_image);
    input_image = resize(input_image);
   
    
    if(part == "2.1")
      { 
	//read 2 images as input
	CImg<double>input1(inputFile.c_str());
	input1 = gray(input1);
	CImg<double>input2(outputFile.c_str());
        input2 = gray(input2);
        input2.resize(input1.width(), input1.height(),1,1,0);
        CImg<double>low = filter(input1,gaussianFilter(5.0));
        CImg<double>low2 = filter(input2,gaussianFilter(5.0));
        
	low.save("low1.jpg");
	low2.save("low2.jpg");
	CImg<double>result(input1, "xyzc", 0);
	for(int x = 0; x <input1.width(); x++){
	  for(int y = 0; y<input1.height(); y++){
	    result(x,y,0,0) = low(x,y,0,0) + input2(x,y,0,0)- low2(x,y,0,0);
	  }
	}
	result.normalize(0,255);    
	result.save("morphing.png");
      }
    if(part == "3.1")
      { 
       	CImg<double>real(input_image, "xyzc", 0);
        CImg<double>imag(input_image, "xyzc", 0);
        fft(input_image, real, imag);

        CImg<double>magnitude= fft_magnitude(real, imag);
        magnitude.normalize(0,255);
        magnitude.save(outputFile.c_str());
      }
    if(part == "3.2")
      { //input = noise1.png 
       	input_image = remove_interference(input_image);
        input_image.save(outputFile.c_str());
      }
    if(part == "3.3")
      { 
	CImg<double>filter = gaussianFilter(5.0);
        filter.resize(input_image.width(),input_image.height(),1,1,0);
        
        CImg<double> output(input_image, "xyzc", 0);
	output = fft_filter(input_image, filter);
        output.save(outputFile.c_str());

      }
    if(part == "3.4")
      { //input from the result of 3.3
	CImg<double>filter = gaussianFilter(5.0);
        filter.resize(input_image.width(),input_image.height(),1,1,0);
        
        CImg<double> output(input_image, "xyzc",0);
	output = fft_defilter(input_image, filter);
	output.save(outputFile.c_str());
      }
    // add more if conditions for other parts
    int N = rand()%input_image.width()+1;

    if(part == "4.1"){

      input_image = mark_image(input_image, N);
      input_image.save(outputFile.c_str());

    }

    if(part == "4.2"){

      input_image = check_image(input_image, N);

    }
  } 
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}








