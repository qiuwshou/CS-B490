/*
   B490/B659 Project 1 Skeleton Code    (1/2015)

   Be sure to read over the project document and this code (should you choose to use it) before 
   starting to program. 


Compiling:
A simple console command 'make' will excute the compilation instructions
found in the Makefile bundled with this code. It will result in an executable
named p1.

Running:
The executable p1 takes commands of the form:
./p1 problem_ID input_File ouput_File additional_Arguments

Some examples:

./p1 2.1 input.png out.png 

This runs the 'averageGrayscale' function defined in this file and described
in the project doc Part 2, problem 1 on the input. The output is saved as out.png. 

----

./p1 4.2b input.jpg output.png 0.5

This runs the Gaussian filter function (gaussianFilter) with sigma = 0.5 on input.jpg.
This is problem 2b from Part 4 in the project documentation. 

*/


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;


//Part 2 - Basic image operations
CImg<double> averageGrayscale(CImg<double> input);
CImg<double> simpleBW(CImg<double> input);
CImg<double> advancedBW(CImg<double> input);

//Part 3 - Adding noise
CImg<double> uniformNoise(CImg<double> input);
CImg<double> gaussianNoise(CImg<double> input, double sigma);
CImg<double> saltAndPepperNoise(CImg<double> input);

//Part 4 - Filtering
CImg<double> filter(CImg<double> input, CImg<double> filter);
CImg<double> meanFilter(CImg<double> input, int filterSize);
CImg<double> gaussianFilter(CImg<double> input, double sigma);
CImg<double> medianFilter(CImg<double> input, int size);

double meanSquaredError(CImg<double> original, CImg<double> other);

int main(int argc, char **argv){


   if(argc < 4){
      cout << "Insufficent number of arguments. Please see documentation" << endl;
      cout << "p1 problemID inputfile outputfile" << endl;
      return -1;
   }


   char* inputFile = argv[2];
   char* outputFile = argv[3];
   cout << "In: " << inputFile <<"  Out: " << outputFile << endl;

   CImg<double> input(inputFile);
   CImg<double> output; 

   if(!strcmp(argv[1], "2.1")){  
      cout << "# Problem 2.1 - Average Grayscale" << endl;
      if(input.spectrum() != 3){ cout << "INPUT ERROR: Input image is not a color image!" << endl;return -1;}
      output = averageGrayscale(input);
      output.save(outputFile);
   }
   else if(!strcmp(argv[1], "2.2a")){
      cout << "# Problem 2.1a - Simple Threshold Black and White" << endl;
      if(input.spectrum() != 1){ cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
      output = simpleBW(input);
      output.save(outputFile);
   }
   else if(!strcmp(argv[1], "2.2b")){
      cout << "# Problem 2.2b - Advanced Threshold Black and White" << endl;
      if(input.spectrum() != 1){ cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
      output = advancedBW(input);
      output.save(outputFile);
   }
   else if(!strcmp(argv[1], "3.1")){
      cout << "# Problem 3.1 - Uniform Noise" << endl;
      if(input.spectrum() != 1){ cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
      output = uniformNoise(input);
      output.save(outputFile);
   }
   else if(!strcmp(argv[1], "3.2")){
      cout << "# Problem 3.2 - Gaussian Noise" << endl;
      if(input.spectrum() != 1){ cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
      if(argc != 5){ cout << "INPUT ERROR: Provide sigma as additional argument!" << endl;return -1;}
      output = gaussianNoise(input, atof(argv[4]));
      output.save(outputFile);
   }
   else if(!strcmp(argv[1], "3.3")){
      cout << "# Problem 3.3 - Salt & Pepper Noise" << endl;
      if(input.spectrum() != 1){ cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
      output = saltAndPepperNoise(input);
      output.save(outputFile);
   }
   else if(!strcmp(argv[1], "4.2a")){
      cout << "# Problem 4.2a - Mean Filter Noise" << endl;
      if(input.spectrum() != 1){ cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
      if(argc != 5){ cout << "INPUT ERROR: Provide filter size as additional argument!" << endl;return -1;}
      output = meanFilter(input, atoi(argv[4]));
      output.save(outputFile);
   }
   else if(!strcmp(argv[1], "4.2b")){
      cout << "# Problem 4.2b - Gaussian Filter Noise" << endl;
      if(input.spectrum() != 1){ cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
      if(argc != 5){ cout << "INPUT ERROR: Provide sigma as additional argument!" << endl;return -1;}
      output = gaussianFilter(input, atof(argv[4]));
      output.save(outputFile);
   }
   else if(!strcmp(argv[1], "4.3")){
      cout << "# Problem 4.3 - Median Noise" << endl;
      if(input.spectrum() != 1){ cout << "INPUT ERROR: Input image is not a grayscale image!" << endl;return -1;}
      if(argc != 5){ cout << "INPUT ERROR: Provide filter size as additional argument!" << endl;return -1;}
      output = medianFilter(input, atoi(argv[4]));
      output.save(outputFile);
   }
   else if(!strcmp(argv[1], "5")){

      cout << "# Problem 5 - Noise Removal Analysis" << endl;
      //FIXME You will need to implement this section yourself
     
      CImg<double> cat_gray("cat_gray.jpg");
       
      // Add noise to grayscale image
      CImg<double> in_uniformNoise = uniformNoise(input);
      CImg<double> in_gaussianNoise = gaussianNoise(input, 2.5);
      CImg<double> in_saltAndPepperNoise = saltAndPepperNoise(input);
     
      in_uniformNoise.save("noise_uniform.jpg");
      in_gaussianNoise.save("noise_gaussian.jpg");
      in_saltAndPepperNoise.save("noise_SP.jpg");
 
      int meanFilterSize = 5;
      double gaussianSigma = 1.8;
      int  medianFilterSize = 5;

      CImg<double> images[9];
      images[0] = meanFilter(in_uniformNoise, meanFilterSize);
      images[1] = meanFilter(in_gaussianNoise, meanFilterSize);
      images[2] = meanFilter(in_saltAndPepperNoise, meanFilterSize),

      images[3] = gaussianFilter(in_uniformNoise, gaussianSigma);
      images[4] = gaussianFilter(in_gaussianNoise, gaussianSigma);
      images[5] = gaussianFilter(in_saltAndPepperNoise, gaussianSigma);

      images[6] = medianFilter(in_uniformNoise, medianFilterSize);
      images[7] = medianFilter(in_gaussianNoise, medianFilterSize);
      images[8] = medianFilter(in_saltAndPepperNoise, medianFilterSize);
 
      printf("-----------------------------------------------------------------------------\n");

      char noise_names[3][50] = {"UNIFORM ", "GAUSSIAN", "S & P   "};
      char filter_names[3][50] = {"MEAN    ", "GAUSSIAN", "MEDIAN  "};
      char out_names[9][9] = {"out1.jpg","out2.jpg","out3.jpg",
                              "out4.jpg","out5.jpg","out6.jpg",
                              "out7.jpg","out8.jpg","out9.jpg"};
      for (int i = 0; i < 9; i++) {
         char *noise_name, *filter_name;
         if (i < 3) filter_name = filter_names[0];
         else if (i < 6) filter_name = filter_names[1];
         else filter_name = filter_names[2];

         noise_name = noise_names[i % 3];
         double error = meanSquaredError(input, images[i]);
         printf("[%s]   %s   filter over   %s   noise   ERROR: %f\n", out_names[i], filter_name, noise_name, error);

         images[i].save(out_names[i]);
      
      } 
         
      printf("-----------------------------------------------------------------------------\n");
   
   } else if(!strcmp(argv[1], "6.1")){
      cout << "# Problem 6.1 - Separable Kernel Convolutions" << endl;

      std::clock_t start = std::clock();
      //FIXME You will need to implement this section yourself
      std::cout << "Time Elapsed: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

   }
   else if(!strcmp(argv[1], "6.2")){
      cout << "# Problem 6.2 - Dynamic Box Filter" << endl;

      std::clock_t start = std::clock();
      //FIXME You will need to implement this section yourself
      std::cout << "Time Elapsed: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

   }
   else if(!strcmp(argv[1], "6.3")){
      cout << "# Problem 6.3 - Fast Gaussian Smoothing" << endl;

      std::clock_t start = std::clock();
      //FIXME You `will need to implement this section yourself
      std::cout << "Time Elapsed: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

   }
   else{
      cout << "Unknown input command" << endl; 
   }



   return 0;
}


// PART 2 - BASIC IMAGE OPERATIONS


// Problem 2.1
CImg<double> averageGrayscale(CImg<double> input){
   //Creates a new grayscale image with same size as the input initialized to all 0s (black)
   CImg<double> output(input.width(), input.height(), 1, 1); 

   //FIXME fill in new grayscale image
   for (int i=0; i<input.width(); i++) {
      for (int j=0; j<input.height(); j++) {
         int r = input(i,j,0,0);
         int g = input(i,j,0,1);
         int b = input(i,j,0,2);
         int avg = (r+g+b) / 3;
         output(i,j,0,0) = avg;
      }
   }

   return output;
}

// Problem 2.2a
CImg<double> simpleBW(CImg<double> input){
   //Creates a new grayscale image with same size as the input initialized to all 0s (black)
   CImg<double> output(input, "xyzc", 0);  

   //FIXME Do stuff
   for (int i=0; i<input.width(); i++) {
      for (int j=0; j<input.height(); j++) {
         int v = input(i,j,0,0);
         if (v > 127)
            output(i,j,0,0) = 255;
         else
            output(i,j,0,0) = 0;
      }
   } 

   return output;
}

// Problem 2.2b
CImg<double> advancedBW(CImg<double> input){
   //Creates a new grayscale image with same size as the input initialized to all 0s (black)
   CImg<double> output(input.width(), input.height(), 1, 1); 

   //FIXME Do stuff
   for (int i=0; i < output.width(); i++) {
      for (int j=0; j < output.height(); j++) {

         int v = input(i,j,0,0);
         if (v > 127) v = 255; else v = 0;
         int e = input(i,j,0,0) - v;

         if (i+1 < output.width()) { // looking at next row
            if (j-1 > -1) // bottom left
               input(i+1,j-1,0,0) = input(i+1,j-1,0,0) + 0.2 * e;
            if (j+1 < output.height()) // bottom right
               input(i+1,j+1,0,0) = input(i+1,j+1,0,0) + 0.1 * e;
            input(i+1,j,0,0) = input(i+1,j,0,0) + 0.3 * e; // bottom middle
         }
         if (j+1 < output.height()) // same row, next column
            input(i,j+1,0,0) = input(i,j+1,0,0) + 0.4 * e;

         output(i,j,0,0) = v;
      }
   }

   return output;
}



//Part 3 - Adding noise
CImg<double> uniformNoise(CImg<double> input){

   //Creates a new image with same size as the input initialized to all 0s (black)
   CImg<double> output(input, "xyzc", 0); 

	for (int i = 0; i < input.width(); i++) {
		for (int j = 0; j < input.height(); j++) {
			double s = rand() / (double)RAND_MAX*21;
			double add = rand() / (double)RAND_MAX*2;
			double out_pixel;
			if (add) {
				out_pixel = input(i,j,0,0) + s;
				if (out_pixel > 255) out_pixel = 255.0;
			} else {
				out_pixel = input(i,j,0,0) - s;
				if (out_pixel < 0) out_pixel = 0.0;
			}

			output(i,j,0,0) = out_pixel;
		}
	}
	
   return output;
}

// Compliments of Wikipedia: "Box-Muller Transform"
#define TWO_PI 6.2831853071795864769252866
double generateGaussianNoise(const double mu, const double sigma)
{
	using namespace std;
	static bool haveSpare = false;
	static double rand1, rand2;
 
	if(haveSpare)
	{
		haveSpare = false;
		return (sigma * sqrt(rand1) * sin(rand2)) + mu;
	}
 
	haveSpare = true;
 
	rand1 = rand() / ((double) RAND_MAX);
	if(rand1 < 1e-100) rand1 = 1e-100;
	rand1 = -2 * log(rand1);
	rand2 = (rand() / ((double) RAND_MAX)) * TWO_PI;
 
	return (sigma * sqrt(rand1) * cos(rand2)) + mu;
}

CImg<double> gaussianNoise(CImg<double> input, double sigma){

   //Creates a new image with same size as the input initialized to all 0s (black)
   CImg<double> output(input, "xyzc", 0); 

	for (int i = 0; i < input.width(); i++) {
		for (int j = 0; j < input.height(); j++) {
			double noise = generateGaussianNoise(0.0, sigma);
			double out_pixel = input(i,j,0,0) + noise;
			
			// floor/ceiling
			if (out_pixel > 255)
				out_pixel = 255.0;
			else if (out_pixel < 0)
				out_pixel = 0.0;

			output(i,j,0,0) = out_pixel;
		}
	}

   return output;
}

CImg<double> saltAndPepperNoise(CImg<double> input){

   //Creates a new image with same size as the input initialized to all 0s (black)
   CImg<double> output(input, "xyzc", 0); 

   double p=0.005;
   int humber = input.width()*input.height();

   for(int x = 0; x<input.width();x++){
      for(int y = 0; y<input.height();y++){
         int random = rand()%500;
         if(random == 5){
            int random2 =rand()%2;
            if(random2 == 1){
               output(x,y,0,0)=255;
            }
            else{
               output(x,y,0,0)=0;
            }
         }
         else{
            output(x,y,0,0)=input(x,y,0,0);
         }
      }
   }

   return output;
}

//Part 4 - Filtering

//Problem 4.1
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


// Problem 4.2a
CImg<double> meanFilter(CImg<double> input, int filterSize){

   //Creates a new grayscale image (just a matrix) to be our filter
   CImg<double> H(filterSize, filterSize, 1, 1); 
   
   //FIXME Fill filter values
   
   //printf("MEAN KERNEL VALUES:\n");

   for (int i = 0; i < filterSize; i++) {
      for (int j = 0; j < filterSize; j++) {
         H(i,j,0,0) = 1.0/(filterSize * filterSize); // already normalized
         //printf("\t%f", H(i,j,0,0));
      }
      //printf("\n");
   }
   
   //Convole with filter and return
   return filter(input, H);

}


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
   printf("GAUSSIAN KERNEL VALUES:\n");
   for (int i = 0; i < filterSize; i++) {
      for (int j = 0; j < filterSize; j++) {
         H(i,j,0,0) = H(i,j,0,0) / sum;
         printf("\t%f", H(i,j,0,0));
      }
      printf("\n");
   }

   // Convolve with filter and return
   return filter(input, H);
}


int cmp_double(const void *x, const void *y) {
   double xx = *(double*)x;
   double yy = *(double*)y;
   if (xx < yy) return -1;
   if (xx > yy) return 1;
   return 0;
}

CImg<double> medianFilter(CImg<double> input, int size) {

   // Creates a new image with the same size as the input initialized to all 0s (black)
   CImg<double> output(input, "xyzc", 0);

   //FIXME Do median filtering
   
   int s = size/2;
   for (int i = s; i < input.width() - s; i++) {
      for (int j = s; j < input.height() - s; j++) {
         
         double neighbors[size*size];
         for (int u = 0 - s; u <= s; u++) {
            for (int v = 0 - s; v <= s; v++) {
               neighbors[(u+s)*size + (v+s)] = input(i-u, j-v, 0, 0);
            }
         }

         int n_len = sizeof(neighbors)/sizeof(neighbors[0]); // neighbors length
         qsort(neighbors, n_len, sizeof(neighbors[0]), cmp_double);
         
         double median;
         if (n_len % 2)
            median = (neighbors[n_len/2 -1] + neighbors[n_len/2]) / 2.0;
         else
            median  = neighbors[n_len/2];
      
         output(i,j,0,0) = median;   
      }
   }
   return output; 
}

double meanSquaredError(CImg<double> original, CImg<double> other) {
   if (original.width() != other.width() || original.height() != other.height())
      return -1;
   
   double error = 0.0;
   for (int i = 0; i < original.width(); i++) {
      for (int j = 0; j < original.height(); j++) {
         error += pow(original(i,j,0,0) - other(i,j,0,0), 2);
      }
   }

   return error / (original.width() * original.width());
   
}
