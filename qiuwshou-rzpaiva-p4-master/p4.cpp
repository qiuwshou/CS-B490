// B490/B659 Project 4 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

//part2
CImg<double> warping(CImg<double> input){
  //use inverse matrix
  CImg<double>matrix(3,3,1,1);
  matrix(0,0,0,0) = 1.12; matrix(1,0,0,0) = -0.315; matrix(2,0,0,0) = 222;
  matrix(0,1,0,0) = 0.109; matrix(1,1,0,0) = 0.685; matrix(2,1,0,0) = -19.9;
  matrix(0,2,0,0) = 0.000265; matrix(1,2,0,0) = -0.000597; matrix(2,2,0,0) = 1.08;
  //use a larger image for output
  CImg<double>output(2*input.width(),2*input.height(),1,3);

  for(int x = 0; x < output.width(); x++){
    for(int y = 0 ; y < output.height(); y++){

      double x_prime = matrix(0,0,0,0)*x + matrix(1,0,0,0)*y + matrix(2,0,0,0);
      double y_prime = matrix(0,1,0,0)*x + matrix(1,1,0,0)*y + matrix(2,1,0,0);
      double w_prime = matrix(0,2,0,0)*x + matrix(1,2,0,0)*y + matrix(2,2,0,0);
      int x2 = (x_prime/w_prime);
      int y2 = (y_prime/w_prime);


      if( x2 >= 0 and x2 <= input.width() and y2 >=0 and y2 <= input.height()){
        output(x,y,0,0) = input(x2,y2,0,0);
        output(x,y,0,1) = input(x2,y2,0,1);
        output(x,y,0,2) = input(x2,y2,0,2);
      }
      else{
        output(x,y,0,0) = 255;
        output(x,y,0,1) = 255;
        output(x,y,0,2) = 255;
      }
    }
  }
  return output;
}
// part 3 distance
float D(vector<float> & d1, vector<float> &d2)
{
	float a =0.0;
	float error =0.0;
	for(int i =0; i<128;i++)
	{
		a = (float)d1[i]-(float)d2[i];
		a *= a;
		error+=a;
	}
	error = sqrt(error);
	return error;

}

//part 3 compare 2 descriptors
typedef struct
{
	int d1,d2;
	float vec_x, vec_y;
}desc_pairs;
void compareDescriptors2(vector<SiftDescriptor>& v1, vector<SiftDescriptor>& v2,  float ratio,vector<desc_pairs>&pairs)
{
	float lowest=10000.0;
	float lowest2=10000.0;
	int lowestid=-1,lowestid2=-1,count=0;
	float dis =0.0;
	desc_pairs a;
	for(int i = 0; i< v1.size();i++)
	{
		lowest = 10000.0;
		lowest2 = 10000.0;
		lowestid=-1;
		lowestid2=-1;
		for(int j = 0; j< v2.size();j++)
		{
			dis = D(v1[i].descriptor,v2[j].descriptor);
			if(dis<lowest){
				lowest2 = lowest;
				lowestid2 = lowestid;
				lowest = dis;
				lowestid =j;

			}
		}
		if(lowestid !=-1 && lowestid2!=-1 )
		{
			
			float dis = lowest/lowest2;
			if(dis<ratio)
			{
				a.d1 = i; a.d2 = lowestid;
				a.vec_x = (v2[lowestid].col-v1[i].col);
				a.vec_y = (v2[lowestid].row-v1[i].row);
				pairs.push_back(a);
			}
	//					cout<< "best match for #"<<i<<" is #"<<lowestid<<" with"<<lowest;
	//					cout<< "second is #"<<lowestid2<<" with"<<lowest2<<endl;
		}
	}


}
//compare descriptors with threshold
void compareDescriptors(vector<SiftDescriptor>& v1, vector<SiftDescriptor>& v2,float threshold, vector<desc_pairs>&pairs)
{
	float lowest=10000.0;
	int lowestid=-1,count=0;
	float dis =0.0;
	desc_pairs a;
	for(int i = 0; i< v1.size();i++)
	{
		lowest = 10000.0;
		lowestid=-1;
		for(int j = 0; j< v2.size();j++)
		{
			dis = D(v1[i].descriptor,v2[j].descriptor);
			if(dis<lowest){
				lowest = dis;
				lowestid =j;

			}
		}
		if(lowestid !=-1 && lowest< threshold)
		{
			a.d1 = i; a.d2 = lowestid;
			pairs.push_back(a);
			//			cout<< "best match for #"<<i<<" is #"<<lowestid<<" with"<<lowest<<endl;
		}
	}


}
typedef struct
{
	int count;
	string filename;
}image_count;

void getTranslation2(vector<desc_pairs> &pairs, float &vec_x, float &vec_y)
{
	srand(time(NULL));
	int best_e= -1;
	int best_count =0;
	int count = 0;
	for(int k = 0;k<floor(pairs.size()*0.3); k++)
	{
		count = 0;
		int e =  rand()%pairs.size();
		vec_x=pairs[e].vec_x;
		vec_y=pairs[e].vec_y;
		for(int i=0;i<pairs.size();i++)
		{
			if(abs(vec_x - pairs[i].vec_x) < 2 && abs(vec_y - pairs[i].vec_y )<2)
			{	
				count ++;	
			}
		}
		if(count > best_count)
		{
			best_count = count;
			best_e = e;
		}
	}
	if(best_e !=-1)
	{
		vec_x = pairs[best_e].vec_x;
		vec_y = pairs[best_e].vec_y;
		cout<<" best : "<<best_e<< " "<<vec_x<<" "<<vec_y<<endl;
	}
	
}
void getTranslation(vector<desc_pairs> &pairs, float &vec_x, float &vec_y)
{
	vec_x=vec_y=0;
	for(int i=0;i<pairs.size();i++)
	{
		vec_x += pairs[i].vec_x;
		vec_y += pairs[i].vec_y;
	}
	vec_x/=pairs.size();
	vec_y/=pairs.size();
	
}

bool compare(const image_count&a,const image_count &b)
{
	return a.count>b.count;
}
int main(int argc, char **argv)
{
	try {

		if(argc < 2)
		{
			cout << "Insufficent number of arguments; correct usage:" << endl;
			cout << "    p4 part_id ..." << endl;
			return -1;
		}

		string part = argv[1];
		string inputFile = argv[2];

		if(part == "part2")
		{
		  CImg<double> input(inputFile.c_str());
		  CImg<double> output(input,"xyzc",0);
		  output = warping(input);
		  output.save("2.png");
		}
		else if(part == "part3")
		{
			// This is just a bit of sample code to get you started, to
			// show how to use the SIFT library.
			if(argc < 5)
			{
				cout << "Insufficent number of arguments; correct usage:" << endl;
				cout << "    p4 part_id threshold query img1..." << endl;
				if(argc >= 4)
					cout << "Missing at least one other img " << endl;
				return -1;
			}
			float threshold = atof(argv[2]);
			string inputFile = argv[3];

			CImg<double> query_image(inputFile.c_str());
			// convert image to grayscale
			CImg<double> query_gray = query_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> query_descriptors = Sift::compute_sift(query_gray);
			vector<SiftDescriptor> descriptors = vector<SiftDescriptor>();
			cout<< "query file "<<inputFile<<endl; 
			vector<image_count> images= vector<image_count>();
			vector<desc_pairs> final = vector<desc_pairs>();
			vector<desc_pairs> temp = vector<desc_pairs>();
			time_t start , stop;
			for(int i= 4; i< argc; i++)
			{
				cout<< "comparing with : "<<inputFile; 
				time(&start);
				inputFile = argv[i];

				CImg<double> input_image(inputFile.c_str());
				if(input_image.spectrum()!=3)
				{
					cout<<"skipping. not rgb"<<endl;continue;
				}
				CImg<double> gray; // convert image to grayscale
				gray = input_image.get_RGBtoHSI().get_channel(2);
				descriptors.clear(); temp.clear();
				descriptors = Sift::compute_sift(gray);
				compareDescriptors(query_descriptors,descriptors,threshold,temp);
				if( final.size() < temp.size())
				{
					final.clear();
					final = temp;
				}
				image_count c; c.count = temp.size(); c.filename =inputFile;
				images.push_back(c);
				time(&stop);
				cout<<" - done! "<<(double) difftime(stop,start)<<" sec ("<<descriptors.size()<<")"<<endl;

			}
			sort(images.begin(),images.end(),compare);
			for(int i=0;i<images.size();i++)
				cout<<"#"<<i<<" "<<images[i].filename<<" count "<<images[i].count<<endl;

			cout<<"Computing output with: "<<images[0].filename<<endl;

			CImg<double> input_image(images[0].filename.c_str());
			CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
			descriptors = Sift::compute_sift(gray);
			int height = input_image.height();
			if(query_image.height()>height)
				height = query_image.height();
			CImg<double> output(query_image.width()+input_image.width(),height,1,3,0);
			cout<<"_---___---_"<<output.width()<<" "<<output.height()<<endl;
			cout<<"_---___---_"<<query_image.width()+input_image.width()<<" "<<height<<endl;
			int max = (descriptors.size()>query_descriptors.size())?descriptors.size():query_descriptors.size();
			int a,b=0;
			input_image.normalize(0,255);
			
			query_image.normalize(0,255);
			for(int i=0; i<query_image.width();i++)
				for(int j=0;j<query_image.height();j++)
				{
					output(i,j,0,0)=query_image(i,j,0,0);
					output(i,j,0,1)=query_image(i,j,0,1);
					output(i,j,0,2)=query_image(i,j,0,2);
				}
			for(int i=query_image.width(), x=0; i<output.width();x++,i++)
				for(int j=0;j<input_image.height();j++)
				{
					output(i,j,0,0)=input_image(x,j,0,0);
					output(i,j,0,1)=input_image(x,j,0,1);
					output(i,j,0,2)=input_image(x,j,0,2);
				}
			const double color[] ={255.0,128.0,64.0};
			for(int i=0; i<final.size(); i++)
			{
				//	cout << "Descriptor #" << i << ": x=" << descriptors[i].col << " y=" << descriptors[i].row << " descriptor=(";
				/*
				   for(int l=0; l<128; l++)
				   cout << descriptors[i].descriptor[l] << "," ;
				   cout << ")" << endl;
				   */
				a = final[i].d1;
				b = final[i].d2;
				output.draw_line(query_descriptors[a].col,query_descriptors[a].row,descriptors[b].col+query_image.width(),descriptors[b].row,color);
				for(int j=0; j<5; j++)
					for(int k=0; k<5; k++)
						if(j==2 || k==2)
							for(int p=0; p<3; p++)
							{
								//cout<<query_descriptors[i].col+k<<" "<<query_descriptors[i].row+j<<endl;
								output((int)query_descriptors[a].col+k,(int) query_descriptors[a].row+j, 0, p)=0;
								//cout<<descriptors[i].col+k+query_image.width()<<" "<<descriptors[i].row+j<<endl;
								output(descriptors[b].col+k+query_image.width(), descriptors[b].row+j, 0, p)=0;
							}

			}
			output.get_normalize(0,255).save("sift.png");


		}
		else if(part == "part4")
		{
			if(argc < 6)
			{
				cout << "Insufficent number of arguments; correct usage:" << endl;
				cout << "    p4 part_id output img1 img2 ratio" << endl;
				return -1;
			}
			// do something here!
			string inputFile1 = argv[3];
			string inputFile2 = argv[4];
			string outputFile = argv[2];
			float ratio = atof(argv[5]);
			CImg<double> input1(inputFile1.c_str());
			CImg<double> input2(inputFile2.c_str());
			int width = input1.width()+input2.width();
			int height= input1.height()+input2.height();
			cout<<"stitching "<<inputFile1<<" ("<<input1.width()<<"x"<<input1.height()<<") with "<<inputFile2<<" ("<<input2.width()<<"x"<<input2.height()<<")"<<endl;
			CImg<double> output(width,height,1,3,-1);
			cout<<" output width "<< width<<" output height "<<height<<endl;
			cout<<"computing descriptors"<<endl;
			vector<SiftDescriptor> descriptors1 = vector<SiftDescriptor>();
			vector<SiftDescriptor> descriptors2 = vector<SiftDescriptor>();
			vector<desc_pairs> temp = vector<desc_pairs>();
			
			CImg<double> gray; // convert image to grayscale
			gray = input1.get_RGBtoHSI().get_channel(2);
			descriptors1 = Sift::compute_sift(gray);
			
			gray = input2.get_RGBtoHSI().get_channel(2);
			descriptors2 = Sift::compute_sift(gray);
			compareDescriptors2(descriptors1,descriptors2,ratio,temp);
			cout<<"img1 descriptors"<<descriptors1.size()<<" img2 descriptors: "<<descriptors2.size()<<endl;
			cout<<"matches "<<temp.size()<<endl;
			float vec_x =0,vec_y=0;
			getTranslation2(temp, vec_x, vec_y);

			cout<< "Vector x: "<<vec_x<<" Vector y: "<<vec_y<<endl;
			input1.normalize(0,255); input2.normalize(0,255);

			for(int i=0;i<input1.width();i++)
				for(int j=0;j<input1.height();j++)
				{
					output(i+output.width()/2-input1.width()/2,j+output.height()/2-input1.height()/2,0,0) = input1(i,j,0,0);
					output(i+output.width()/2-input1.width()/2,j+output.height()/2-input1.height()/2,0,1) = input1(i,j,0,1);
					output(i+output.width()/2-input1.width()/2,j+output.height()/2-input1.height()/2,0,2) = input1(i,j,0,2);
				}

int take_mean =0;
			for(int i=0;i<input2.width();i++)
				for(int j=0;j<input2.height();j++)
				{
					take_mean=1;
					if(output(i+output.width()/2-input2.width()/2-ceil(vec_x),j+output.height()/2-input2.height()/2-ceil(vec_y),0,0)==-1 || output(i+output.width()/2-input2.width()/2-ceil(vec_x),j+output.height()/2-input2.height()/2-ceil(vec_y),0,1) ==-1 || output(i+output.width()/2-input2.width()/2-ceil(vec_x),j+output.height()/2-input2.height()/2-ceil(vec_y),0,2)==-1)
					take_mean=0;
					
					output(i+output.width()/2-input2.width()/2-ceil(vec_x),j+output.height()/2-input2.height()/2-ceil(vec_y),0,0) += input2(i,j,0,0)+1;

					output(i+output.width()/2-input2.width()/2-ceil(vec_x),j+output.height()/2-input2.height()/2-ceil(vec_y),0,1) += input2(i,j,0,1)+1;
					output(i+output.width()/2-input2.width()/2-ceil(vec_x),j+output.height()/2-input2.height()/2-ceil(vec_y),0,2) += input2(i,j,0,2)+1;
					if(take_mean==1)
					{

						output(i+output.width()/2-input2.width()/2-ceil(vec_x),j+output.height()/2-input2.height()/2-ceil(vec_y),0,0) /=2; 
						output(i+output.width()/2-input2.width()/2-ceil(vec_x),j+output.height()/2-input2.height()/2-ceil(vec_y),0,1) /=2; 
						output(i+output.width()/2-input2.width()/2-ceil(vec_x),j+output.height()/2-input2.height()/2-ceil(vec_y),0,2) /=2; 
					}
				}
			output.save(outputFile.c_str());
			// feel free to add more conditions for other parts (e.g. more specific)
			//  parts, for debugging, etc.
		}
	} 
	catch(const string &err) {
		cerr << "Error: " << err << endl;
	}
}








