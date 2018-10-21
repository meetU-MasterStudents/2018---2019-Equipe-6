#include<iostream>
#include <stdlib.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <map>

using namespace std; 

//Insert all BLOSUM matrices here as constant general variable 2D arrays! this general variable is Ok for the multithreading!

//char AmAc[] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X'};
const char AmAc[] = "ABC\0"; //just for test


struct MultiAlignedSequences
{
  int familySize; 
  int seqLength;
  char** maSequences; 
};


class Profile
{
//Function Members
public:
  //Constructor
  Profile(char*, int);
  //Destructor
  ~Profile();
  //Calculate PSSM
  int PSSMCalculator();
  int CallBLAST();
  void printname();

//Private function members
private:
//some private functions to help building PSSM
int CalculateFrequencyMatrix();
int DisplayFrequencyMatrix();
  
//Data members
private:
  char* _querySequence;
  float** _PSSM;
  //const int _nAmAc = 21;
  const int _nAmAc = strlen(AmAc);
  //const int _nCols = 21;
  char* _profileName;
  MultiAlignedSequences* _mAlignedSequences;  //given by Psi-BLAST!
  
  //Frequency matrix (first step to create substitution matrix)
  map<char,int*> _frequencyMatrix;
};

//Operator overloading for matrix operations!
//A class for errors and error messages!
//A class for events and event messages!
//Multithread: each query one thread!