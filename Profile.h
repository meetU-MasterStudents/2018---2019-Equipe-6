#include<iostream>


//Insert all BLOSUM matrices here as constant general variable 2D arrays! this general variable is Ok for the multithreading!

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
  
//Data members
private:
  char* _querySequence;
  char** _alignedSequences; //given by Psi-BLAST!
  float** _PSSM;
  const int _nRows = 21;
  const int _nCols = 21;
  char* _profileName;
};

//Operator overloading for matrix operations!
//A class for errors and error messages!
//A class for events and event messages!
