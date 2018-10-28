#include<iostream>
#include <stdlib.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <map>
#include <math.h>
#include <tuple>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std; 

//Insert all BLOSUM matrices here as constant general variable 2D arrays! this general variable is Ok for the multithreading!

//char AmAc[] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X'};
const char AmAc[] = "ABC\0"; //just for test

const string BlastOutputName = "QueryBestHits.fasta";
const string MuscleOutputName = "AlignedQueryBestHits.fasta";

struct MultiAlignedSequences
{
  int familySize; 
  int seqLength;
  char** maSequences; 
};

typedef map<char,int*> frequencyMatrix;
typedef tuple<char, char> residuesPair;
typedef map<residuesPair,int*> pairFrequencyMatrix;

class Profile
{
//Function Members
public:
  //Constructor
  Profile(string);
  //Destructor
  ~Profile();
  //Calculate PSSM
  int PSSMCalculator();
  int CallBLAST();
  int CallMUSCLE();
  int ProfileName();
  int DisplayFrequencyMatrix();

//Private function members
private:
unsigned long long CalculateFactorial(int);
float CalculatePermutation(bool,int,int);
float CalculateCombination(bool,int,int);
int CalculateFrequencyMatrix();
int CalculatePairFrequencyMatrix();
int ReadHits(string,vector<string>* );
int WriteHits(string,string,vector<string>);
int DisplayFasta(string);

  
//Data members
private:
  string _queryFile;
  float** _PSSM;
  //const int _nAmAc = 21;
  const int _nAmAc = strlen(AmAc);
  //const int _nCols = 21;
  string _profileName;
  MultiAlignedSequences* _mAlignedSequences;  //given by Psi-BLAST!
  
  //Frequency matrix (first step to create substitution matrix)
  frequencyMatrix _frequencyMatrix;
  
  //Pair frequency matrix (second step to create substitution matrix)
  pairFrequencyMatrix _pairFrequencyMatrix;
  
  //Pair residues
  int _numPairResidues;
};

//Operator overloading for matrix operations!
//A class for errors and error messages!
//A class for events and event messages!
//Multithread: each query one thread!