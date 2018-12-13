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
#include <math.h>

using namespace std; 

//char AmAc[] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-'};
const char AmAc[] = "ARNDCQEGHILKMFPSTWYV-\0";

const string HomstradResultPath = "HomstradResults//";
const string QueryResultPath = "QueryResults//";

struct MultiAlignedSequences
{
  int familySize; 
  int seqLength;
  string* maSequences; 
};

typedef map<char,int*> frequencyMatrix;
typedef map<char,float*> nFrequencyMatrix;
typedef tuple<char, char> residuesPair;
typedef map<residuesPair,int*> pairFrequencyMatrix;
typedef tuple<string,int,int> hitInformation;

class Profile
{
//Function Members
public:
  //Constructors
  Profile(string,string,string,string,bool,bool,bool);
  Profile(string,string);
  //Destructor
  ~Profile();
  //Calculate PSSM
  int PSSMCalculator();
  int CallBLAST();
  int CallMUSCLE();
  int ProfileName();
  int DisplayFrequencyMatrix();
  int WriteFrequencyMatrix();
  int ReadMSA();

//Private function members
private:
	unsigned long long CalculateFactorial(int);
	float CalculatePermutation(bool,int,int);
	float CalculateCombination(bool,int,int);
	int CalculateFrequencyMatrix();
	int CalculatePairFrequencyMatrix();
	int ReadHits(string,vector<hitInformation>* );
	int WriteHits(string,string,vector<hitInformation>);
  int DownloadFromPDB(vector<hitInformation>);
  int DownloadFromUniProt(vector<hitInformation>);
  int LocalDownloadFromBlastDB(vector<hitInformation>);
	int DisplayFasta(string);

  
//Data members
private:
  string _queryFile;
  string _queryName;
  //float** _PSSM;
  //const int _nAmAc = 21;
  const int _nAmAc = strlen(AmAc);
  //const int _nCols = 21;
  string _profileName;
  MultiAlignedSequences* _mAlignedSequences;
  
  //Frequency (Occurrences) matrix: PSSM (first step to create substitution matrix)
  frequencyMatrix _frequencyMatrix;

  //Normalized frequency matrix
  nFrequencyMatrix _nFrequencyMatrix;

  //Log-odd ratio matrix
  nFrequencyMatrix _lorMatrix;
  
  //Pair frequency matrix (second step to create substitution matrix)
  pairFrequencyMatrix _pairFrequencyMatrix;
  
  //Pair residues
  int _numPairResidues;
  
  string _blastOutputName;
  string _fastasAcquiredPath;
  string _muscleOutputName;
  string _pssmOutputName;

  string _eValue = "1e-5";
  string _dataBase = "swissprot";
  bool _applyWeights = false;
  bool _multiThreadBlast = false;
  bool _remoteProcess = false;
  
  bool homstradIrregularity = false;
  int _sp = 0;
};

//Operator overloading for matrix operations!
//A class for errors and error messages!
//A class for events and event messages!
//Multithread: each query one thread
//Use exit(-1) for std errors
//Create a log report! (for example for bad characters found outside AmAc)