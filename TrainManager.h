#include<iostream>

class TrainManager
{

//Function Members
public:
  //Constructor
  TrainManager(char*);
  //Destructor
  ~TrainManager();
  //Calculate PSSM
  void PSSMCalculator(char*);
  void printname();

//Data members
private:
  char* _querySequence;
  float** _PSSM;
  const int _nRows = 20;
  const int _nCols = 20;
  char* _name;
};

//Operator overloading for matrix operations!
