#include<iostream>

#include "TrainManager.h"

using namespace std;

TrainManager::TrainManager(char* querySequence)
{
	_querySequence = querySequence;
	_PSSM = new float*[_nRows];
	for(int i = 0; i < _nRows; i++)
	{
    		_PSSM[i] = new float[_nCols];
	}
	_name = "seq1";
}

TrainManager::~TrainManager()
{
	//dont delete query sequence here because it is initialized somewhere else
	for (int i = 0; i < _nRows; i++)
  	{
  		delete[] _PSSM[i];
  	}
	delete[] _PSSM;
}

void TrainManager::PSSMCalculator(char* querySequence)
{
	for(int i = 0; i < _nRows; i++)
	{
		for(int j = 0; j < _nCols; j++)
		{
			_PSSM[i][j] = 0;
		}
	}
}

void TrainManager::printname()
{
	cout<<_name<<"\n";
}
