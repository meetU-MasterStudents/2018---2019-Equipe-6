#include<iostream>
#include "Profile.h"

using namespace std;

Profile::Profile(char* querySequence, int blusomNumber)
{
	_querySequence = querySequence;
	_PSSM = new float*[_nRows];
	for(int i = 0; i < _nRows; i++)
	{
    		_PSSM[i] = new float[_nCols];
	}
	_profileName = "profile1";
	
	//choose a appropriate BLOSUM matrix based on the given number!
}

Profile::~Profile()
{
	//dont delete query sequence here because it is initialized somewhere else
	for (int i = 0; i < _nRows; i++)
  	{
  		delete[] _PSSM[i];
  	}
	delete[] _PSSM;
}

int Profile::PSSMCalculator()
{
	for(int i = 0; i < _nRows; i++)
	{
		for(int j = 0; j < _nCols; j++)
		{
			_PSSM[i][j] = 0;
		}
	}
	return 0;
}

int Profile::CallBLAST()
{
	return 0;
}

void Profile::printname()
{
	cout<<_profileName<<"\n";
}
