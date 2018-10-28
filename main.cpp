#include<iostream>
#include "Profile.h"

using namespace std;

int main()
{
	Profile* _profile = new Profile("query.fasta");
	_profile->CallBLAST();
	_profile->CallMUSCLE();
	_profile->PSSMCalculator();
	_profile->ProfileName();
	delete _profile;
	return 0;
}
