#include<iostream>
#include "Profile.h"

using namespace std;

int main()
{
	Profile* _profile = new Profile("abcd", 61);
	//_profile->CallBLAST();
	_profile->PSSMCalculator();
	_profile->printname();
	delete _profile;
	return 0;
}
