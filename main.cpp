#include<iostream>
#include "Profile.h"

using namespace std;

int main()
{
	Profile* _profile = new Profile("abcd");
	_profile->CallBLAST();
	_profile->CallMUSCLE();
	_profile->PSSMCalculator();
	_profile->printname();
	delete _profile;
	return 0;
}
