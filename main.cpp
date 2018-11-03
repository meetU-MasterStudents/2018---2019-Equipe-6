#include <dirent.h>
#include "Profile.h"

using namespace std;

int ProcessHomstrad(string homstradPath)
{
	Profile* _profile;
    DIR* mainDirHandler = NULL;
	DIR* subDirHandler = NULL;
    struct dirent* familyFolder = NULL; 
	struct dirent* subFolder = NULL;
	ifstream ifHandler;
	string line;
	char path[512];
    mainDirHandler = opendir(homstradPath.c_str());
	
    if (mainDirHandler == NULL)
	{
		cout<<"Problem in reading main HOMASTRAD directory."<<endl;
        exit(1);
	}
    cout<<"HOMASTRAD dataset is open."<<endl;
	
    while ((familyFolder = readdir(mainDirHandler)) != NULL) 
	{
		if(strcmp(familyFolder->d_name,".") == 0 || strcmp(familyFolder->d_name,"..") == 0)
		{
			continue;
		}
		cout<<"Processing family "<<familyFolder->d_name<<" ... "<<endl;
		
		strcpy(path,homstradPath.c_str());
		strcat(path,familyFolder->d_name);
		subDirHandler = opendir(path);

		strcpy(path,homstradPath.c_str());
		strcat(path,"//");
		strcat(path,familyFolder->d_name);
		strcat(path,"//");
		strcat(path,familyFolder->d_name);
		strcat(path,".map");
		
		_profile = new Profile(familyFolder->d_name, path);
		_profile->PSSMCalculator();	
		/*
		ifHandler.open((path));
		if (ifHandler.is_open())
		{
			while (getline (ifHandler, line))
			{
				cout<<line<<endl;
			}
		}
		ifHandler.close();
		*/
		delete _profile;
		closedir(subDirHandler);
	}	
	closedir(mainDirHandler);
    return 0;
}


int ProcessQuery(string queryPath)
{
	Profile* _profile = new Profile(queryPath);
	_profile->CallBLAST();
	_profile->CallMUSCLE();
	_profile->PSSMCalculator();	
	_profile->ProfileName();
	delete _profile;
	return 0;
}

int main()
{
	//ProcessQuery("query.fasta");
	ProcessHomstrad("..//HOMSTRAD");
	return 0;
}
