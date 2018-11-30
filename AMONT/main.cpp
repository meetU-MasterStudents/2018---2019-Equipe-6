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
		cout<<"Problem in reading main HOMSTRAD directory."<<endl;
        exit(1);
	}
    cout<<"HOMSTRAD dataset is open."<<endl;
	
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


int ProcessQuery(string queryPath,string queryName, string eValue, string dataBase)
{
	Profile* _profile = new Profile(queryPath, queryName, eValue, dataBase);
	_profile->CallBLAST();
	_profile->CallMUSCLE();
	_profile->PSSMCalculator();	
	_profile->ProfileName();
	delete _profile;
	return 0;
}

static void ShowUsage(string exeFile)
{
    cerr << "Usage: " << exeFile << " <option(s)> "
         << "Options:\n"
         << "\t-h,--help\t\tShow this help message\n"
		 << "\t-t,--prochoms\t\t<Homstrad path>\t\tHomstrad dataset path\n"
         << "\t-s,--procquery\t\tProcess a query\n"
		 << "\t-q,--query\t\t<Query path>\t\tQuery file path\n"
		 << "\t-e,--evalue\t\t<e-Value>\t\te-Value for PSI-Blast\n"
		 << "\t-d,--database\t\t<database>\t\tDatabase for PSI-Blast"
         << endl;
}

int main(int argc, char* argv[])
{
	bool HomsOrQuery = false;
	string queryFilePath;
	string eValue;
	string dataBase;
	string homstradDatasetPath;

   // Check the number of parameters
    if (argc < 2) 
	{
		ShowUsage(argv[0]);
        return -1;
    }

	string arg;
    for (int i = 1; i < argc; i++) 
	{
        arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) 
		{
            ShowUsage(argv[0]);
            return 0;
        } 
		else if ((arg == "-p") || (arg == "--procquery")) 
		{
			HomsOrQuery = true;
		}
		else if ((arg == "-t") || (arg == "--prochoms")) 
		{
			HomsOrQuery = false;
            if (i + 1 < argc) 
			{
                homstradDatasetPath = argv[++i];
            } 
			else 
			{
                cerr << "--query option requires one argument." << endl;
                return -1;
            } 
		}
		else if ((arg == "-q") || (arg == "--query")) 
		{
            if (i + 1 < argc) 
			{
                queryFilePath = argv[++i];
            } 
			else 
			{
                cerr << "--query option requires one argument." << endl;
                return -1;
            }  
        }
		else if ((arg == "-e") || (arg == "--evalue")) 
		{
            if (i + 1 < argc) 
			{
                eValue = argv[++i];
            } 
			else 
			{
                cerr << "--evalue option requires one argument." << endl;
                return 1;
            }  
        }
		else if ((arg == "-d") || (arg == "--database")) 
		{
            if (i + 1 < argc) 
			{
                dataBase = argv[++i];
            } 
			else 
			{
                cerr << "--database option requires one argument." << endl;
                return 1;
            }  
        }
    }

	if(HomsOrQuery)
	{
		string command;
		string queryFileName;

		int lastPos = queryFilePath.find_last_of("/");
		queryFileName = queryFilePath.substr(lastPos+1, queryFilePath.length());
		queryFileName.replace(queryFileName.find('.'),6,"");

		command = "rm -rf QueryResults//"+queryFileName+"//";
		system (command.c_str());
		command = "mkdir QueryResults//"+queryFileName+"//";
		system (command.c_str());

		command = "mkdir QueryResults//"+queryFileName+"//Fastas//";
		system (command.c_str());

		ProcessQuery(queryFilePath,queryFileName,eValue,dataBase);
	}
	else
	{
		system ("rm -rf HomstradResults");
		system ("mkdir HomstradResults");
		ProcessHomstrad(homstradDatasetPath);
	}
	return 0;
}
