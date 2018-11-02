#include "Profile.h"
#include <stdlib.h>
#include <stdio.h>
/* stdlib pour exit(), et stdio pour puts() */
#include <string.h>
#include <dirent.h>
/* Pour l'utilisation des dossiers */


#include <iostream>
#ifndef WIN32
    #include <sys/types.h>
#endif


using namespace std;

Profile::Profile(string queryFile)
{
	_queryFile = queryFile;
	
	_mAlignedSequences = new MultiAlignedSequences();
	_mAlignedSequences->familySize = 0;
	_mAlignedSequences->seqLength = 0;	
	
	_profileName = "profile";
}

Profile::~Profile()
{
	//Memory deallocation for the frequency matrix
	for (int i = 0; i < _nAmAc; i++)
  	{
  		delete[] _frequencyMatrix[AmAc[i]];
  	}
		
	//Memory deallocation for the pair frequency matrix
	for (pairFrequencyMatrix::iterator it = _pairFrequencyMatrix.begin(); it!=_pairFrequencyMatrix.end(); ++it)
	{
		delete[] _pairFrequencyMatrix[it->first];
	}
	
	//Memory deallocation for the multi-sequence alignement (output of MUSCLE)
	delete[] _mAlignedSequences->maSequences;
	
	//Delete the MSA structure
	delete _mAlignedSequences;
}

//Calculate factorial
unsigned long long Profile::CalculateFactorial(int n)
{
	if(n > 1)
	{
		return n * CalculateFactorial(n - 1);
	}
    else
	{
        return 1;
	}
}

//Calculate permutations with or without repetition
float Profile::CalculatePermutation(bool repetition, int n, int k)
{
	if(repetition == true)
	{
		return pow (n, k);
	}
    else
	{
        return CalculateFactorial(n) / CalculateFactorial(n-k);
	}
}

//Calculate combinations with or without repetition
float Profile::CalculateCombination(bool repetition, int n, int k)
{
	if(repetition == true)
	{
		return CalculateFactorial(n+k-1) / (CalculateFactorial(k) * CalculateFactorial(n-1));
	}
    else
	{
        return CalculateFactorial(n) / (CalculateFactorial(k) * CalculateFactorial(n-k));
	}
}

//Calculate the frequency matrix
int Profile::CalculateFrequencyMatrix()
{
	//Memory allocation for the frequency matrix: PSSM (and zero initializing)
	for(int i = 0; i < _nAmAc; i++)
	{
    		_frequencyMatrix[AmAc[i]] = new int[_mAlignedSequences->seqLength]();
	}
	
	for(int i=0; i < _mAlignedSequences->familySize; i++)
	{
		for(int j=0; j < _mAlignedSequences->seqLength; j++)
		{
			_frequencyMatrix[_mAlignedSequences->maSequences[i][j]][j]++;
		}
	}
	return 0;
}

//Calculate the pair frequency matrix
int Profile::CalculatePairFrequencyMatrix()
{
	residuesPair pair;
	char r1,r2;
	
	//Memory allocation for the pair residues
	_numPairResidues = CalculatePermutation(false,_nAmAc,2);
	
	//Memory allocation for the pair frequency matrix (and zero initializing)
	int count = 0;
	residuesPair pair1,pair2;
	for(int i = 0; i < _nAmAc; i++)
	{
		for(int j = 0; j < _nAmAc; j++)
		{
			get<0>(pair1) = AmAc[i]; get<1>(pair1) = AmAc[j];
			get<0>(pair2) = AmAc[j]; get<1>(pair2) = AmAc[i];
			if(_pairFrequencyMatrix.count(pair1) == 0 && _pairFrequencyMatrix.count(pair2) == 0)
			{	
				_pairFrequencyMatrix[pair1] = new int[_mAlignedSequences->seqLength]();
				count++;
			}
		}
	}
	
	
	for (pairFrequencyMatrix::iterator it = _pairFrequencyMatrix.begin(); it!=_pairFrequencyMatrix.end(); ++it)
	{
		r1 = get<0>(it->first);
		r2 = get<1>(it->first);
		for(int j=0; j < _mAlignedSequences->seqLength; j++)
		{
			if(r1 != r2)
			{
				_pairFrequencyMatrix[it->first][j] = _frequencyMatrix[r1][j] * _frequencyMatrix[r2][j];
			}
			else
			{
				_pairFrequencyMatrix[it->first][j] = CalculateCombination(false,_frequencyMatrix[r1][j],2);
			}
		}
	}
	
	return 0;
}

int Profile::WriteFrequencyMatrix()
{
	ofstream ofHandler;
	ofHandler.open((PSSMOutputName));
	for(int i = 0; i < _nAmAc; i++)
	{
		ofHandler<<AmAc[i]<<": ";
		for(int j=0; j < _mAlignedSequences->seqLength; j++)
		{
			ofHandler<<_frequencyMatrix[AmAc[i]][j]<<"\t";
		}
		ofHandler<<endl;
		ofHandler.flush();
	}
	ofHandler.close();
}

int Profile::DisplayFrequencyMatrix()
{
	for(int i = 0; i < _nAmAc; i++)
	{
		cout<<AmAc[i]<<": ";
		for(int j=0; j < _mAlignedSequences->seqLength; j++)
		{
			cout<<_frequencyMatrix[AmAc[i]][j]<<"\t";
		}
		cout<<endl;
	}
	
	for (pairFrequencyMatrix::iterator it = _pairFrequencyMatrix.begin(); it!=_pairFrequencyMatrix.end(); ++it)
	{
		std::cout << get<0>(it->first) << " " << get<1>(it->first) <<": ";
		for(int j=0; j < _mAlignedSequences->seqLength; j++)
		{
			cout<<it->second[j]<<"\t";
		}
		cout<<endl;
	}
	
	return 0;
}

int Profile::PSSMCalculator()
{
	ReadMSA();
	CalculateFrequencyMatrix();
	//CalculatePairFrequencyMatrix();
	//DisplayFrequencyMatrix();
	WriteFrequencyMatrix();
	return 0;
}

//This version is not final and needs modifications!
int Profile::CallBLAST()
{
	string command;
	cout<<"Executing PSI-Blast ...";
	if (system(NULL))
	{
		puts ("Ok");
		cout<<"Computer is thinking very hard ..."<<endl;
		cout<<"Do you want to hear a joke meanwhile???"<<endl;
		cout<<"No!!??"<<endl;
		cout<<"What did one bioinformatician say to another after a get together?"<<endl;
		cout<<"Dude, that party was a BLAST!"<<endl;
	}
	else
	{
		exit (EXIT_FAILURE);
	}
	
	//Executing blast command according to https://blast.ncbi.nlm.nih.gov/Blast.cgi
	command = "psiblast -query ";
	command += _queryFile;
	command += " -db pdb -out blastOut.txt -evalue 1e-4 -outfmt 6 -remote "; //For multithread output file name should be modified
	system(command.c_str()); 
	
	vector<string> hits;
	//Extracting best hits (protein pdb codes!)
	ReadHits("blastOut.txt",&hits);
	
	//Downloading protein sequences from pdb website in fasta format 
	for(int i=0; i <hits.size(); i++)
	{		
		//curl
		command = "curl -s -o ";
		command += hits[i].substr (0,4); //XXXXXXX fix it!!
		command += " 'https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=";
		command += hits[i].substr (0,4);  //XXXXXXX fix it!! 
		command += "&compressionType=uncompressed' > curlLog";
		//Execute command
		system (command.c_str());
		
		DisplayFasta(hits[i].substr (0,4));
	}
	
	WriteHits(BlastOutputName,_queryFile,hits);
	
	hits.clear();
	return 0;
}

//Executes MUSCLE algorithm to get multiple sequence alignement
int Profile::CallMUSCLE()
{
	string command;
	command = "muscle -in ";
	command += BlastOutputName;
	command += " -out ";
	command += MuscleOutputName;
	//Execute command
	system (command.c_str());
	return 0;
}

//Creates unaligned fasta file
int Profile::WriteHits(string fileName,string queryFile, vector<string> hits)
{
	string line;
	ofstream ofHandler;
	ofHandler.open((fileName));
	
	//Begin with fasta of query
	{
		ifstream ifHandler((queryFile));
		if(ifHandler.good())
		{
			while (getline (ifHandler, line).good())
			{
				ofHandler<<line<<endl;
				ofHandler.flush();
				if (line.find(">") != string::npos) 
				{
					_profileName = line.replace(0,1,"");
				}
			}
		}
	}
	
	for(int i=0; i <hits.size(); i++)
	{
		{
			ifstream ifHandler((hits[i].substr(0,4))); //fixt it!
			if(ifHandler.good())
			{
				while (getline (ifHandler, line).good())
				{
					ofHandler<<line<<endl;
					ofHandler.flush();
				}
			}
		//ifHandler.close();
		//ifHandler.clear();
		}
	}
	ofHandler.close();
	return 0;
}

//Display content of a fasta file
int Profile::DisplayFasta(string fileName)
{
 
    ifstream input(fileName);
    if(!input.good())
	{
        cerr << "Error opening '"<<fileName<<"'. Bailing out." << endl;
        return -1;
    }
 
    string line, name, content;
    while(getline( input, line ).good() )
	{
		cout<<line<<endl;
        if(line.empty() || line[0] == '>' )
		{
            if( !name.empty() )
			{ 
                cout << name << " : " << content << endl;
                name.clear();
            }
            if( !line.empty() )
			{
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() )
		{
            if( line.find(' ') != string::npos )
			{
                name.clear();
                content.clear();
            } 
			else
			{
                content += line;
            }
        }
    }
    if( !name.empty() )
	{
        cout << name << " : " << content << endl;
    }
 
    return 0;
}

//Reading the protein names of the best hits given by Psi-BLAST
int Profile::ReadHits(string fileName,vector<string>* hits)
{
        string line;
	ifstream fHandler;
	fHandler.open((fileName));
	if (fHandler.is_open())
	{
		cout<<"The hits are as follow:"<<endl;
		while (getline (fHandler, line))
		{
			cout<<line<<endl;
			stringstream ss(line);
			if (getline(ss, line, '	'))
			{
				if (getline(ss, line, '	'))
				{
					(*hits).push_back(line);
				}
			}
		}
		for(int i=0; i <(*hits).size(); i++)
		{
			cout << (*hits)[i] <<endl;
		}
	}
	fHandler.close();
	return 0;
}

int Profile::ReadMSA()
{
	bool firstSeq = false;
	string line;
	int i=-1;
	ifstream ifHandler((MuscleOutputName));
	if(ifHandler.good())
	{
		while (getline (ifHandler, line).good())
		{
			if (line.find(">") != string::npos) 
			{
				if(_mAlignedSequences->familySize == 0)
				{
					firstSeq = true;
				}
				else
				{
					firstSeq = false;
				}
				_mAlignedSequences->familySize++;
			}
			else if(firstSeq)
			{
				_mAlignedSequences->seqLength += line.size();
			}
		}
	}
		
	//We have reached the end of the file (ifHandler) and the EOF flag is already set. So we have to return to the beginning
	ifHandler.clear();
	ifHandler.seekg(0, ios::beg);
	
	cout<<"Number of sequences: "<<_mAlignedSequences->familySize<<endl;
	cout<<"Length of each sequence: "<<_mAlignedSequences->seqLength<<endl;
	
	
	//Memory allocation for the multi-sequence alignement (output of MUSCLE)
	_mAlignedSequences->maSequences = new string[_mAlignedSequences->familySize]();
		
	if(ifHandler.good())
	{
		while (getline (ifHandler, line).good())
		{
			if (line.find(">") != string::npos) 
			{
				i++;
			}
			else
			{
				_mAlignedSequences->maSequences[i] += line;
			}
		}
	}
	ifHandler.close();
	return 0;
}



int Profile::ReadHomstrad()
{
    DIR* hom = NULL;
    struct dirent* metafold = NULL; /* Déclaration d'un pointeur vers la structure dirent. */
    char chemin[500] = "../2018---2019-partage/Data/HOMSTRAD";
    hom = opendir(chemin); /* Ouverture du dossier qui contient HOMSTRAD .map */

    if (hom == NULL) /* Si le dossier n'a pas pu être ouvert */
        exit(1); /* (mauvais chemin par exemple) */
    puts("Le dossier HOMSTRAD a ete ouvert avec succes");

    while ((metafold = readdir(hom)) != NULL) 

      {
      printf("Le fichier lu s'appelle '%s'\n", metafold->d_name);
      DIR* fold_hom = NULL;
      struct dirent* map = NULL;
      fold_hom = opendir(strcat(chemin, metafold->d_name));
      fichier_map=readdir(fold_hom)
       ifstream myReadFile;
      myReadFile.open(strcat(metafold->d_name, '.map'));
      char output[100];
      if (myReadFile.is_open()) 
	{
	while (!myReadFile.eof()) 
	  {
          // problem here... 

	  }
      	}
      }


    if (closedir(hom) == -1) /* S'il y a eu un souci avec la fermeture */
        exit(-1);
    puts("Le dossier HOMASTRAD a ete ferme avec succes");

    return 0;
}


int Profile::ProfileName()
{
	cout<<"This is profile for query: "<<_profileName<<endl;
	return 0;
}
