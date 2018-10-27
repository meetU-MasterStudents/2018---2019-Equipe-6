#include "Profile.h"

using namespace std;

Profile::Profile(char* querySequence)
{
	_querySequence = querySequence;
	
	_mAlignedSequences = new MultiAlignedSequences();
	
	//Memory allocation for the substitution matrix (and zero initializing)
	_PSSM = new float*[_nAmAc];
	for(int i = 0; i < _nAmAc; i++)
	{
    		_PSSM[i] = new float[_nAmAc]();
	}	
	
	
	//Memory allocation for the multi-aligned sequences (output of PSI-Blast)
	_mAlignedSequences->familySize = 4; //just for test
	_mAlignedSequences->seqLength = 5;  //just for test
	_mAlignedSequences->maSequences = new char*[_mAlignedSequences->familySize];
	for(int i = 0; i < _mAlignedSequences->familySize; i++)
	{
    		_mAlignedSequences->maSequences[i] = new char[_mAlignedSequences->seqLength];
	}
	
	//just for test
	_mAlignedSequences->maSequences[0][0] = 'A';
	_mAlignedSequences->maSequences[0][1] = 'B';
	_mAlignedSequences->maSequences[0][2] = 'C';
	_mAlignedSequences->maSequences[0][3] = 'C';
	_mAlignedSequences->maSequences[0][4] = 'B';
	
	_mAlignedSequences->maSequences[1][0] = 'A';
	_mAlignedSequences->maSequences[1][1] = 'C';
	_mAlignedSequences->maSequences[1][2] = 'A';
	_mAlignedSequences->maSequences[1][3] = 'B';
	_mAlignedSequences->maSequences[1][4] = 'A';
	
	_mAlignedSequences->maSequences[2][0] = 'C';
	_mAlignedSequences->maSequences[2][1] = 'B';
	_mAlignedSequences->maSequences[2][2] = 'C';
	_mAlignedSequences->maSequences[2][3] = 'C';
	_mAlignedSequences->maSequences[2][4] = 'B';
	
	_mAlignedSequences->maSequences[3][0] = 'A';
	_mAlignedSequences->maSequences[3][1] = 'C';
	_mAlignedSequences->maSequences[3][2] = 'B';
	_mAlignedSequences->maSequences[3][3] = 'B';
	_mAlignedSequences->maSequences[3][4] = 'B';
	
	
	//Memory allocation for the frequency matrix (and zero initializing)
	for(int i = 0; i < _nAmAc; i++)
	{
    		_frequencyMatrix[AmAc[i]] = new int[_mAlignedSequences->seqLength]();
	}
	
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
	
	
	
	_profileName = "profile1";
}

Profile::~Profile()
{
	//dont delete query sequence here because it is initialized somewhere else
	
	//Memory deallocation for the substitution matrix
	for (int i = 0; i < _nAmAc; i++)
  	{
  		delete[] _PSSM[i];
  	}
	delete[] _PSSM;
	
	
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
	
	//Memory deallocation for the multi-aligned sequences (output of PSI-Blast)
	for (int i = 0; i < _mAlignedSequences->familySize; i++)
  	{
  		delete[] _mAlignedSequences->maSequences[i];
  	}
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
	CalculateFrequencyMatrix();
	CalculatePairFrequencyMatrix();
	DisplayFrequencyMatrix();
	
	for(int i = 0; i < _nAmAc; i++)
	{
		for(int j = 0; j < _nAmAc; j++)
		{
			_PSSM[i][j] = 0;
		}
	}
	return 0;
}


//This version is not final and needs modifications!
int Profile::CallBLAST()
{
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
	
	//Amino acide sequence (protein)
	//number of sequences:
	system("grep -c '^>' ..//blastdb//cow.1.protein.faa");
	system("head -6 ..//blastdb//cow.1.protein.faa > ..//blastdb//cow.small.faa");
	
	system("psiblast -query ..//blastdb//cow.small.faa -db pdb -out cow_blast_results.txt -evalue 1e-5 -outfmt 6 -remote "); //Please check this website: https://blast.ncbi.nlm.nih.gov/Blast.cgi
	
	vector<string> hits;
	//Extracting best hits (protein pdb name!)
	ReadHits("cow_blast_results.txt",&hits);
	
	string command;
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
		
	WriteHits(BlastOutputName,"..//blastdb//cow.small.faa",hits);
	
	hits.clear();
	return 0;
}

//Execute MUSCLE algorithm to get multiple sequence alignement
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

//Create unaligned fasta file
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

void Profile::printname()
{
	cout<<_profileName<<"\n";
}