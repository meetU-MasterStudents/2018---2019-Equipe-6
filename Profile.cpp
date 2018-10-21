#include "Profile.h"

using namespace std;

Profile::Profile(char* querySequence, int blusomNumber)
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
	
	
	_profileName = "profile1";
	//choose an appropriate BLOSUM matrix based on the given number!
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
	
	//Memory deallocation for the multi-aligned sequences (output of PSI-Blast)
	for (int i = 0; i < _mAlignedSequences->familySize; i++)
  	{
  		delete[] _mAlignedSequences->maSequences[i];
  	}
	delete[] _mAlignedSequences->maSequences;
	
	//Delete the MSA structure
	delete _mAlignedSequences;
}

//Calculates the frequency matrix
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
	return 0;
}


int Profile::PSSMCalculator()
{
	CalculateFrequencyMatrix();
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

int Profile::CallBLAST()
{
	cout<<"Executing PSI-Blast ...";
	if (system(NULL))
	{
		puts ("Ok");
	}
	else
	{
		exit (EXIT_FAILURE);
	}
	//Nucleotide sequence
	system ("blastdbcmd -db refseq_rna.00 -entry nm_000122 -out test_query.fa");
	system("blastn -query test_query.fa -db refseq_rna.00 -task blastn -dust no -outfmt \"7 qseqid sseqid evalue bitscore\" -max_target_seqs 2");
	
	//Amino acide sequence (protein)
	//number of sequences:
	system("grep -c '^>' ..//blastdb//cow.1.protein.faa");
	system("head -6 ..//blastdb//cow.1.protein.faa > ..//blastdb//cow.small.faa");
	system("makeblastdb -in ..//blastdb//human.1.protein.faa -dbtype prot");
	//system("blastp -query cow.small.faa -db human.1.protein.faa");
	system("blastp -query ..//blastdb//cow.small.faa -db human.1.protein.faa -out cow_vs_human_blast_results.txt");
	//system("less cow_vs_human_blast_results.txt");
	
	return 0;
}

void Profile::printname()
{
	cout<<_profileName<<"\n";
}