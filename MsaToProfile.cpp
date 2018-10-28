#include<iostream>
#include "Profile.h"

// pour l'automate
#include <stdio.h>
#include <stdlib.h>



using namespace std;

int main()
{
	//Automate
     
#define S_init 0
#define S_identifiant 1
#define S_commentaire 2
#define S_sequence 3
#define S_end 4
#define S_error 5

	
	char cc;
	int state = S_init;
	FILE* file;
	//string* query = new string();
	string query;
	int compteur=0;
	
	file = fopen("/home/roselyne/Bureau/BIM-info/MEETU/2018---2019-Equipe-6/example_MSA.fasta", "r");

	while (state<S_end) {
	  cc=getc(file);
	  switch(state){

	  case S_init:

	    switch(cc){
	    case '\t':
	    case ' ':
	      break;
	    case '>':
	      state=S_identifiant;
	      break;
	    default:
	      state=S_error; 
	    }
	    break;	    
	    
	  case S_identifiant:
	    switch(cc){
	    case '\t':
	    case ' ':
	      state=S_commentaire;
	      break;
	    case '\n':
	      state=S_sequence;
	      break;
	    case EOF:
	      state=S_error;
	    }
	    break;	    
	    
	  case S_commentaire:
	    switch(cc){
	    case '\n':
	      state=S_sequence;
	      break;
	    case EOF:
	      state=S_error;
	    }
	    break;

	  case S_sequence:
	    switch(cc){
	    case 'A':
	    case 'R':
	    case 'D':
	    case 'N':
	    case 'C':
	    case 'E':
	    case 'Q':
	    case 'G':
	    case 'H':
	    case 'I':
	    case 'L':
	    case 'K':
	    case 'M':
	    case 'F':
	    case 'P':
	    case 'S':
	    case 'T':
	    case 'W':
	    case 'Y':
	    case 'V':
	    case '-':
	    case '*':
	    case 'X':
	      if (compteur==0){
		  query += cc;
		}
	      break;
	    case '\n':
	      break;	  
	    case '>':
	      compteur++;
	      state=S_identifiant;
	      break;
	    case EOF:
	      state=S_end;
	      break;
	    default:
	      state=S_error;
	      }
	    break;

	  }

	}
	fclose(file);
	printf(" %s \n", query.c_str());
	//





	return 0;
}


