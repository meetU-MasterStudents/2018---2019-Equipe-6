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



int main(void)
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
	  myReadFile >> output;
	  std::cout<<output;

	}
      }
      }



    if (closedir(hom) == -1) /* S'il y a eu un souci avec la fermeture */
        exit(-1);
    puts("Le dossier HOMASTRAD a ete ferme avec succes");

    return 0;
}


