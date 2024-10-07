#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "structures.h"
#include "globals.h"
#include "rwa_functions.h"
#include "twotheta_functions.h"
#include "allocate_memory.h"
#include "hkl_functions.h"

int main(int argc, char * argv[])
{

  /*declaring filenames*/
  char RFLCfilename [200];
  char MJINfilename [200];
  double min_twotheta; 
  double max_twotheta;
  FileCabinet fcab;
  FermiSphere fsph;
  MottJonesConditions mjc;
  TwoTheta tth;
  Initialize_MottJonesConditions(&mjc);
  Initialize_TwoTheta(&tth);

  if (argc == 2) {
      //copy file string to inupt files
    strcpy(MJINfilename, argv[1]);
    min_twotheta = 0.0;
    max_twotheta =  60.0;
  }
  else if (argc == 4) {
    strcpy(MJINfilename, argv[1]);
    min_twotheta = atof(argv[2]);
    max_twotheta = atof(argv[3]);
  }
  else {
    printf("USAGE: find_reflections <*.mjin> ");
	exit(0);
  }

  /*Read the mj file and store variables*/
  read_mjin_header(MJINfilename, &mjc, &fcab, &fsph);

  /*open log file and print header*/
  strcpy(RFLCfilename, fcab.MJOUTfilename);
  strcat(RFLCfilename, ".rflc");

  /*voidsearch reflections for inside min to max 2theta*/
  search_reflection(RFLCfilename, &mjc, &tth, min_twotheta, max_twotheta);

  /*print reflection indices into mjin file*/
  append_mjin(MJINfilename, &tth);

  return(0);
}
