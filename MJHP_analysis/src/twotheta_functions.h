#ifndef twotheta_functions_H
#define twotheta_functions_H  

  void calculate_powder_pattern(FILE *flog, TwoTheta *TTH, BinaryGrid *BIN, NumberGrid *GRD, UnitCell *UC, Symmetry *SYM); 
  
  void fold_reflections_toBZ(FILE *flog, TwoTheta *TTH); 
  
  void symmetry_folded_reflections(FILE *flog, TwoTheta *TTH, Symmetry *SYM); 

  void print_reflections(char filename[200], TwoTheta *TTH);

  void read_reflections(char filename[200], TwoTheta * TTH);

  void concatinate_twotheta_potential(FILE * flog, EnergyContribution * ECON, EnergyStep * ESTP, TwoTheta *TTH); 

#endif
