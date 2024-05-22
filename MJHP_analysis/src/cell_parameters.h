#ifndef cell_parameters_H
#define cell_parameters_H 

  void Determine_CellParameters(FILE * flog, UnitCell* UC, NumberGrid* GRD); 

  void Calculate_FermiRadius(UnitCell* UC, FermiSphere* FS); 

  void Calculate_FermiDegree(UnitCell* UC, FermiSphere* FS); 

  void symmorphic_symmetry(Symmetry* SYM); 

#endif
