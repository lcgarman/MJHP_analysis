#ifndef hkl_functions_H
#define hkl_functions_H  

  typedef struct {
    int H;
    int K;
    int L;
  } Conventional;

  void transform_HKL(MottJonesConditions * MJC); 

  void find_symmetric_hkl(VectorIndices* VECT, Symmetry* SYM, NumberGrid* GRD); 

  void find_MJregion(VectorIndices *VECT, UnitCell *UC); 

  void concatinate_HKL_potential(EnergyContribution * ECON, EnergyStep * ESTP); 

  void integrate_HKL_potential(EnergyContribution * ECON, EnergyStep * ESTP, AtomicVariables * ATM);

  void integrate_Erange_potential(EnergyContribution * ECON, EnergyStep * ESTP, VectorIndices * VECT);

  Conventional HKL_convert_toP(MottJonesConditions * MJC, int H, int K, int L);

  void combine_potential_grid(NumberGrid* GRD, BinaryGrid* LOCAL, BinaryGrid* NONLOC, BinaryGrid* POT);

#endif
