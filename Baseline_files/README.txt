The occupation matrix of a molecule dummy begins with the line,
  Mulliken occupations of core levels of molecule dummy
And ends with the line,
  End of occupation matrix.

The matrix is designated as follows: 
1. The first element of first row and first column is the molecule name
  dummy
2. The first element of each column (except first column) are the baseline KS energies obtained 
  dummy  E_1  . . . E_n
3. The first element of each row (except first row) corresponds to the heavy atom number in the geometry file
  dummy   E_1 . . . E_n
  atom_1
  .
  .
  atom_n
4. Each row are the Mulliken occupations for that atom corresponding to KS energies of 1s core levels:
Suppose for the 1st atom in a system containing 3 CONF atoms, the occupations corresponding to atom number 1 in the geometry file are in the first row:
dummy   E_1  E_2  E_3
atom_1  0.0  0.1  0.9
