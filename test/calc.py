from cebeconf import calc_be

# Run with default settings
# calc_be('glucose_chain_UFF.xyz')

# Default maximum neightbors for large system
# Max_neigh={"MaxN_C":16, "MaxN_N":12, "MaxN_O":8, "MaxN_F":4}

# Change the number of neighbors for F

#for FN in range(1,12): 
#    Max_neigh={"MaxN_C":16, "MaxN_N":12, "MaxN_O":8, "MaxN_F":FN} 
#    print("===",Max_neigh,"===")
#    calc_be('1.xyz', 'direct' **Max_neigh)

#Max_neigh={"MaxN_C":23, "MaxN_N":12, "MaxN_O":8, "MaxN_F":8} 
#calc_be('1.xyz', 'direct', **Max_neigh)
#calc_be('1.xyz', 'delta', **Max_neigh)
#calc_be('bigqm7w_002951.xyz', 'direct', **Max_neigh)
#calc_be('bigqm7w_002951.xyz', 'delta', **Max_neigh)
#calc_be('benzene.xyz', 'direct', **Max_neigh)

#calc_be('benzene.xyz', 'direct')
#calc_be('benzene.xyz', 'delta')


#for CH4.xyz  delta_scf_ML_10.xyz  glucose_chain_UFF.xyz  glucose_ring_UFF.xyz  H2O.xyz  HF.xyz  NH3.qxyz  test.xyz
#calc_be('test.xyz', 'direct', 'ACM')
#calc_be('test.xyz', 'delta', 'ACM')
#for i in "CH4.xyz", "NH3.xyz", "H2O.xyz", "HF.xyz":
#out=calc_be("1H_2pyramidinone_0_0_0.xyz", 'delta', 'AtmEnv')
#out=calc_be("CH4.xyz", 'direct', 'ACM')
#out=calc_be("benzene.xyz", 'delta', 'AtmEnv')
#calc_be(i, 'delta', 'AtmEnv')

#out=calc_be("1H_2pyramidinone_2_1_2.xyz", 'direct', 'ACM')
#print(out)
#out=calc_be("1H_2pyramidinone_2_1_2.xyz", 'direct', 'AtmEnv')
#print(out)
out=calc_be("1.xyz", 'direct', 'ACM')
print(out)
out=calc_be("1.xyz", 'delta', 'ACM')
print(out)
