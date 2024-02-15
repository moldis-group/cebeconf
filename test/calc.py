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


calc_be('test.xyz', 'direct', 'ACM')
calc_be('test.xyz', 'delta', 'ACM')
calc_be('test.xyz', 'direct', 'AtmEnv')
calc_be('test.xyz', 'delta', 'AtmEnv')
