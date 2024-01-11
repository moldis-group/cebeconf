from cebeconf import calc_be

calc_be('glucose_chain_UFF.xyz')


# Default maximum neightbors for large system
Max_neigh={"MaxN_C":16, "MaxN_N":12, "MaxN_O":8, "MaxN_F":4}
calc_be('glucose_chain_UFF.xyz', **Max_neigh)

for N_F in range(1,5):
    Max_neigh={"MaxN_C":16, "MaxN_N":12, "MaxN_O":8, "MaxN_F":N_F}
    print(Max_neigh)
    calc_be('glucose_chain_UFF.xyz', **Max_neigh)
