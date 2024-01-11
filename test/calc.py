from cebeconf import calc_be

# Run with default settings
# calc_be('glucose_chain_UFF.xyz')

# Default maximum neightbors for large system
# Max_neigh={"MaxN_C":16, "MaxN_N":12, "MaxN_O":8, "MaxN_F":4}

# Change the number of neighbors for F
Max_neigh={"MaxN_C":16, "MaxN_N":12, "MaxN_O":8, "MaxN_F":2}
calc_be('glucose_chain_UFF.xyz', **Max_neigh)
