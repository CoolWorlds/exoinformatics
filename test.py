import numpy as np
import exoinformatics

# Inputs: radii, arbitrary units, rank-ordered in terms of semi-major axis
# Outputs: W's = occupancy of microstate occuped (S=logW)

radii=np.random.rand(5)
print 'random radii = ',radii
N = len(radii)
    
# Calculate the tally
tally=[0 for i in range(N-1)]      # define an empty list
tallyhistory=[0 for i in range(N)] # why have one empty list when you can have two
tallyhistory[0]=0
for i in range(N-1):
    tally[i] = 2*exoinformatics.heaviside(radii[i+1]-radii[i]) - 1
    tallyhistory[i+1] = int(np.sum(tally))
T = int(np.sum(tally))

# Calculate the integral
integral=[0 for i in range(N-1)]     # define another empty list, they're great
for i in range(N-1):
    integral[i] = 0.5*(tallyhistory[i+1]+tallyhistory[i])
area = np.sum(integral)

# Calculate the change score
delta=[0 for i in range(N-2)]     # define yet another empty list, we need more damn it!
for i in range(N-2):
    delta[i] = 1 - exoinformatics.kronecker(tally[i],tally[i+1])
D = int(np.sum(delta))

# Calculate M0, microstate index
M0 = exoinformatics.M0(N,T)
W0 = exoinformatics.W0(N,M0)

# Find out what k you're in, the sub-microstate index
k = exoinformatics.findk(N,M0,area)

print 'W_T, W_I, W_C = ',W0,exoinformatics.WI(N,M0,k),exoinformatics.WC(N,M0,D)