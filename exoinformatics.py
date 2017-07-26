import numpy as np
import scipy.special as special
import math

# ==============================================================================
def heaviside(x):
  if x < 0:
      y = 0
  else:
      y = 1
  return y
# ==============================================================================

# ==============================================================================
def kronecker(x,y):
  if x == y:
      z = 1
  else:
      z = 0
  return z
# ==============================================================================

# ==============================================================================
def eulerian(n,q):
  A = 0
  for j in range(q+1):
      A = A + (-1)**j*int(special.binom(n+1,j))*(q-j+1)**n
  return A
# ==============================================================================

# ==============================================================================
def M0(N,T):
  M = int(0.5*(N+T+1))
  return M
# ==============================================================================

# ==============================================================================
def W0(N,M):
  W = eulerian(N,M-1)
  return W
# ==============================================================================

# ==============================================================================
def Omega0(N):
  Omega = N
  return Omega
# ==============================================================================

# ==============================================================================
def omegai(N,M):
  if N>1:
      omega = (M-1)*(N-M)+1
  elif N==1:
      omega = 1
  
  return omega
# ==============================================================================

# ==============================================================================
def omegac(N,M):
  if N>1:
      omega = max(2*min(M-1,N-M),1) - kronecker(0.5*(N-1),M-1.0)
  elif N==1:
      omega = 1
  
  return omega
# ==============================================================================

# ==============================================================================
def OmegaI(N):
  Omega = ( 15*N**2 - 2*N**3 -7*N )/6
  return Omega
# ==============================================================================

# ==============================================================================
def OmegaC(N):
  Omega = int(math.floor(((N-1)**2+3)/2))
  return Omega
# ==============================================================================

# ==============================================================================
def findk(N,M,I):

    length=omegai(N,M)
    absmax=0.5*(N-1)*(N-1) - (np.minimum(M-1,N-M))**2
    if M <= 0.5*N:
        Iscore=[(-absmax+2*k) for k in range(length)]
    else:
        Iscore=[(absmax-2*k) for k in range(length)]
    Iscore=np.sort(Iscore) #sort
    Iscore=np.flipud(Iscore)    #reverse
    Iscore=Iscore.tolist() # back to list
    for k in range(length):
        if abs(Iscore[k]-I) < 0.1:
            bestk = k + 1
    #print 'best k = ',bestk
      
    return bestk
# ==============================================================================

# ==============================================================================
def WI(N,M,k):
  if M == 1 or M == N:
      W = 1
  elif k == 1:
      W = int(special.binom(N-1,M-1))
  elif M == 2 or M == N-1:
      W = int(special.binom(omegai(N,M)+1,k))-1
  elif N == 5:
      if M == 3:
          if k == 3:
              W = 22
          else: #elif k == 2 or 4
              W = 16
  elif N == 6:
      if M == 3 or M == 4:
          if k == 4:
              W = 80
          elif k == 3 or k == 5:
              W = 66
          else: #elif k == 2 or k == 6:
              W = 35
  elif N == 7:
      if M == 4:
          if k == 5 or k == 6:
              W = 494
          elif k == 4 or k == 7:
              W = 382
          elif k == 3 or k == 8:
              W = 222
          else: #elif k == 2 or k == 9:
              W = 90
      elif M == 3 or M == 5:
          if k == 5:
              W = 269
          elif k == 4 or k == 6:
              W = 233
          elif k == 3 or k == 7:
              W = 149
          else: #elif k == 2 or k == 8:
              W = 64
  elif N == 8:
      if M == 4 or M == 5:
          if k == 7:
              W = 2785
          elif k == 6 or k == 8:
              W = 2540
          elif k == 5 or k == 9:
              W = 1918
          elif k == 4 or k == 10:
              W = 1175
          elif k == 3 or k == 11:
              W = 560
          else: #elif k == 2 or k == 12:
              W = 189
      elif M == 3 or M == 6:
          if k == 6:
              W = 855
          elif k == 5 or k == 7:
              W = 765
          elif k == 4 or k == 8:
              W = 540
          elif k == 3 or k == 9:
              W = 288
          else: #elif k == 2 or k == 10:
              W = 105
  elif N == 9:
      if M == 5:
          if k == 9:
              W = 23402
          elif k == 8 or k == 10:
              W = 21916
          elif k == 7 or k == 11:
              W = 17956
          elif k == 6 or k == 12:
              W = 12764
          elif k == 5 or k == 13:
              W = 7754
          elif k == 4 or k == 14:
              W = 3918
          elif k == 3 or k == 15:
              W = 1568
          else: # if k == 2 or k == 16:
              W = 448
      elif M == 4 or M == 6:
          if k == 8 or k == 9:
              W = 13536
          elif k == 7 or k == 10:
              W = 11736
          elif k == 6 or k == 11:
              W = 8767
          elif k == 5 or k == 12:
              W = 5561
          elif k == 4 or k == 13:
              W = 2913
          elif k == 3 or k == 14:
              W = 1198
          else: #elif k == 2 or k == 15:
              W = 350
      elif M == 3 or M == 7:
          if k == 7:
              W = 2632
          elif k == 6 or k == 8:
              W = 2400
          elif k == 5 or k == 9:
              W = 1806
          elif k == 4 or k == 10:
              W = 1091
          elif k == 3 or k == 11:
              W = 503
          else: #elif k == 2 or k == 12:
              W = 160
  elif N == 10:
      if M == 5 or M == 6:
          if k == 11:
              W = 171224
          elif k == 10 or k == 12:
              W = 162826
          elif k == 9 or k == 13:
              W = 139841
          elif k == 8 or k == 14:
              W = 108019
          elif k == 7 or k == 15:
              W = 74485
          elif k == 6 or k == 16:
              W = 45297
          elif k == 5 or k == 17:
              W = 23836
          elif k == 4 or k == 18:
              W = 10521
          elif k == 3 or k == 19:
              W = 3690
          else: # if k == 2 or k == 20:
              W = 924
      elif M == 4 or M == 7:
          if k == 10:
              W = 63770
          elif k == 9 or k == 11:
              W = 60213
          elif k == 8 or k == 12:
              W = 50592
          elif k == 7 or k == 13:
              W = 37597
          elif k == 6 or k == 14:
              W = 24427
          elif k == 5 or k == 15:
              W = 13610
          elif k == 4 or k == 16:
              W = 6299
          elif k == 3 or k == 17:
              W = 2295
          else: #elif k == 2 or k == 18:
              W = 594
      elif M == 3 or M == 8:
          if k == 8:
              W = 7934
          elif k == 7 or k == 9:
              W = 7332
          elif k == 6 or k == 10:
              W = 5756
          elif k == 5 or k == 11:
              W = 3776
          elif k == 4 or k == 12:
              W = 2005
          elif k == 3 or k == 13:
              W = 817
          else: #elif k == 2 or k == 14:
              W = 231
  else:
      W = -1
  return W
# ==============================================================================

# ==============================================================================
def WC(N,M,D):
  if M == 1 or M == N:
      W = 1
  elif D == 1:
      W = 2*int(special.binom(N-1,M-1)) - kronecker(M,1) - kronecker(M,N)
  elif D == 2 and omegac(N,M) == 2:
      W = W0(N,M) - ( 2*int(special.binom(N-1,M-1)) - kronecker(M,1) - kronecker(M,N) )
  elif N == 5:
      if D == 2:
          W = 22
      else: #elif D == 3
          W = 32
  elif N == 6:
      if D == 2:
          W = 71
      elif D == 3:
          W = 150
      else: #elif D == 4
          W = 61
  elif N == 7:
      if D == 2:
          if omegac(N,M) == 4:
              W = 200
          else: #elif omegac(N,M) == 5:
              W = 220
      elif D == 3:
          if omegac(N,M) == 4:
              W = 482
          else: #elif omegac(N,M) == 5:
              W = 888
      elif D == 4:
          if omegac(N,M) == 4:
              W = 479
          else: #elif omegac(N,M) == 5:
              W = 724
      else: #elif D == 5
          W = 544
  elif N == 8:
      if D == 2:
          if omegac(N,M) == 4:
              W = 521
          else: #elif omegac(N,M) == 6:
              W = 629
      elif D == 3:
          if omegac(N,M) == 4:
              W = 1316
          else: #elif omegac(N,M) == 6:
              W = 3472
      elif D == 4:
          if omegac(N,M) == 4:
              W = 2414
          else: #elif omegac(N,M) == 6:
              W = 4897
      elif D == 5:
          W = 5166
      else: #elif D == 6
          W = 1385
  elif N == 9:
      if D == 2:
          if omegac(N,M) == 4:
              W = 1290
          elif omegac(N,M) == 6:
              W = 1734
          else: #elif omegac(N,M) == 7:
              W = 1794
      elif D == 3:
          if omegac(N,M) == 4:
              W = 3292
          elif omegac(N,M) == 6:
              W = 11240
          else: #elif omegac(N,M) == 7:
              W = 16032
      elif D == 4:
          if omegac(N,M) == 4:
              W = 9970
          elif omegac(N,M) == 6:
              W = 25568
          else: #elif omegac(N,M) == 7:
              W = 32250
      elif D == 5:
          if omegac(N,M) == 6:
              W = 30552
          else: #elif omegac(N,M) == 7:
              W = 58860
      elif D == 6:
          if omegac(N,M) == 6:
              W = 19028
          else: #elif omegac(N,M) == 7:
              W = 31242
      else: #elif D == 7
          W = 15872
  elif N == 10:
      if D == 2:
          if omegac(N,M) == 4:
              W = 3083
          elif omegac(N,M) == 6:
              W = 4659
          else: #elif omegac(N,M) == 8:
              W = 4999
      elif D == 3:
          if omegac(N,M) == 4:
              W = 7818
          elif omegac(N,M) == 6:
              W = 32682
          else: #elif omegac(N,M) == 8:
              W = 60030
      elif D == 4:
          if omegac(N,M) == 4:
              W = 32867
          elif omegac(N,M) == 6:
              W = 115053
          else: #elif omegac(N,M) == 8:
              W = 173526
      elif D == 5:
          if omegac(N,M) == 6:
              W = 144972
          else: #elif omegac(N,M) == 8:
              W = 408438
      elif D == 6:
          if omegac(N,M) == 6:
              W = 157658
          else: #elif omegac(N,M) == 8:
              W = 359838
      elif D == 7:
          W = 252750
      else: #elif D == 8
          W = 50521
  else:
      W = -1
  return W
# ==============================================================================

# ==============================================================================
def randomswap(radii):
    
    lastelement = len(radii) - 1
    
    # Choose a random element = swap1
    swap1 = np.random.randint(0,lastelement+1)
    
    # Calculate swap2
    if swap1 == 0:
        swap2 = swap1 + 1
    elif swap1 == lastelement:
        swap2 = swap1 - 1
    else:
        swap2 = swap1 + int(2.0*(np.random.randint(0,2) - 0.5))
    
    # Create the swap array
    swaps = np.array([swap1,swap2])
    #print 'swaps = ',swaps[0],swaps[1]
    
    # Conduct the swap
    newradii=[0 for i in range(len(radii))]
    for i in range(len(newradii)):
        if i == swaps[0]:
            newradii[i] = radii[swaps[1]]
        elif i == swaps[1]:
            newradii[i] = radii[swaps[0]]
        else:
            newradii[i] = radii[i]
    
    return newradii
# ==============================================================================