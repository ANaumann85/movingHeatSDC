'''
Saves quadrature weights for spectral deferred corrections for different types of nodes for M=2 to M=8 nodes.
Generates a collection of ASCII files each with M+2 x M entries arranged in the following format
S[0,0] .... S[0,M-1]
S[1,0] .... S[1,M-1]
  .
  .
  .
S[M-1,0] .. S[M-1,M-1]
t[0] ...... t[M-1]
q[0] ...... q[M-1]
Here, S[m,n] are the node-to-node quadrature weights, t[n] the nodes and q[n] the weights for quadrature over the full interval.
All values are computed for the unit interval [0,1] but can easily be scaled to an arbitrary interval [a,b] by multiplying the S and q with (b-a)
and t --> a + (b-a)*t
'''

import Collocation
import CollocationClasses
import os

types = [ ["CollGaussLegendre", "legendre"], ["CollGaussLobatto", "lobatto"], ["CollGaussRadau_Right", "radau_right"], ["CollGaussRadau_Left", "radau_left"], ["EquidistantNoLeft", "equi_noleft"] ]

precision = "%19.16f"
separator = "    "
dir = os.path.join(".", "sdc_quad_weights")
if not os.path.exists(dir):
  os.makedirs(dir)

for t in types:
  for M in range(2,9):
    coll = getattr(CollocationClasses, t[0])(M, 0.0, 1.0)
    filename = t[1]+"-M"+str(M)+".dat"
    file = open(os.path.join(dir,filename), 'w')
    Smat = coll.Smat
    Smat = Smat[1:,1:]
    for m in range(M):
      for n in range(M):
        file.write((precision % Smat[m,n]))
        file.write(separator)
      file.write("\n")
    for m in range(M):
      file.write((precision % coll.nodes[m]))
      file.write(separator)
    file.write("\n")
    for m in range(M):
      file.write((precision % coll.weights[m]))
      file.write(separator)