#!/usr/bin/env python
# I create cubature using Gauss-Legendre or Gauss-Hermite quadrature.
# The first input chooses which.
# The second input is the number of dimensions over which to integrate.
# The next inputs are for GL the bounds of the UQ space supports,
# for GH they are the average and standard deviation of the Gaussian
# in each domain.
# The final line is the order of quadrature to use in each dimension.
# Output is ordered weights followed by ordered abscissas.
import numpy as np

def gen_Gauss_Legendre(a,b,Q):
    GLx,GLw = np.polynomial.legendre.leggauss(Q) #get Gauss-Legendre wgts and nodes
    GLw = GLw * (b-a) / 2.0                   #normalize weights
    GLx = GLx * (b-a) / 2.0 + (b+a) / 2.0     #map nodes
    return GLw,GLx

def gen_Gauss_Hermite(mu,sig,Q):
    GHx,GHw = np.polynomial.hermite.hermgauss(Q) #get Gauss-Hermite wgts and nodes
    GHw = GHw / np.sqrt(np.pi)                #normalize weights
    GHx = GHx *np.sqrt(2.0)*sig + mu          #map nodes
    return GHw,GHx

def increment_SCpt(SCpt,Q,n):
    #increment to next point in scheme
    SCpt[n] += 1
    if SCpt[n]==Q[n]:
        SCpt[n] = 0
        increment_SCpt(SCpt,Q,n+1)

def gen_wgts_nodes(Q,abms,quadtype):
    #generate cubature weights and nodes (as lists) in corresponding order
    #solve wgts and nodes for quadrature in each dimension j
    for j in range(0,len(Q)):
        if j==0:                      #initialize/append next data slots
            wgts  = [[0] * Q[j]]
            nodes = [[0] * Q[j]]
        else:
            wgts.append( [0]*Q[j])
            nodes.append([0]*Q[j])
        if quadtype=='GL':            #get new data
            twgts,tnodes = gen_Gauss_Legendre(abms[j][0],abms[j][1],Q[j])
        elif quadtype=='GH':
            twgts,tnodes = gen_Gauss_Hermite(abms[j][0],abms[j][1],Q[j])
        for k in range(0,len(twgts)): #store new data
            wgts[j][k] = twgts[k]
            nodes[j][k]= tnodes[k]
    
    #use quadrature wgts and nodes to create cubature wgts and nodes
    for i in range(0,np.prod(Q)):
        if i==0:                      #initialize/append next data slots
            SCpt =  [0]*len(Q)
            cwgts= [0]
            cnodes=[[0]*len(Q)]
        else:
            cwgts.append(0)
            cnodes.append([0]*len(Q))
        cwgts[i] = 1.0
        for j in range(0,len(Q)):     #store new data
            cwgts[i]    = cwgts[i]*wgts[j][SCpt[j]]
            cnodes[i][j]= nodes[j][SCpt[j]]
        if i == np.prod(Q)-1:         #exit when data for last point stored
            return cwgts,cnodes
        increment_SCpt(SCpt,Q,0)      #get indices for new point

def order_wgts_nodes(cwgts,cnodes):
    #orders weights and nodes in decreasing order of weight
    for i in range(0,len(cwgts)-1):     #make largest remaining wgt first
        maxwgt = 0.0
        maxii  = i
        for ii in range(i,len(cwgts)):  #find largest weight in remainder
            if cwgts[ii]>maxwgt:
                maxwgt = cwgts[ii]
                maxii  = ii
        if not i==maxii:                #put largest at beginning or rem.
            twgts        = cwgts[i]
            cwgts[i]     = cwgts[maxii]
            cwgts[maxii] = twgts
            tnodes       = cnodes[i]
            cnodes[i]    = cnodes[maxii]
            cnodes[maxii]= tnodes

#read and format input
text_file1 = open("auxiliary/cubature/Gen_Ord_Cubature.inp.txt")
quadtype = text_file1.readline()
quadtype = str(quadtype).strip('\n')
quadtype = str(quadtype).strip('\r')

numD = text_file1.readline()
numD = str(numD).strip('\n')
numD = int(numD)

abms = []
for k in range(0,numD):
    text = text_file1.readline()
    text = text.split(' ')
    abms.append([float(i) for i in text])

Q = text_file1.readline()
Q    = Q.split(' ')
Q    = [int(i) for i in Q]

text_file1.close()

cwgts,cnodes = gen_wgts_nodes(Q,abms,quadtype) #generate wgts and nodes
order_wgts_nodes(cwgts,cnodes)                    #put in descending wgt order

#output
text_file = open("auxiliary/cubature/Gen_Ord_Cubature.w.x.txt","w")
for w in cwgts:
    text_file.write("%15.13f\n"%w)
text_file.write("\n")
for pt in cnodes:
    for d in pt:
        text_file.write("%15.13f  "%d)
    text_file.write("\n")
text_file.close()
