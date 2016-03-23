#!/usr/bin/env python
# I read the covariance type (1 or 2), the slab thickness, process average,
# process variance, process correlation length, the number of eigenmodes, and
# the number of points for a numerical eigenvector solve.  I solve the
# Fredholm integral equation which yields eigenvalues and eigenvectors
# associated with a certain covariance.  I am currentely enabled
# to do this for an exponential covariance or the covariance required
# for the KL expansion when performing a lognormal transformation
# of an exponential covariance.
import numpy as np
import time

def expcov(x1,x2,lamc):                    #exponential covariance
    return np.exp(-abs(x1-x2)/lamc)

def Gausscov(x1,x2,procave,procvar,lamc):  #Gaussian cov for LN translation
    num = np.log(expcov(x1,x2,lamc)*procvar/procave**2+1.0)
    den = np.log(procvar/procave**2+1.0)
    return num/den

#read input
covtypecode,s,procave,procvar,lamc,maxnumeigs,numNystrom = np.loadtxt('auxiliary/Nystrom/Nystrominp.txt')
covtypecode = int(covtypecode)
assert covtypecode==1 or covtypecode==2
if covtypecode==1:
    covtype = 'numerical-exp'
if covtypecode==2:
    covtype = 'numerical-Gauss'
maxnumeigs = int(maxnumeigs)
numNystrom = int(numNystrom)
#print 'covtypecode:',covtypecode
#print 's      :',s
#print 'procave:',procave
#print 'procvar:',procvar
#print 'lamc   :',lamc
#print 'maxnumeigs:',maxnumeigs
#print 'numNystrom:',numNystrom


#Numerical solve of Fredholm integral equation using Nystrom method
print 'starting Fredholm solve'
tstart = time.time()
#setup autocov matrix (K)
K = [[0.0]*numNystrom for i in range(numNystrom)]
step = s/(numNystrom+1)
Nyxvec = np.linspace(step/2.0,s-step/2.0,numNystrom)
#print 'covtype:',covtype
for i,x1 in enumerate(Nyxvec):
    for j,x2 in enumerate(Nyxvec):
        if covtype=='numerical-exp':
            K[i][j] = expcov(x1,x2,lamc)
        if covtype=='numerical-Gauss':
            K[i][j] = Gausscov(x1,x2,procave,procvar,lamc)
#setup weight matrices
W = np.zeros_like(K)
for i in range(0,numNystrom): #Simpson's rule quadrature
    W[i][i] = 1.0/numNystrom*s
Whalf = np.sqrt(W)
Whalfinv = np.linalg.inv(Whalf)
W = np.diag(W)
#setup "Nystrom matrix" (B)
B = np.dot(Whalf,np.dot(K,Whalf))
#Nystrom eigenvalue solve
numereigs,eigvecs_star = np.linalg.eig(B)
#keep eigs that will use in KL truncation (rest used in spatial eigfuncs)
eigs = numereigs[:maxnumeigs]
#eigenvectors from Nystrom eigenvectors
eigvecs = np.zeros_like(eigvecs_star)    #initialize eigenvectors
eigvecs = np.dot(Whalfinv,eigvecs_star)  #matrix-wise translation
for i in range(0,numNystrom):            #anchor all correctly (+/-)
    if eigvecs[0,i]<0:
        eigvecs[:,i] = -eigvecs[:,i]
print 'finished Fredholm solve, took',time.time()-tstart,' seconds'

#print 'eigs:',eigs
#print 'eigvecs:',eigvecs

#print output
#print len(eigs),len(eigvecs[0][:])
f1 = open('auxiliary/Nystrom/Nystromout.txt','w+')
print >>f1, 'eigenvalues:'
for i in range(0,maxnumeigs):
    print >>f1, eigs[i]
for i in range(0,maxnumeigs):
    print >>f1, ' '
    print >>f1, 'eigenvector',i
    for ix in range(0,numNystrom):
        print >>f1, eigvecs[ix][i]
f1.close()
