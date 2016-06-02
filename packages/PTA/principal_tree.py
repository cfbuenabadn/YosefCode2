import numpy as np
import scipy.sparse as sparse
from scipy.sparse.csgraph import minimum_spanning_tree
from numpy.matlib import repmat
import matplotlib.pyplot as plt


class principalTree:
    def __init__(self,xs,params,num_m):
        self.x = xs
        self.params = params
        self.k = num_m
        self.MU = None
        self.d = len(xs)
        self.n = len(xs[0]) 
        self.history = History()
        self.model = Model()
        
        if 'MU' in self.params:
            self.MU = self.params['MU']
        else:
            self.MU = np.zeros((self.d,num_m))
            index = np.random.permutation(self.n)
            for i in range(num_m):
                self.MU[:,i]= self.x[:,index[i]]
            
        self.Lambda = self.params['lambda']*self.n
        
        if 'bandwidth' not in self.params:
            distsqX = sqdist(X,X)
            self.sigma = 0.01 * np.sum(distsqX)/(self.n*self.n)
        else:
            self.sigma = self.params['bandwidth']
           
        self.maxIter = 100000
        if 'maxIter' in self.params:
            self.maxIter = self.params['maxIter']
        
    def fit(self):
        # train the principal tree structure
        maxIter = self.maxIter
        MU = self.MU
        history = self.history
        model = self.model
        X = self.x
        K = self.k
        sigma = self.sigma
        epsilon = 0.001
        for i in range(maxIter):
            distsqMU = sqdist(MU,MU)
            stree = minimum_spanning_tree(sparse.csr_matrix(sparse.tril(distsqMU))).toarray()
            stree = stree + stree.transpose()
            e = np.not_equal(stree,0)

            # store history
            history.MU.append(MU)
            history.stree.append(stree)
            
            # update data assignment matrix
            distMUX = sqdist(MU,X)
            min_dist = repmat(np.min(distMUX,0),K,1)
            tmp_distMUX = distMUX - min_dist
            tmp_R = np.exp(-tmp_distMUX.transpose()/sigma)
            R = tmp_R / repmat(np.sum(tmp_R,1), K,1).transpose()
            
            #compute objective function
            obj1 = - sigma * np.sum(np.log(np.sum(np.exp(-tmp_distMUX/sigma),0) ) - min_dist[0,:]/ sigma)
            reg = np.sum(np.sum(stree))
            obj = (obj1 + 0.5 * self.Lambda * reg)/self.n
            history.objs.append(obj)
            
            #projected mean square error
            #mse = sum(sum( R* distMUX.transpose()))
            #mse = obj1 + sigma *self.n* np.log(K)

            projd = np.min(distMUX,0)
            mse = np.mean(projd)

            history.mse.append(mse)
    
            #length of the structure
            history.length.append(reg)
            
            #terminate condition
            if i > 0:
                if abs((obj - old_obj)/old_obj) < epsilon:
                    #print '# of iters: ', i
                    break
                    
            L = np.diag(sum(e,0)) - e
            MU = np.dot(np.dot(X,R),np.linalg.inv(self.Lambda * L + np.diag(np.sum(R,0))))
            old_obj = obj
            
        self.history = history
        model.MU = MU
        model.stree = stree
        model.sigma = self.sigma
        model.Lambda = self.Lambda
        self.model = model
        
        return (history,model)
               
class History:
    def __init__(self):
        self.MU = []
        self.stree = []
        self.objs = []
        self.mse = []
        self.length = []
        
class Model:
    def __init__(self):
        self.MU = None
        self.stree = None
        self.sigma = None
        self.Lambda = None
                
        
def sqdist(a,b):
    aa = np.sum(a**2,0)
    bb = np.sum(b**2,0)    
    ab = np.dot(a.transpose(),b)
    d = abs(repmat(aa.transpose(),bb.shape[0],1).transpose() + repmat(bb.transpose(),aa.shape[0],1) - 2*ab)
    return d







