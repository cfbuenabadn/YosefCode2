import numpy as np
import principal_tree as pt


def gapstats(x, params,lower,upper,num_m):
    num_lambda = 10
    N_permutation = 10
    dim, N = x.shape
    Lambda = np.logspace(np.log10(lower),np.log10(upper),num_lambda)
        
    Obj_orginal = np.zeros((1,len(Lambda)))
    obj_perm  = np.zeros((len(Lambda),N_permutation))
    for i in range(len(Lambda)):
        params['lambda'] = Lambda[i]
        ptree= pt.principalTree(x,params,num_m)
        history, model = ptree.fit()
        Obj_orginal[0][i] = history.mse[-1]
    
    # start permutations    
    for s in range(N_permutation): 
        #print '<<<<< doing %d of %d permutation tests...' % (s, N_permutation)
        P_data = np.zeros((dim,N))
        for n in range(dim):
            P_data[n,:] = x[n,np.random.permutation(N)]
            
        for i in range(len(Lambda)):
            params['lambda'] = Lambda[i]
            pModel= pt.principalTree(P_data,params,num_m)
            history,a = pModel.fit()
            obj_perm[i,s] = history.mse[-1]



    Gap = np.mean(np.log(obj_perm),1) - np.log(Obj_orginal)

    Gap_std = np.std(np.log(obj_perm), axis=1)*np.sqrt(1+1.0/N_permutation)
    Gap_minusSE = Gap-Gap_std
    arry = Gap_minusSE[1:].tolist()
    arry.append(Gap_minusSE[-1])
    arry = np.asarray(arry)
    index = np.nonzero(Gap > arry.reshape(Gap.shape))[0]
 

    max_gap = np.max(Gap)
    max_idx = np.argmax(Gap)

    Lambda_opt = np.asarray([Lambda[max_idx],Lambda[index[0]]])   
    
    print 'gapstats done!'
    return Lambda_opt

