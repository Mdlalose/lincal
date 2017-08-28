import numpy as np
qq=np.load('Block_eigevec_beamvar_all_freq.npy')
nblock=qq[0].shape[0]


#first, find out how many eigenmodes we need to keep in the largest block
#search for the maximum eigenvalues within a block, nmax = number of eigenvalues above eigenvalue_threst.

edges=np.zeros(nblock+1,dtype='int')
q_all_freq=[]
for nu_i in range(len(qq)):
	nmax=0
	thresh=1e-8
	nvis=0
	for i in range(nblock):
	    	myeig=qq[nu_i][i][0]
    	    	nkept=np.sum(myeig>thresh*myeig.max())
    	    	nvis=nvis+myeig.size
    		edges[i+1]=edges[i]+myeig.size
    		if nkept>nmax:
        		nmax=nkept

	print 'total number of visibilities is ' + repr(nvis) + ' with ' + repr(nmax) + ' kept eigenvalues at max.'

	# computeting q = vl^1/2, l_above_threshold
	q=np.zeros([2*nvis,nmax*2])
	myeig_vec =[]
	for i in range(nblock):
    		myeig=np.real(qq[nu_i][i][0])
    		myvecs=np.real(qq[nu_i][i][1])
    		ind=myeig>thresh*myeig.max() # picking up an index
    		#print 'index max eiegen', ind
    		myeig_use=myeig[ind] # pict max eigenvalue
    		#myeig_vec.append(myeig_use)
   		myvecs_use=myvecs[:,ind] # pick  corresponding max eigenvec
   		for j in range(len(myeig_use)):
        		myvecs_use[:,j]= myvecs_use[:,j]*np.sqrt(myeig_use[j])
        		q[2*edges[i]:2*edges[i+1]:2,2*j]= 1.0 # myvecs_use[:,j]
        		q[(2*edges[i]+1):2*edges[i+1]:2,2*j+1]=  0.0  #myvecs_use[:,j]

	q_all_freq.append(q)	
        
