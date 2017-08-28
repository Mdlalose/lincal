#simulating visibilities
from astropy.io import fits
from astropy import units as u
import numpy as np
import math
import anta_pos as ant
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import healpy as hp
from astropy.io import fits
#import generate_gaussian_data as gdata
# read map
#map_I = hp.read_map('GLEAM_EGC.fits')
#hp.mollview(map_I, coord=['G','E'], title='Histogram equalized Ecliptic', unit='mK', norm='hist', min=-1,max=1, xsize=2000)

RA_src, Dec_src = np.load('src_pos.npy')[0], np.load('src_pos.npy')[1]
src_flux = np.load('src_flux.npy')[1]

xdim = int(sys.argv[1])
ydim = int(sys.argv[2])
long_0= 21.4278 # deg
lalt_0=  -30.7224 # deg
ant_positions= ant.get_antenna_pos(xdim,ydim,long_0,lalt_0)
#x,y = ant.get_antenna_pos(xdim,ydim,long_0,lalt_0)
zen_vec = [0.0,1.0,0.0] #np.round(ant_positions[0]/np.linalg.norm(ant_positions[0]))





x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))
xyz = np.concatenate((x,y,z))
mydot = x*zen_vec[0]+ y*zen_vec[1] + z*zen_vec[2]

theta = np.arccos(mydot)
sigma = 0.5*(1.22*(0.5/6.0))
beam = np.exp(-0.5*(theta/sigma)**2)
ii =  beam > 0.1

RA_src= RA_src[ii]
Dec_src= Dec_src[ii]
src_flux = src_flux[ii]

# ants epsilon

ant_eps = np.load('Ant_epsolon0.05.npy')[0]
blVectors, blPairs = ant.get_blvectors(ant_positions)
	#np.save('redund_blVectors.npy',[blVectors,blPairs])
	

ublDict = {}
for b in range(len(blVectors)):
     	# grouping togather all blPairs with same blVector
    	if ublDict.has_key(blVectors[b]): ublDict[blVectors[b]].append(blPairs[b])
    	else: ublDict[blVectors[b]] = [blPairs[b]]

ublIndexDict = {antPair: i for i,antPairs in enumerate(ublDict.values()) for antPair in antPairs }
ublVectors = np.array([ant_positions[antList[0][0]]-ant_positions[antList[0][1]] for antList in ublDict.values()])



print "With", len(ant_positions), "antennas there are", len(ublDict.items()), "unique baselines."
#return [ublVectors,ublIndexDict]

	
zen = zen_vec 
ant1, ant2 =[],[]
for k in range(len(blPairs)):
		ant1.append(np.array(blPairs)[k][0])
		ant2.append(np.array(blPairs)[k][1])

vis_map =[]
for i in range(len(ant1)):
		vis_map.append(ublIndexDict[(ant1[i],ant2[i])])

x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))
xyz = np.concatenate((x,y,z))
mydot = x*zen_vec[0]+ y*zen_vec[1] + z*zen_vec[2]

theta = np.arccos(mydot)
beam = np.exp(-0.5*(theta/sigma)**2)
			 
			
plt.plot(theta,beam,'.')
plt.ylabel('Power Response')
plt.xlabel(r'$\theta$ [radians]')





vis_ = np.load('obstrue_vis_8x8redundant_data.npy')
vis_b_src = np.load('obstrue_vis_8x8redundant_data_all_bright_src.npy')
vis = vis_-vis_b_src
			

def unique_blocks(data,blVectors,ublVectors,nu,ant_eps):
        block_vis_nu =[]
	block_corr_nu=[]
	
	for nu_i in range(nu.size):
		Red_block=[]
		Red_block_vis=[]
		ant_eps_i = ant_eps[nu_i]
		data_i = data[nu_i]
    		for u_uniq in ublVectors:
        		tmp = []
			tmp_1 =[]
        		for i in range(len(blVectors)):
           	 		if blVectors[i][0] == u_uniq[0] and blVectors[i][2]==u_uniq[2] :
					sigma_ij = 1.0/(sigma*(1.0 + 0.0*ant_eps_i[ant1[i]]))**2 + 1.0/(sigma*(1.0 +0.0*ant_eps_i[ant2[i]]))**2
					
               		  		tmp.append(sigma_ij)
					tmp_1.append(data_i[i])	
		 	Red_block.append(tmp)
			Red_block_vis.append(tmp_1)
		block_corr_nu.append(Red_block)
		block_vis_nu.append(Red_block_vis)
		
    	return [block_corr_nu,block_vis_nu]

nu = np.arange(170.0,230,1.0)
# computing visibility correlation for 170MHz-230MHz

N_BLOCK_all_freq=[]
N_BLOCK_all_freq_vis=[]
block_matrix = unique_blocks(vis,blVectors,ublVectors,nu,ant_eps)


for nu_i in range(nu.size):

	N_block_corr=[]
	N_block_corr_vis=[]
	for n_block in range(len(block_matrix[0][nu_i])):

		#correlation matrix
        	corr_sigma = np.zeros((len(block_matrix[0][nu_i][n_block]),len(block_matrix[0][nu_i][n_block])))
		corr_vis_ = np.zeros((len(block_matrix[0][nu_i][n_block]),len(block_matrix[0][nu_i][n_block])),dtype='complex')
		
		for vis_i in range(len(block_matrix[0][nu_i][n_block])):
			for vis_j in range(len(block_matrix[0][nu_i][n_block])):
				
				#computing correlation for vis i and j
				corr_sigma[vis_i][vis_j] = np.pi*np.power(1.0/2.0*block_matrix[0][nu_i][n_block][vis_i] + 1.0/2.0*block_matrix[0][nu_i][n_block][vis_j],-1)
				

				corr_vis_[vis_i][vis_j] = np.conj(block_matrix[1][nu_i][n_block][vis_i])*block_matrix[1][nu_i][n_block][vis_j]
				
		N_block_corr.append(corr_sigma)
		N_block_corr_vis.append(corr_vis_)
	N_BLOCK_all_freq.append(N_block_corr)
	N_BLOCK_all_freq_vis.append(N_block_corr_vis)




#computing gamma

gamma_all_freq=[]
for nu_i in range(nu.size):

	sigma_ab_sum =[]
	vis_ab_sum = []
	for k_block in range(len(N_BLOCK_all_freq_vis[nu_i])):
		sigma_ab_sum.append(np.sum(np.matrix.flatten(np.matrix(N_BLOCK_all_freq_vis[nu_i][k_block]))))
		vis_ab_sum.append(np.sum(np.matrix.flatten(np.matrix(N_BLOCK_all_freq[nu_i][k_block]))))
	
	gamma_all_freq.append(np.sum(vis_ab_sum)/np.sum(sigma_ab_sum))



#N_BLOCK_all_freq,gamma_freq = np.load('block_corr_vis_gamma_ij_8x8.npy')[0],np.load('block_corr_vis_gamma_ij_8x8.npy')[2]
Block_Eigenvec=[]
for nu_i in range(len(N_BLOCK_all_freq)):
	Eigenvect_nu_i=[]
	for block_i in range(len(N_BLOCK_all_freq[nu_i])):
		c_ii = np.real(np.power(gamma_all_freq[nu_i],-1))*N_BLOCK_all_freq[nu_i][block_i]
		Eigenvect_nu_i.append(np.linalg.eig(c_ii))
	Block_Eigenvec.append(Eigenvect_nu_i)	
	





