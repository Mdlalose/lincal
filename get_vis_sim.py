#simulating visibilities
from astropy.io import fits
from astropy import units as u
import numpy as np
import numpy
import math
import anta_pos as ant
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import healpy as hp
from astropy.io import fits
import heapq

# source catalogue (RA,Dec) & source flux
RA_src, Dec_src = np.load('src_pos.npy')[0], np.load('src_pos.npy')[1]
src_flux = np.load('src_flux.npy')[1]

xdim = int(sys.argv[1])
ydim = int(sys.argv[2])
#location of the array, SKA site, Karoor, South Africa
long_0= 21.4278 # deg
lalt_0=  -30.7224 # deg
ant_positions= ant.get_antenna_pos(xdim,ydim,long_0,lalt_0) 
#x,y = ant.get_antenna_pos(xdim,ydim,long_0,lalt_0)
zen_vec = zen_vec = [np.cos(np.radians(lat_0))*np.cos(np.radians(long_0)), np.cos(np.radians(lat_0))*np.sin(np.radians(long_0)), np.sin(np.radians(lat_0))] #np.round(ant_positions[0]/np.linalg.norm(ant_positions[0]))

"""
# plot the antenna position and redundant baselines
plt.figure()
plt.plot(x,y,'*')
plt.plot([x[0][0],x[1][0]],[y[-1][0],y[-1][0]],'k')
plt.text(x[0][0],y[-1][0],'0')
plt.text(x[1][0],y[-1][0],'1')
plt.text(x[2][0],y[-1][0],'2')
plt.plot([x[0][0],x[1][0]],[y[-1][0],y[-1][0]],'k')
plt.text(x[0][0],y[3][0],'3')
plt.text(x[1][0],y[3][0],'4')
plt.text(x[2][0],y[3][0],'5')
plt.text(x[0][0],y[0][0],'6')
plt.text(x[1][0],y[0][0],'7')
plt.text(x[2][0],y[0][0],'8')

plt.plot([x[0][0],x[2][0]],[y[-1][0],y[-1][0]])
plt.plot([x[0][0],x[0][0]],[y[3][0],y[-1][0]])
plt.plot([x[0][0],x[1][0]],[y[-1][0],y[3][0]])
plt.plot([x[0][0],x[2][0]],[y[-1][0],y[3][0]])
plt.plot([x[0][0],x[0][0]],[y[-1][0],y[0][0]])
plt.plot([x[0][0],x[1][0]],[y[-1][0],y[0][0]])
plt.plot([x[0][0],x[2][0]],[y[-1][0],y[0][0]])
plt.plot([x[0][0],x[2][0]],[y[0][0],y[3][0]])
plt.plot([x[0][0],x[1][0]],[y[0][0],y[3][0]])
plt.plot([x[0][0],x[2][0]],[y[0][0],y[-1][0]])
plt.plot([x[0][0],x[1][0]],[y[0][0],y[-1][0]])


plt.xlabel('East-West Antenna Position [m]')
plt.ylabel('North-South Position [m]')
plt.grid()


"""

# selecting source within 32 deg2

x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))
xyz = np.concatenate((x,y,z))
mydot = x*zen_vec[0]+ y*zen_vec[1] + z*zen_vec[2]

theta = np.arccos(mydot)

sigma = 0.5*(1.22*(0.5/6.0))
beam = np.exp(-0.5*(theta/sigma)**2)
beam[theta > np.pi/2.0] =0
ii= beam > 0.1
x=x[ii]
y=y[ii]
z=z[ii]
src_flux = src_flux[ii]
xyz = np.concatenate((x,y,z))
mydot = x*zen_vec[0]+ y*zen_vec[1] + z*zen_vec[2]

theta = np.arccos(mydot)
beam = np.exp(-0.5*(theta/sigma)**2)

#selecting bright sources
myflux = beam*src_flux
flux_10b = np.sort(myflux)[-10:]
myflux_sorted = np.sort(myflux)
ff=np.cumsum(myflux_sorted)
ff=ff/ff[-1]
nsrc=np.sum(ff>0.75)
myvar=np.sort(myflux**2)

flub = np.sort(myflux)

print 'Total number of sources', nsrc
nsim=100
v1=numpy.zeros(nsim,dtype=complex)
v2=numpy.zeros(nsim,dtype=complex)
v3=numpy.zeros(nsim,dtype=complex)
for j in range(nsim):
    myphase=numpy.exp(2*numpy.pi*numpy.random.rand(myflux.size)*numpy.complex(0,1))
    v1[j]=numpy.sum(myphase*flub)
    v2[j]=numpy.sum(myphase[-10:]*flub[-10:])
    v3[j]=v1[j]-v2[j]
vv=numpy.real(numpy.asarray([numpy.mean(numpy.conj(v1)*v1),numpy.mean(numpy.conj(v2)*v2),numpy.mean(numpy.conj(v3)*v3)]))


#print numpy.mean(numpy.conj(v1)*v1),numpy.mean(numpy.conj(v2)*v2),numpy.mean(numpy.conj(v3)*v3)
"""
print numpy.sqrt(vv)

RA_b, Dec_b, src_b =[], [], []
for src in range(len(flux_10b)):
	for i in range(len(src_flux)):
		if np.round(myflux[i],5) == np.round(flux_10b[src],5):
			# np.round(flux_10b[src],5), np.round(flux_[i],5)
			RA_b.append(RA_src[i])
			Dec_b.append(Dec_src[i])
			src_b.append(src_flux[i])

		else:
			pass


np.save('bright_src_data.npy',[RA_b, Dec_b, src_b])
"""

def get_V_uvw(blVectors,RA_src,Dec_src,src_flux,l_wave,dia_m,zen,nu_min,nu_max,ant1,ant2,ant_eps):
	nu = np.arange(nu_min,nu_max,1.0)
	nu_0 = (nu_max-nu_min)/2.0
	V_nu = []
	myvis=[]

	expec_err =[]

	sigma = 0.5*(1.22*(l_wave/dia_m))
	for nu_i in range(len(nu)):
		#print (nu[nu_i]/nu_0)**-2.5
		V_uvw =   np.zeros(len(blVectors), dtype = "complex")
		#V_uvw_error =  np.zeros(len(blVectors), dtype = "complex")
		
		
		#x,y, z = np.random.randn(nsrc), np.random.randn(nsrc), np.random.randn(nsrc)
		x, y, z = np.cos(np.radians(Dec_src))*np.cos(np.radians(RA_src)), np.cos(np.radians(Dec_src))*np.sin(np.radians(RA_src)), np.sin(np.radians(Dec_src))
		r=np.sqrt(x**2 + y**2 + z**2)
		mydot = (x/r)*zen[0] + (y/r)*zen[1] + (z/r)*zen[2]
  	        theta = np.arccos(mydot)
  	        ant_eps_i = ant_eps[nu_i]
		#print ant_eps_i.size, nu[nu_i]

	
		for k in range(len(blVectors)):
			 
			
		         
			 beam12 = np.exp(-(0.5*theta**2)*(1.0/(sigma*(1.0 + ant_eps_i[ant1[k]]))**2 + 1.0/(sigma*(1.0 +ant_eps_i[ant2[k]]))**2))
			 beam = np.exp(-(theta/sigma)**2)
			 """
			 #plt.title('0.001$\epsilon$')
			 plt.plot(theta,beam,'.')
			 plt.ylabel('Power Response')
			 plt.xlabel(r'$\theta$ [radians]')
			 #plt.plot(theta,beam12,'.',label='Imperfection Beams')
			 #plt.legend(loc='best')
			 plt.grid()
	
			 
                         
                         """
  			 r=np.sqrt(x**2 + y**2 + z**2)
  			 
			 mydotxyz =  (x/r)*blVectors[k][0] + (y/r)*blVectors[k][1] + (z/r)*blVectors[k][2]
			 exp_phase_shift = np.exp(-2j*np.pi*(mydotxyz/l_wave))
			 tmp= np.sum(beam*src_flux*np.power(nu[nu_i]/nu_0,-2.5)*exp_phase_shift)
			 V_uvw[k]=tmp
			 # first order error in visibility due to imperfections in antenna locations
			 #delta_blVectors= np.array(blVectors_err[bl]) - np.array(blVectors[bl])
			 #mydotxyz_err =  (x/r)*delta_blVectors[0] + (y/r)*delta_blVectors[1] + (z/r)*delta_blVectors[2]
			 #tmp_err = -1j*np.sum(beam*src_flux*((2.0*np.pi*mydotxyz_err)/l_wave))

			 #V_uvw_error[bl] = tmp_err
		
			 
		V_nu.append(V_uvw)
	

	return np.array(V_nu)
			 
		        	


Ant_eps = np.load('Ant_epsolon0.05.npy')[0] # primary beam errors
blVectors, blPairs = ant.get_blvectors(ant_positions)
	#np.save('redund_blVectors.npy',[blVectors,blPairs])
	
#computing baselines and unique base vectors
ublDict = {}
for b in range(len(blVectors)):
     	# grouping togather all blPairs with same blVector
    	if ublDict.has_key(blVectors[b]): ublDict[blVectors[b]].append(blPairs[b])
    	else: ublDict[blVectors[b]] = [blPairs[b]]

ublIndexDict = {antPair: i for i,antPairs in enumerate(ublDict.values()) for antPair in antPairs }
ublVectors = np.array([ant_positions[antList[0][0]]-ant_positions[antList[0][1]] for antList in ublDict.values()])



print "With", len(ant_positions), "antennas there are", len(ublDict.items()), "unique baselines."
#return [ublVectors,ublIndexDict]

	
   
# computing antenna indeces and visibility map indeces    
ant1, ant2 =[],[]
for k in range(len(blPairs)):
		ant1.append(np.array(blPairs)[k][0])
		ant2.append(np.array(blPairs)[k][1])

vis_map =[]
for i in range(len(ant1)):


		vis_map.append(ublIndexDict[(ant1[i],ant2[i])])

#Vis_sim_data = get_V_uvw(blVectors,RA_src, Dec_src, src_flux,0.5,6.0,zen_vec,170,230,ant1,ant2,Ant_eps)

"""
Vis_bright=[]
for src_i in range(len(RA_b)):
	Vis_sim_data = get_V_uvw(blVectors,RA_b[src_i], Dec_b[src_i], src_b[src_i],0.5,6.0,zen_vec,170,230,ant1,ant2,Ant_eps)
	Vis_bright.append(Vis_sim_data)
	
"""
#np.save('obstrue_vis_' + repr(xdim) + 'x' + repr(ydim)  +'beamvar_data_all_freq_10_bright_src'+'.npy',Vis_bright)


"""
#RA_src, Dec_src, src_flux = np.load('10_brightest_src.npy')[0], np.load('10_brightest_src.npy')[1], np.load('10_brightest_src.npy')[2]
#sim_vis_n_sample =[]
for n_i in range(1):

	visTrue = get_V_uvw(blVectors,RA_src,Dec_src,src_flux,0.5,6.0,zen_vec,170,230,ant1,ant2,Ant_eps[n_i])
	sim_vis_n_sample.append(visTrue)
	




#np.save('obstrue_vis_' + repr(xdim) + 'x' + repr(ydim)  + '20_data_redundant'+ '.npy',[sim_vis_n_sample,ant1,ant2,vis_map])
"""

#simulating visibility position errors
"""
error = float(sys.argv[3])
n_sample = int(sys.argv[4])
ant_1=[]
sg=[]
bl_=[]
v1 = np.load('obstrue_vis_8x80.0.npy')[0]
#bl_red = np.load('redund_blVectors.npy')[0]
for n in range(n_sample):

	ant_positions_err= np.random.randn(ant_positions.shape[0],ant_positions.shape[1])
	#print ant_positions_err[0][2]

	ant_positions_err[:,1] = 0.0
	ant_positions_err = ant_positions  + error*ant_positions_err
        ant_p = np.array(ant_positions_err)
	nAntennas = len(ant_p)

	
       

	blVectors, blPairs = ant.get_blvectors(ant_positions_err)
	#np.save('redund_blVectors.npy',[blVectors,blPairs])
	

	ublDict = {}
	for b in range(len(blVectors)):
     		# grouping togather all blPairs with same blVector
    		if ublDict.has_key(blVectors[b]): ublDict[blVectors[b]].append(blPairs[b])
    		else: ublDict[blVectors[b]] = [blPairs[b]]

	ublIndexDict = {antPair: i for i,antPairs in enumerate(ublDict.values()) for antPair in antPairs }
	ublVectors = np.array([ant_positions_err[antList[0][0]]-ant_positions_err[antList[0][1]] for antList in ublDict.values()])



	print "With", len(ant_positions_err), "antennas there are", len(ublDict.items()), "unique baselines."
	#return [ublVectors,ublIndexDict]

	
        
	ant1, ant2 =[],[]
	for k in range(len(blPairs)):
		ant1.append(np.array(blPairs)[k][0])
		ant2.append(np.array(blPairs)[k][1])

	vis_map =[]
	for i in range(len(ant1)):
		vis_map.append(ublIndexDict[(ant1[i],ant2[i])])
	
	
	visTrue = get_V_uvw(blVectors,RA_src,Dec_src,src_flux,0.5,6.0,zen_vec,170,230,ant1,ant2,ant_eps)
	
	print np.std(visTrue-v1, axis=1)/np.std(v1, axis=1), np.exp(2j*np.pi*error*(0.05/0.5)).imag
	np.save('obstrue_vis_' + repr(xdim) + 'x' + repr(ydim) + repr(error) + 'beam' +'.npy',[visTrue,ant1,ant2,vis_map])
	
	#print np.linalg.norm(v1-visTrue[0]),  np.exp(2j*np.pi*error*(0.05/0.5)).imag
	if  np.mean(np.std(visTrue-v1, axis=1)/np.std(v1, axis=1)) <= np.exp(2j*np.pi*error*(0.05/0.5)).imag*1.3:
		break

	#print blVectors[2]

        #sg.append(np.std(visTrue-v1)/np.std(v1))
	#bl_.append(blVectors[6][0])
        
	print " successfully calculated visibilities"
"""
	
"""
for pos in range(len(ant_positions)):
	plt.title('$\mathbf{r} + \delta \mathbf{r}(\sigma=0.05)$')
	plt.plot(ant_positions[pos][0],ant_positions[pos][1], '*', label = 'Perfect Redundant')
	#plt.plot(ant_positions_err[pos][0],ant_positions_err[pos][1], '.', label = 'Near-Redundant')
	plt.xlabel('E-W [m]')
	plt.ylabel('N-S [m]')
	plt.grid()

"""


"""

# grouping visibilities according to redudndant blocks
block =np.load('block_index_edges_8x8.npy')[0]
Vis_src = Vis_bright
Vis_10_brght_src =[]
for src_j in range(len(RA_b)):
	vis_brgt_grp_frq=[]

	for freq in range(len(Vis_src[src_j])):
	        tmp_1 =[]
	 	for bl_i in range(len(block)):
			
			tmp_1.append(Vis_src[src_j][freq][block[bl_i]])
	
	 	
	 	vis_brgt_grp_frq.append(tmp_1)
	Vis_10_brght_src.append(vis_brgt_grp_frq)


	

# grouping antenna indeces according to redundant blocks
ant_1_grp=[]
ant_2_grp =[]
ant1 = np.array(ant1)
ant2 = np.array(ant2)
for bl_i in range(len(block)):
		ant_1_grp.append(ant1[block[bl_i]])
		ant_2_grp.append(ant2[block[bl_i]])



Vis_all_src= Vis_10_brght_src

#Separating visibility  as real and imaginary [r1,i1,r2,i2...]
Block_all_10_src=[]
for src_k in range(len(Vis_10_brght_src)):
	block_freq_vis=[]
	for j in range(len(Vis_all_src[src_k])):

		block_vis=[]
		for p in range(len(Vis_all_src[src_k][j])):
			tmp_1 =[]
			for vis_i in range(len(Vis_all_src[src_k][j][p])):

		
				tmp_1.append(np.real(Vis_all_src[src_k][j][p][vis_i]))
				tmp_1.append(np.imag(Vis_all_src[src_k][j][p][vis_i]))
			
		
			block_vis.append(tmp_1)	

		block_freq_vis.append(np.concatenate(block_vis).ravel())
	Block_all_10_src.append(block_freq_vis)

	

"""



				
