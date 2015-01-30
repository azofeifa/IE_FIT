import  model
import multiprocessing as mp
import numpy as np

global no_plot
LST 	= list()
try:
	import matplotlib.pyplot as plt
	no_plot=False
except:
	no_plot=True 
def accumulateResults(result):
	LST.append(result)

def wrapper(h,clf,i):
	center 	= min(h.X)
	X 	= [x-center for x in h.X]
	H 		= clf.fit(X,weights=h.Y)
	return H,i
def run(H,np=8,maxBIC=None, penality=None,rt=3,binSize=200):
	
	pool = mp.Pool(processes=np) 
	for t,i in enumerate(H.values()):
		if maxBIC:
			clf 	= model.NU(bic=True,rt=rt,alpha=2,beta=500,BIC_PEN=penality, maxBIC=maxBIC,hist=binSize)
		else:
			clf 	= model.NU(bic=False,rt=rt,alpha=2,beta=500,hist=binSize)
			
		pool.apply_async(wrapper, args=( i,clf,t  ), callback=accumulateResults)
	pool.close()
	pool.join()

	for p in LST:
		h,t 		= p	
		rvs 		= h["rvs"]
		I 			= H.values()[t]
		I.rvs 		= rvs
	return H