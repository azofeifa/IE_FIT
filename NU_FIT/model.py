import math as m
import numpy as np
np.seterr(all="ignore")
global no_plot
try:
	import matplotlib.pyplot as plt
	no_plot=False
except:
	no_plot=True 

def LOG(x):
	if x > 0:
		return m.log(x)
	return -np.Inf
def testRealData(FILE="/Users/joeyazo/Desktop/Lab/gro_seq/files/bed_files/DMSO2_3.sorted.fiveprime.pos.BedGraph",chrom="chr1",
	start=84533886, stop=84704619):
	FH 	= open(FILE)
	X,Y =list(),list()
	F 	= False
	for line in FH:
		c,st, sp, coverage 	= line.strip("\n").split("\t")
		st, sp, coverage 	= int(st), int(sp), int(coverage)
		if chrom == "chr1" and st > start and sp < stop:
			F 	= True
			X+=[st for i in range(st, sp)]
			Y+=[coverage for i in range(st, sp)]
		elif F:
			break
	FH.close()
	return X,Y


class uniform:
	def __init__(self, a, b, w=1.0):
		self.a 		= a
		self.b 		= b
		self.w 		= w
		self.type	= "uniform"
	def __str__(self):
		return "uniform, a: " + str(self.a) + ", b: " + str(self.b)+ ", weight: " + str(self.w)
	def pdf(self, x):
		return self.w*(1.0 / (self.b - self.a))*int(self.a<=x<=self.b)
class normal:
	def __init__(self, mu, sigma,w=1.0):
		assert sigma > 0, "sigma value must be positive"
		self.mu 	= mu
		self.sigma 	= sigma
		self.w 		= w
		self.type 	= "normal"
	def __str__(self):
		return "normal, mu: " + str(self.mu) + ", sigma: " + str(self.sigma)+ ", weight: " + str(self.w)
	def pdf(self, x):
		return (self.w / (self.sigma*m.sqrt(2*m.pi)))*m.exp(-m.pow(x-self.mu,2)/(2*self.sigma**2))
class NU:
	def __init__(self, k=1, ct=0.0001, 
		mt=100, rt = 1, bic=False,
		hist=500,m=0, kappa=0,alpha=-2,beta=0, BIC_PEN=10, maxBIC=3, split =False,gibbs=False ):
		self.k			= k
		self.rvs 		= None
		self.ct 		= ct
		self.mt 		= mt
		self.rt 		= rt
		self.bic 		= bic
		self.hist 		= hist
		self.maxBIC 	= maxBIC
		self.rvs 		= None
		self.LL 		= None
		self._params 	= None
		self.converged 	= False
		self.BIC_PEN 	= BIC_PEN
		self.split 		= split
		#priors, default cancel each other and give MLE solution
		self.alphas 	= [100 for i in range(0, k*2)]
		self.m 			= m #for mu
		self.kappa 		= kappa #for mu
		self.alpha 		= alpha #for sigma
		self.beta 		= beta #for sigma
		self.gibbs 		= gibbs





	def _estimate(self, X, Y,K,rev=False):
		if K == 0: #assume just a uniform distribution
			a,b 	= min(X), max(X)
			ll 		= sum([LOG(1. / (b-a))*y for y in Y if y])
			return ll, [uniform(a, b)],True

		prevLL,t,converged 		= 0,0, False
		#==============================
		#set the dirchlet prior
		self.alphas 			= [1 for k in range(0,2*K) ]
		#==============================
		# random initialize
		#==============================
		maxX 					= max(X)
		minX 					= min(X)
		#pick i's
		IS 	= [np.random.uniform(min(X), maxX) for i in range(0, K)]
		IS.sort()
		#pick s's 
		SS 	= [1000 for i in range(0, K)]
		#pick w's 
		initW	= 1.0 / (K*2)
		rvs = [normal(IS[i],SS[i],w=initW) for i in range(0, K)]
		if rev:
			rvs+= [uniform(minX, IS[i],w=initW) for i in range(0, K)] 	
		else:
			rvs+= [uniform(IS[i], maxX,w=initW) if i == 0 else uniform(IS[i], np.random.uniform(IS[i], maxX),w=initW)  for i in range(0, K)] 
		w 	= np.zeros((X.shape[0],K*2))
		prevLL 	= -np.Inf
		while not converged and t< self.mt:
			#compute weight
			for k,rv in enumerate(rvs):
				w[:,k] 	= map(lambda x: rv.pdf(x), X)
			#compute likelihood
			#normalize weights
			for i in range(0, w.shape[0]):
				if sum(w[i,:]): #don't want to divide by zero
					w[i,:] /= sum(w[i,:])
				
			#compute new pi, weight sample wean, weighted sample variance
			for k, rv in enumerate(rvs):
				rv.w 		= (w[:,k].dot(Y) + self.alphas[k] - 1) / (sum(Y) + sum(self.alphas) - K)
				if hasattr(rv, "mu"):
					uRVa 	= [RV for RV in rvs if hasattr(RV, "a") and getattr(RV, "a") == rv.mu]
					rv.mu 	= ((w[:,k]*X).dot(Y) + self.kappa*self.m) / (sum(w[:,k]*Y)+self.k)
					if uRVa:	uRVa[0].a 	= rv.mu
					rv.sigma= m.sqrt(( (w[:,k]*(X-rv.mu)**2).dot(Y) + 2*self.beta + self.kappa*pow(rv.mu-self.m,2))   / ((sum(w[:,k]*Y)) + self.alpha + 2))
					if rv.sigma < 0.0001: #one of the components blew up ):
						return -np.Inf, rvs,False
			

			LL 		= sum([LOG(sum([rv.pdf(x) for rv in rvs]))*y for x,y in zip(X,Y)])
			if self.split:
				#============
				# instead of searching through all Ls make it gready...
				# #pick new Ls, run accross data pick best split based on max loglikelihood
				uniforms 	= [rv for rv in rvs if rv.type == "uniform"]
				prevBs 		= [rv.b for rv in rvs if rv.type == "uniform"]
				maxUniL 	= LL
				argU 		= None
				for u in uniforms:
					if u.b != maxX:
						trys 	= np.linspace(u.b-10, u.b+10,5)
						for l in trys:
							u.b 	= l
							ull 	= sum([LOG(sum([rv.pdf(x) for rv in rvs]))*y for x,y in zip(X,Y)])
							if ull > maxUniL:
								maxUniL = ull
								argU 	= [rv.b for rv in rvs if rv.type == "uniform"]
				if argU:
					for i, argl in enumerate(argU):
						uniforms[i].b 	= argl
					LL 			= maxUniL
				else:
					for i, b in enumerate(prevBs):
						uniforms[i].b 	= b
				if abs(LL - prevLL) < self.ct:
					print "converged"
					return prevLL, rvs,True
			prevLL 	= LL	
			t+=1
		return prevLL, rvs,converged
	def _estimate_gibbs(self,X,Y,k,rev=False):
		ll 			= -np.Inf
		rvs 		= list()
		converged 	= False

		return ll,rvs, converged
		


	
	def fit(self, X, weights = None,rev=False):
		if self.hist is not None:
			counts,edges 	= np.histogram(X, bins=self.hist,weights=weights)
			X 				= edges[1:]
			Y 				= counts
		else:
			X 				= np.array(X)
			Y 				= np.array([1. for x in X])
		for ITER in range(0, self.rt):
			if not self.bic:
				k 							= self.k
				maxLL, maxRVs,maxC 	= -np.Inf, None,False
				for t in range(0, self.rt):
					if not self.gibbs:
						LL, rvs,converged 	= self._estimate(X, Y, k,rev=rev)
					else:
						LL,rvs, converged 	= self._estimate_gibbs(X,Y,k,rev=rev)

					if maxLL < LL or maxLL == -np.Inf:
						maxLL 	= LL
						maxRVs 	= rvs 
						maxC 	= converged
			else:
				n 			= sum(Y)
				R 			= list()
				for k in range(0, self.maxBIC):
					maxLL, maxRVs,maxC 	= -np.Inf, None,False
					for t in range(0, self.rt):
						if not self.gibbs:
							LL, rvs,converged 	= self._estimate(X, Y, k,rev=rev)
						else:
							LL,rvs, converged 	= self._estimate_gibbs(X,Y,k,rev=rev)

						if maxLL < LL or maxRVs is None:
							maxLL 	= LL
							maxRVs 	= rvs 
							maxC 	= converged
					
					R.append((maxLL, maxRVs, k, maxC))
				bic_scores 			= np.array([-2. * ll + ((k+1)*(k*self.BIC_PEN)) *LOG(n) for ll,rvs, k,c in R])
				argmin 				= bic_scores.argmin()
				maxLL, maxRVs, k,maxC 	= R[argmin]
			
			self.rvs 	= maxRVs
			self.LL 	= maxLL
			self.k 		= k
			self.converged = maxC
			H 			= {"weights": [rv.w for rv in self.rvs],
			"mu": [rv.mu for rv in self.rvs if hasattr(rv, "mu")],
			"sigma": [rv.sigma for rv in self.rvs if hasattr(rv, "mu")], "converged": maxC, "rvs": maxRVs}
			self._params = H
			return H

	def display(self, X, weights=None):
		if not no_plot:
			counts,edges 	= np.histogram(X, bins=200,weights=weights, normed=1)
			X 				= (edges[1:] + edges[:-1]) / 2.
			Y 				= counts
			xs 				= np.linspace(min(X)-10, max(X)+10, 1000)
			for rv in self.rvs:
				print rv
			plt.bar(X,Y,alpha=0.5,width=(max(X) - min(X)) / len(X))
			p1, 	= plt.plot(xs, map(lambda x: sum([rv.pdf(x) for rv in self.rvs ]),  xs), linewidth=2.5, linestyle = "--", color='green')
			plt.show()


	def predict(self, x):
		P 		= np.array([rv.pdf(x) for rv in self.rvs])
		P 		/= sum(P) 
		return P.argmax()
		
	def _func(self, x):
		return sum([rv.pdf(x) for rv in self.rvs])
if __name__ == "__main__":
	D 	 = [i for i in np.random.normal(0,1,1000)] + [i for i in np.random.uniform(0,100,1000)]
	D 	+= [i for i in np.random.normal(20,1,500)] + [i for i in np.random.uniform(20,75,1000)]

	Y,X = np.histogram(D, bins=1000)
	X 	= (X[:-1] + X[1:])/2.
	clf = NU(hist=1000,gibbs=False, k=2, split=True)
	clf.fit(X, weights=Y)
	clf.display(D)





