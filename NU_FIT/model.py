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
	def __init__(self, k=1, ct=0.001, 
		mt=200, rt = 1, bic=False,
		hist=200,m=0, kappa=0,alpha=-2,beta=0, BIC_PEN=10, maxBIC=3 ):
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
		#priors, default cancel each other and give MLE solution
		self.alphas 	= [1 for i in range(0, k*2)]
		self.m 			= m #for mu
		self.kappa 		= kappa #for mu
		self.alpha 		= alpha #for sigma
		self.beta 		= beta #for sigma






	def _estimate(self, X, Y,K):
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
		
		#pick i's
		IS 	= [np.random.uniform(min(X), max(X)) for i in range(0, K)]
		IS.sort()
		#pick s's 
		SS 	= [np.random.gamma(self.alpha, self.beta) for i in range(0, K)]
		#pick w's 
		initW	= 1.0 / (K*2)
		rvs = [normal(IS[i],SS[i],w=initW) for i in range(0, K)]
		rvs+= [uniform(IS[i], max(X),w=initW) for i in range(0, K)] 

		w 	= np.zeros((X.shape[0],K*2))

		while not converged and t< self.mt:
			#compute weight
			for k,rv in enumerate(rvs):
				w[:,k] 	= map(lambda x: rv.pdf(x), X)
			#compute likelihood
			LL 		= sum([LOG(sum([rv.pdf(x) for rv in rvs]))*y for x,y in zip(X,Y)])
			if abs(LL - prevLL) < self.ct:
				return prevLL, rvs,True
			prevLL 	= LL
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
			t+=1
		return prevLL, rvs,converged
	
	def fit(self, X, weights = None):
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
					LL, rvs,converged 	= self._estimate(X, Y, k)

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
						LL, rvs,converged 	= self._estimate(X, Y, k)
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
	#test bayesian priors
	#simulate
	D 		= [x for x in np.random.normal(10,1,100)]+ [x for x in np.random.uniform(10,100,100)]
	clf 	= NU(alpha=2,beta=2, m=0, kappa=0.01,bic=True)
	clf.fit(D)
	plt.hist(D, bins=100, normed=1)
	xs 	= np.linspace(min(D), max(D), 1000)
	plt.plot(xs, [clf._func(x) for x in xs], linewidth=4., linestyle="--")
	plt.show()










