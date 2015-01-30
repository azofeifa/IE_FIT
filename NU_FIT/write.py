import numpy as np
def predict(rvs, x):
	norm 	= sum([rv.pdf(x) for rv in rvs])
	probs 	= np.array([rv.pdf(x) / norm for rv in rvs])
	return rvs[probs.argmax()]



def writeIGV(H, OUT,strand,D):
	ID 		= "IE_" + str(D["-chr"]) +"_" +str(D["-BIC"]) + "_"+str(D["-rt"])
	header 	= "track name=" + ID+ "2015-01-08 17:38:10 visibility=2 useScore=2 cgGrades=50 cgColour1=white cgColour2=yellow cgColour3=red height=30\n"
	FHW 	= open(OUT, "w")
	FHW.write(header)
	for I in H.values():
		rvs 	= I.rvs
		mus 	= ",".join([str(rv.mu) for rv in rvs if hasattr(rv, "mu")])
		sigmas 	= ",".join([str(rv.sigma) for rv in rvs if hasattr(rv, "sigma")])
		weights = ",".join([str(rv.w) for rv in rvs if hasattr(rv, "w")])
		params 	= mus + "_" + sigmas + "_" + weights
		center 	= min(I.X)
		centered= [x-center for x in I.X]
		predicts= [(predict(rvs, x),i) for i, x in enumerate(centered)]
		prev 	= ""
		start 	= None
		for TYPE, i in predicts:
			if TYPE!=prev:
				if start is not None:
					if prev.type == "normal":
						ID 		= "initiation"
						score 	= "500"
						RGB 	= "255,0,0"
					else:
						ID 		= "elongation"
						RGB 	= "0,0,255"
					FHW.write(I.chrom+"\t" + str(I.X[start]) + "\t" + str(I.X[i]) + "\t" + ID + "\t" + 
					str(score) + "\t" + strand+ "\t" + str(I.X[start]) + "\t" + str(I.X[i]) + "\t" + RGB + "\t" + params + "\n")
				start 	= i
			prev=TYPE
			prevI=i
		if prev.type == "normal":
			ID 		= "initiation"
			score 	= "500"
			RGB 	= "255,0,0"
		else:
			ID 		= "elongation"
			RGB 	= "0,0,255"
			score 	= "100"
		FHW.write(I.chrom+"\t" + str(I.X[start]) + "\t" + str(I.X[prevI]) + "\t" + ID + "\t" + 
		str(score) + "\t" + strand+ "\t" + str(I.X[start]) + "\t" + str(I.X[prevI]) + "\t" + RGB + "\t" + params+ "\n")
	FHW.close()
