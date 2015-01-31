import numpy as np
import datetime
def predict(rvs, x):
	norm 	= sum([rv.pdf(x) for rv in rvs])
	probs 	= np.array([rv.pdf(x) / norm for rv in rvs])
	return rvs[probs.argmax()]



def writeIGV(H, OUT,strand,D):
	ID 		= "IE_" + str(D["-chr"]) +"_" +str(",".join([str(i) for i in D["-BIC"]])) + "_"+str(D["-rt"]) + "_" + str(D["-bin"]) + "_" + str(D["-time"])
	DATE 	= str(datetime.datetime.now())
	header 	= "track name=" + ID+ " " + DATE  + " visibility=2 useScore=2 cgGrades=50 cgColour1=white cgColour2=yellow cgColour3=red height=30\n"
	FHW 	= open(OUT, "w")
	FHW.write(header)
	for I in H.values():
		rvs 	= I.rvs
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
						RGB 	= "255,0,0"
						score 	= str(prev.w)
						params 		= str(prev.mu+center) + "," + str(prev.sigma) + "," + str(prev.w)+":"+I.name
					else:
						ID 		= "elongation"
						RGB 	= "0,0,255"
						score 	= str(prev.w)
						params 		= str(prev.a+center) + "," + str(prev.b+center) + "," + str(prev.w) + ":"+I.name
					FHW.write(I.chrom+"\t" + str(I.X[start]) + "\t" + str(I.X[i]) + "\t" + ID + "\t" + 
					str(score) + "\t" + strand+ "\t" + str(I.X[start]) + "\t" + str(I.X[i]) + "\t" + RGB + "\t" + params + "\n")
				start 	= i
			prev=TYPE
			prevI=i
		if prev.type == "normal":
			ID 		= "initiation"
			score 	= str(prev.w)
			RGB 	= "255,0,0"
			params 		= str(prev.mu+center) + "," + str(prev.sigma) + "," + str(prev.w)+ ":"+I.name
		else:
			ID 		= "elongation"
			RGB 	= "0,0,255"
			score 	= str(prev.w)
			params 		= str(prev.a+center) + "," + str(prev.b+center) + "," + str(prev.w)+ ":"+I.name
	
		FHW.write(I.chrom+"\t" + str(I.X[start]) + "\t" + str(I.X[prevI]) + "\t" + ID + "\t" + 
		str(score) + "\t" + strand+ "\t" + str(I.X[start]) + "\t" + str(I.X[prevI]) + "\t" + RGB + "\t" + params+ "\n")
	
	FHW.close()



