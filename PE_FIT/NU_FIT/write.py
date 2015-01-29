import numpy as np
def predict(rvs, x):
	norm 	= sum([rv.pdf(x) for rv in rvs])
	probs 	= np.array([rv.pdf(x) / norm for rv in rvs])
	return rvs[probs.argmax()]



def writeIGV(H, OUT,strand):
	header 	= "track name=Pausing_Elongation_Calls 2015-01-08 17:38:10 visibility=2 useScore=2 cgGrades=50 cgColour1=white cgColour2=yellow cgColour3=red height=30\n"
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
						ID 		= "pausing"
						score 	= "500"
						RGB 	= "255,0,0"
					else:
						ID 		= "elongation"
						RGB 	= "0,0,255"
					FHW.write(I.chrom+"\t" + str(I.X[start]) + "\t" + str(I.X[i]) + "\t" + ID + "\t" + 
					str(score) + "\t" + strand+ "\t" + str(I.X[start]) + "\t" + str(I.X[i]) + "\t" + RGB + "\n")
				start 	= i
			prev=TYPE
			prevI=i
		if prev.type == "normal":
			ID 		= "pausing"
			score 	= "500"
			RGB 	= "255,0,0"
		else:
			ID 		= "elongation"
			RGB 	= "0,0,255"
			score 	= "100"
		FHW.write(I.chrom+"\t" + str(I.X[start]) + "\t" + str(I.X[prevI]) + "\t" + ID + "\t" + 
		str(score) + "\t" + strand+ "\t" + str(I.X[start]) + "\t" + str(I.X[prevI]) + "\t" + RGB + "\n")
	FHW.close()
