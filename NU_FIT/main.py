import utils,read, fit_NU,write
def run(argv):
	H 		= {"convertAnnotation", "classify"}
	assert len(argv) > 1, 'need to specify program parameters'
	module 	= argv[1]
	assert module in H, 'need to specify either convertAnnotation or classify'
	if module == "convertAnnotation":
		D 			= utils.userParameters(argv)
		annotFile 	= D["-i"][0]
		s,c,st,sp,n 	= D["-j"]
		out 		= D["-o"][0]		


		utils.makeGeneFormat(annotFile, out,int(s), int(c), int(st), int(sp),int(n))
	if module == "classify":
		D 			= utils.userParameters(argv)
		regionFile 	= D["-j"][0]
		strand 		= D["-s"][0]
		BedGraphFile= D["-i"][0]
		OUT 		= D["-o"][0]
		test 		= False
		if "-t" in D:
			test 	= True
		regions 	= read.readIntervals(regionFile)#read in annotation intervals
		
		H 			= read.insertBedGraphFile(BedGraphFile, regions, strand, test=test) #assert coverage data into each of the annotations
		H 			= fit_NU.run(H)
		write.writeIGV(H, OUT,strand)





