import utils,read, fit_NU,write, os,sys
def run(argv):
	H 		= {"convertAnnotation", "classify", "concatenate"}
	assert len(argv) > 1, 'need to specify program parameters'
	module 	= argv[1]
	assert module in H, 'need to specify either convertAnnotation or classify'
	if module == "convertAnnotation":
		D 			= utils.userParameters(argv)
		annotFile 	= D["-i"][0]
		s,c,st,sp,n = D["-j"]
		out 		= D["-o"][0]		
		verbose 	= D["-v"]


		utils.makeGeneFormat(annotFile, out,int(s), int(c), int(st), int(sp),int(n))
	if module == "classify":
		D 			= utils.userParameters(argv)
		verbose 	= D["-v"]
		if D is None:
			print "exiting..."
			return False 
		regionFile 	= D["-j"]
		if regionFile is None or not os.path.isfile(regionFile[0]):
			print "(-j) not found or file does not exist."
			print "exiting..."
			return False
		regionFile 	= regionFile[0]
		if D["-s"]:
			strand 		= D["-s"][0]
		else:
			strand 		= "+"
		if strand is None:
			print "warning, user did not specify strand assuming in forward strand in the IGV output file"
			strand 	= "+"
		BedGraphFile= D["-i"]
		if BedGraphFile is None or not os.path.isfile(BedGraphFile[0]):
			print "(-i) not found or file does not exist"
			print "exiting..."
			return False
		BedGraphFile= D["-i"][0]
		maxBIC, penality 		= None, None
		rt 						= 1
		if D["-BIC"] is not None:
			if len(D["-BIC"])!= 2:
				print "-BIC command found, but not the right number of parameters (int int)"
				print "exiting..."
				return False
			maxBIC, penality 	= D["-BIC"]
			maxBIC, penality 	= int(maxBIC), float(penality)
		if D["-rt"] is not None:
			if len(D["-rt"]!= 1):
				print "-rt option found, but not the right number of parameters (int)"
				print "exiting..."
				return False
			rt 	= int(D["-rt"])
		
		OUT 		= D["-o"]
		if OUT is None:
			print "(-o) please specify an OUT file name and path"
			print "exiting"
			return False
		OUT 		= D["-o"][0]
		if verbose:
			utils.WELCOME(D)
		if D["-chr"]:
			specChrom 	= D["-chr"][0]
		else:
			specChrom 	= None
		test 		= D["-t"]
		single 		= D["-single"]
		if D["-np"]:
			np 			= int(D["-np"][0])
		else:
			np 			= 8

		if test:
			print "...warning, running test"
		regions 	= read.readIntervals(regionFile,single=single)#read in annotation intervals
		if verbose:
			sys.stdout.flush()
			print "reading in BedGraphFile      :",
			sys.stdout.flush()
		H 			= read.insertBedGraphFile(BedGraphFile, regions, strand, test=test, spec=specChrom) #assert coverage data into each of the annotations
		if verbose:
			print "finished"
			sys.stdout.flush()
			print "running mixture model        :",
		H 			= fit_NU.run(H, np=np,maxBIC=maxBIC, penality=penality,rt=rt)
		if verbose:
			print "finished"
		D["-BIC"] 	= maxBIC, penality
		D["-rt"] 	= rt
		write.writeIGV(H, OUT,strand,D)
	if module == "concatenate":
		pass





