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
		regionFile 	= D["-j"][0]
		if regionFile is None or not os.path.isfile(regionFile):
			print "(-j) not found or file does not exist."
			print "exiting..."
			return False
		if D["-s"]:
			strand 		= D["-s"][0]
		else:
			strand 		= "+"
		if strand is None:
			print "warning, user did not specify strand assuming in forward strand in the IGV output file"
			strand 	= "+"
		BedGraphFile= D["-i"][0]
		if BedGraphFile is None or not os.path.isfile(BedGraphFile):
			print "(-i) not found or file does not exist"
			print "exiting..."
			return False
		OUT 		= D["-o"][0]
		if OUT is None:
			print "(-o) please specify an OUT file name and path"
			print "exiting"
			return False
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
			print "reading in BedGraphFile    :",
			sys.stdout.flush()
		H 			= read.insertBedGraphFile(BedGraphFile, regions, strand, test=test, spec=specChrom) #assert coverage data into each of the annotations
		if verbose:
			print "finished"
			sys.stdout.flush()
			print "running mixture model      :",
		H 			= fit_NU.run(H, np=np)
		if verbose:
			print "finished"
		write.writeIGV(H, OUT,strand)
	if module == "concatenate":
		pass





