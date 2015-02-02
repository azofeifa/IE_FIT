import utils,read, fit_NU,write, os,sys,time
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
		if D is None:
			print "exiting..."
			return False 
		
		verbose 	= D["-v"]
		regionFile 	= D["-j"]
		if regionFile is None or not os.path.isfile(regionFile[0]):
			print "(-j) not found or file does not exist."
			print "exiting..."
			return False
		regionFile 	= regionFile[0]
		if D["-s"] is not None:
			if len(D["-s"]) != 1:
				print "-s command found but not the right number of parameters (int)"
				print "exiting..."
				return False
			else:
				strand 		= D["-s"][0]
		else:
			print "-s command not found, please specify strand"
			print "exiting..."
			return False
		if strand is None:
			print "warning, user did not specify strand assuming in forward strand in the IGV output file"
			strand 	= "+"
		BedGraphFile= D["-i"]
		if BedGraphFile is None or not os.path.isfile(BedGraphFile[0]):
			print "(-i) not found or file does not exist"
			print "exiting..."
			return False
		BedGraphFile 			= D["-i"][0]
		maxBIC, penality 		= None, None
		binSize 				= 200
		rt 						= 1
		interval 				= None
		if D["-int"] is not None:
			if len(D["-int"]) != 1:
				print "-int command found, but not the right number of parameters (int)"
				print "exiting"
				return False
			interval 	= int(D["-int"][0])
		if D["-BIC"] is not None:
			if len(D["-BIC"])!= 2:
				print "-BIC command found, but not the right number of parameters (int int)"
				print "exiting..."
				return False
			maxBIC, penality 	= D["-BIC"]
			maxBIC, penality 	= int(maxBIC), float(penality)
		if D["-rt"] is not None:
			if len(D["-rt"])!= 1:
				print "-rt option found, but not the right number of parameters (int)"
				print "exiting..."
				return False
			rt 	= int(D["-rt"][0])
		
		OUT 		= D["-o"]
		if OUT is None:
			print "(-o) please specify an OUT file name and path"
			print "exiting"
			return False
		OUT 		= D["-o"][0]
		if D["-bin"] is not None:
			if len(D["-bin"])!=1:
				print "-bin option found, but not the right number of parameters (int)"
				print "exiting"
				return False
			binSize = int(D["-bin"][0])
		if verbose:
			utils.WELCOME(D)
		if D["-chr"]:
			specChrom 	= D["-chr"][0]
		else:
			specChrom 	= None
		test 		= D["-t"]
		single 		= D["-single"]
		merge 		= D["-merge"]
		if single and merge:
			print "both -single and -merge options are found however they are mutually exclusive, choose one"
			print "exiting..."
			return False
		if D["-np"]:
			np 			= int(D["-np"][0])
		else:
			np 			= 8

		if test:
			print "...warning, running test"
		regions 	= read.readIntervals(regionFile,single=single, merge=merge, interval=interval)#read in annotation intervals
		if verbose:
			sys.stdout.flush()
			print "reading in BedGraphFile      :",
			sys.stdout.flush()
		H 			= read.insertBedGraphFile(BedGraphFile, regions, strand, test=test, spec=specChrom) #assert coverage data into each of the annotations
		if verbose:
			print "finished"
			sys.stdout.flush()
			print "running mixture model        :",
		start 		= time.clock()
		H 			= fit_NU.run(H, np=np,maxBIC=maxBIC, penality=penality,rt=rt,binSize=binSize,strand=strand)
		if verbose:
			print "finished"
		D["-BIC"] 	= maxBIC, penality
		D["-rt"] 	= rt
		D["-chr"] 	= specChrom
		D["-bin"] 	= binSize
		D["-time"] 	= time.clock()-start
		write.writeIGV(H, OUT,strand,D)
	if module == "concatenate":
		D 			= utils.userParameters(argv)
		DIR,OUT 	= None, None
		if D is None:
			print "(-i) and (-j) commands not found... "
			print "exiting..."
			return False 
		if D["-i"] is None:
			print "-i command not found, need to specify an input directory (path)"
			print "exiting"
			return False
		else:
			if len(D["-i"]) != 1:
				print "-i command found but not the write number of parameters..."
				print "exiting..."
				return False
			else:
				DIR=D["-i"[0]]
		if D["-j"] is None:
			print "-j command not found, need to specify an output file name (path)"
			print "exiting"
			return False
		else:
			if len(D["-j"]) != 1:
				print "-j command found but not the write number of parameters..."
				print "exiting..."
				return False
			else:
				OUT=D["-j"[0]]
		read.readDirIE_OUT(DIR, OUT)







		





