import utils,time,isolate_overlaps,os
def readIntervals(FILE,STRAND, single=False, merge=False, interval=None,pad=(0,0)):
	FH 		= open(FILE)
	D 		= {"+":{}, "-":{}}
	lines	= FH.readlines()

	for i,line in enumerate(lines):
		lineArray 	= line.strip("\n").split("\t")
		assert len(lineArray)==5, "strand\tchrom\tstart\t\stop\tname(ID)\n format please..."
		strand, chrom,start, stop,name 	= lineArray
		if strand not in D:
			D[strand] 	= {}
		if chrom not in D[strand]:
			D[strand][chrom] 		= list()
		D[strand][chrom].append((int(start), int(stop), name))

	#sort
	for strand in D:
		for chrom in D[strand]:
			N 						= len(D[strand][chrom])
			size 					= N /50.
			D[strand][chrom].sort()
			if interval is not None:
				start, stop 			= size*interval,size*(interval+1)
				D[strand][chrom] 		= [d for i,d in enumerate(D[strand][chrom]) if start<=i<=stop]
			if single:
				D[strand][chrom] 	= isolate_overlaps.run(D[strand][chrom])
			elif merge:
				D[strand][chrom] 	= isolate_overlaps.merge(D[strand][chrom])
			else:
				D[strand][chrom] 	= D[strand][chrom]
			if STRAND=="-":
				D[strand][chrom] 		= [(start-pad[1], stop+pad[0], name) for start, stop, name in D[strand][chrom]]
			else:
				D[strand][chrom] 		= [(start-pad[0], stop+pad[1], name) for start, stop, name in D[strand][chrom] ]
			D[strand][chrom] 			= utils.tree(D[strand][chrom])
				
		if STRAND not in D:
			print "user specified strand is not present in annotation file"
			D 	= None 
	return D
class interval:
	def __init__(self, start, stop, name,chrom):
		self.start, self.stop,self.name,self.chrom 	= start, stop, name,chrom
		self.X, self.Y 								= list(),list()
		self.rvs 									= None
	def __str__(self):
		return self.chrom+":"+str(self.start) + "-"+str(self.stop) + ", " + self.name

def insertBedGraphFile(bgFile, D,strand,test=False, spec=None):
	FH 	= open(bgFile)
	D 	= D[strand]
	prev= ""
	H 	= {}
	i 	= 0
	for line in FH:
		lineArray 	= line.strip('\n').split("\t")
		i+=1
		assert len(lineArray)==4, "error in bedgraph file line: + " +str(i) +  ", format: chrom\tstart\tstop\tcoverage\n"
		chrom,start, stop, cov 	= lineArray

		if spec is None or chrom==spec:
			start, stop, cov 		= int(start), int(stop), int(cov)
			F 						= D[chrom].searchInterval((start, stop))
			if (test and start > pow(10,7)) or (chrom != spec and spec is not None):
				break
			if F:
				for st, sp, name in F:
					if name not in H:
						H[name] 	= interval(st, sp, name, chrom)
					H[name].X+=[i for i in range(start, stop)]
					H[name].Y+=[cov for i in range(start, stop)]
	FH.close()
	return H

def readDirIE_OUT(DIR, out):
	FHW 	= open(out, "w")
	for f in os.listdir(DIR):
		FH 	= open(DIR+f)
		for line in FH:
			FHW.write(line)
		FH.close()
	FHW.close()




















