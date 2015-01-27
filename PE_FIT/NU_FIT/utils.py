import time
def makeGeneFormat(FILE, OUT,s, c, st, sp,n):#need to have position of strand, chrom,start, stop of genes
	FH 					= open(FILE)
	D 					= {"+":{}, "-":{}}
	header 				= True
	for line in FH:
		if not header:
			lineArray 		= line.strip("\n").split("\t")
			strand, chrom 	= lineArray[s], lineArray[c]
			start, stop 	= int(lineArray[st]), int(lineArray[sp])
			name 			= lineArray[n]
			if chrom not in D[strand]:
				D[strand][chrom] 	= list()
			D[strand][chrom].append((start, stop,name))
		else:
			header 			= False
	FH.close()
	#sort
	for strand in D:
		for chrom in D[strand]:
			D[strand][chrom].sort()
	#write out to FILE
	FHW 	= open(OUT, "w")
	for strand in D:
		for chrom in D[strand]:
			for start, stop,name in D[strand][chrom]:
				FHW.write(str(strand)+ "\t" + str(chrom) + "\t" + str(start) + "\t" + str(stop) + "\t" + name+"\n" )
	FHW.close()
def userParameters(argv):
	h 	= None
	D 	= {}
	for a in argv:
		if "-"==a[0]:
			h 	= a
			D[h] = list()
		elif h:
			D[h].append(a)
	return D
			
class info:
	def __init__(self, start, stop, info=None):
		self.start 	= start
		self.stop 	= stop
		self.info 	= info

class node:
	def __init__(self,start, stop):
		self.start 				= start
		self.stop 				= stop
		self.intervals 			= list()
		self.right, self.left 	= None, None
	def search(self, interval):
		finds 					= list()
		for info in self.intervals:
			if info[0] < interval[0] <info[1] or info[0] < interval[1] <info[1] or (interval[0] < info[0] and interval[1] > info[1]) :
				finds.append(info)
		return finds
	def checkIntervals(self):
		d 		= list()
		while self.intervals:
			curr 	= self.intervals.pop()
			ADD 	= True
			if not d:
				d.append(curr)
			else:
				for i,prev in enumerate(d):
					start, stop,name 	= prev
					if curr[0]==prev[0] and curr[1]==prev[1]:
						name+="," +curr[2]
						ADD 	= False
					d[i] 	= ( start, stop, name)
			if ADD:
				d.append(curr)
		self.intervals=d


	def __str__(self):
		return str(self.start)+"-"+str(self.stop) + ", " + str(len(self.intervals))
class treeNode:
	def __init__(self):
		self.node 	= None
		self.right 	= None
		self.left 	= None
	def build(self, nodes):
		i 			= len(nodes)/ 2
		self.node 	= nodes[i]

		if i > 0:
			self.left 	= treeNode()
			self.left.build(nodes[:i])
		if i+1 < len(nodes):
			self.right 	= treeNode()
			self.right.build(nodes[i+1:])
	def searchInterval(self, interval):
		start, stop 	= interval
		if self.node.start < start < self.node.stop or self.node.start < stop < self.node.stop or (start < self.node.start and stop > self.node.stop):
			return self.node.search(interval)
		if stop < self.node.start and self.left:
			return self.left.searchInterval(interval)
		if start > self.node.stop and self .right:
			return self.right.searchInterval(interval)
		return None

	def searchPoint(self, point):
		raise TypeError, "searchPoint method not completed yet"


class tree:
	def __init__(self, *args):
		assert len(args) < 2, "either no or 1 argumnet"
		self.root 	 	= None
		if len(args)==1 and args[0]:
			self.build(args[0])
	def assemble(self, LST):
		nodes 	= list()
		while LST:
			i 	= 0
			N 	= len(LST)
			o_st,o_sp 	= LST[i][0],LST[i][1]
			while i < N and (LST[i][0] <= o_sp and LST[i][1]>= o_st) : #want to find where there are no overlaps split on that
				o_st, o_sp 	= min((o_st, LST[i][0])),max((o_sp, LST[i][1]))
				i+=1
			left, right 	=  LST[:i], LST[i:]
			x,y 			= [info[0] for info in left],[info[1] for info in left]
			NODE 			= node(min(x), max(y))
			NODE.intervals 	= left	
			NODE.checkIntervals()	

			nodes.append(NODE)
			LST 			= right
		return nodes
			
	def build(self, LST):
		LST.sort()
		nodes 	= self.assemble(LST)
		#root is the middle of nodes
		i 			= len(nodes)/ 2
		self.root 	= treeNode()
		self.root.build(nodes)
	def searchInterval(self, interval):
		assert self.root, "interval tree has not been built yet"
		return self.root.searchInterval(interval)




