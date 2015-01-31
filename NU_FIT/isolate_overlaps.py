import time
def run(LST):
	new 	= list()
	LST.sort()
	while LST:
		i 			= 0
		N 			= len(LST)
		o_st,o_sp 	= LST[i][0],LST[i][1]
		while i < N and (LST[i][0] <= o_sp and LST[i][1]>= o_st) : #want to find where there are no overlaps split on that
			o_st, o_sp 	= min((o_st, LST[i][0])),max((o_sp, LST[i][1]))
			i+=1
		left, right 	=  LST[:i], LST[i:]
		if len(left)==1:
			new.append(left[0])
		LST 			= right
	return new

def merge(LST):
	new 	= list()
	LST.sort()
	while LST:
		i 			= 0
		N 			= len(LST)
		o_st,o_sp 	= LST[i][0],LST[i][1]
		while i < N and (LST[i][0] <= o_sp) : #want to find where there are no overlaps split on that
			o_st, o_sp 	= min((o_st, LST[i][0])),max((o_sp, LST[i][1]))
			i+=1

		left, right 	=  LST[:i], LST[i:]

		allName 		= ""
		for start, stop, name in left:
			allName+=str(start) + "_" + str(stop) + "_" + name + ","
		new.append((o_st, o_sp, allName))
		prev 			= 0

		LST 			= right
		LST.sort()
	return new




	
