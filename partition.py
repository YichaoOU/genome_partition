import numpy as np

##-------------- below is k-partite graph algorithm ----------------##

class node:

	def __init__(self,name):
		self.name = name
	def display(self):
		print "node name: ", self.name
	
class edge:

	def __init__(self,pre,next,weight):
		self.pre = pre
		self.next = next
		self.weight = weight
		
	def display(self):
		print self.pre,"->",self.next," weight: ", self.weight

class part:

	def __init__(self):
		self.node=[]
		self.edge=[]
		
	def insertNode(self,node):
		self.node.append(node)
		
	def insertEdge(self,edge):
		self.edge.append(edge)
		
	def display(self):
		for node in self.node:
			node.display()
		for edge in self.edge:
			edge.display()

class digraph:

	def __init__(self,num):
		self.part = [0]*num
		self.size = num
		
	def insertPart(self,part,num):
		self.part[num]=part
		
	def findWeight(self,arrowTo,num):
		output_weight = []
		output_name   = []
		for item in self.part[num].edge:
			if item.next == arrowTo:
				output_name.append(item.pre)
				output_weight.append(item.weight)
		return [output_name,output_weight]
		
	def display(self):
		for i in range(len(self.part)):
			print "part",i
			self.part[i].display()

def longestPath(g):
	max_dis = {}
	max_path = {}
	for node in g.part[1].node:
		max_dis[node.name]=0
		max_path[node.name] = ["-1->"+str(node.name)]
	for i in range(2,g.size):
		for node in g.part[i].node:
			temp_dis={}
			[name,weight] = g.findWeight(node.name,i)
			for j in range(len(name)):
				temp_dis[name[j]] = max_dis[name[j]]+weight[j]
			max_dis[node.name] = max(temp_dis.values())
			best_path = []
			for key in temp_dis.keys():
				if temp_dis[key] == max_dis[node.name]:
					best_path.append(str(key)+"->"+str(node.name))
			max_path[node.name] = best_path	
	#print "in loop"
	#print max_path
	return [max_dis,max_path]

def findPath(path,node):
	while(path.has_key(node)):
		print node
		node = int(path[node][0].split("->")[0])

##-------------- above is k-partite graph algorithm ----------------##	
			
class block:

	def __init__(self,start,end):
		self.start 	= start
		self.end	= end
		
	def name(self,name):
		self.name   = name
		
	def display(self):
		print "start: ",self.start," end: ",self.end, " name: ", self.name
		
##------------------- below is compatible test ---------------------##			
def FGT(matrix):
	a =matrix.tolist()
	b=[]
	for item in a:
		b.append(str(item))
	c = list(set(b))	
	if len(c) == 4:
		return True
	else:
		return False

def Compatible(matrix):
	[m,n] = matrix.shape # n=2
	gamete = [0]*9
	for i in range(9):
		gamete[i] = []
	gamete01 = 0
	gamete10 = 0
	gamete11 = 0
	gamete00 = 0
	# 2 : grey color, means the pair must be compatible
	# 0 : green color, must be in-phase
	# 1 : orange color, must be out-phase
	# 3 : blue color, either in or out phase
	# 4 : red color, must be incompatible
	for i in range(m):
		first = int(matrix[i,0])
		second = int(matrix[i,1])
		index = first*3+second
		gamete[index].append(i)
		if first == 0:
			if second == 1:
				gamete01 = 1
			elif second == 0:
				gamete00 = 1
			else:
				gamete00 = 1
				gamete01 = 1
		elif first == 1:
			if second == 1:
				gamete11 = 1
			elif second == 0:
				gamete10 = 1
			else:
				gamete10 = 1
				gamete11 = 1
		elif first == 2:
			if second == 1:
				gamete11 = 1
				gamete01 = 1
			if second == 0:
				gamete10 = 1
				gamete00 = 1	
	set1 = [0]*8
	for i in range(8):
		set1[i] = set(gamete[i])
	if gamete00 + gamete01 + gamete10 + gamete11 == 4:
		return [gamete[8],4]# 
	if len(gamete[8]) == 0:
		return [[],2]#
	if set1[4] | set1[7] | set1[5] != set():
		return [gamete[8],0]# 
	if set1[1] | set1[2] != set() and set1[3] | set1[6] != set():
		return [gamete[8],1]# 
	else:
		return [gamete[8],3]
	
def checkState(a,b):
	sum = a+b
	if sum == 2: # a=1,b=1
		return 0
	elif sum > 2: # one is 0/1 one is 3
		return min(a,b) # 0/1
	elif sum < 0: # one is -5
		return max(a,b) # 0/1/3
	else: # 
		return sum # 
		
##------------------- above is compatible test ---------------------##	
		
##----------- below is the haplotype & genotype class --------------##		
class haplotype:

	def __init__(self,file):
		self.file = file
		self.matrix = None
		self.core = []
		self.CLR = [] # store the blocks retrieved from L to R
		self.CRL = [] # store the blocks retrieved from R to L
		self.CUber = [] # store all the maximal interval
		self.kGroup = [] # store the group divided by core
		self.diGraph = None # store the k-partite structure
		self.best_path = None # store the maximal-k-cover set of blocks
		self.max_dis = None 
		self.flaggingSNP = []
		self.criticalSNP = []
		
	def read_file(self):
		array=[]
		for line in open(self.file).readlines():
			col=[]
			line=line.split()
			for item in line:
				col.append(item)
			array.append(col)
		self.matrix = np.array(array)
		[self.n,self.m] = self.matrix.shape

	def LRScan(self):
		current_start_site = 1
		l=[] # is the current interval
		# l will grow like this : l=[0] -> l=[0,1] -> l=[0,1,2] till the current interval is closed
		for i in range(self.m): # start from the leftmost SNP
			for j in range(len(l)): # perform FGT against every SNP in the current active interval
				if FGT(np.concatenate((self.matrix[:,l[j]:l[j]+1],self.matrix[:,i:i+1]),axis=1)):
					self.CLR.append(block(current_start_site,i)) # because the SNP(i+1) is incompatible, so SNP(i+1) is not belong to that interval, so the end position is i
					current_start_site=i+1 # i is the index of array, its real position is i+1
					l=[]
					break
			l.append(i)
		self.CLR.append(block(current_start_site,self.m))
	
	def displayCLR(self):
		for item in self.CLR:
			print "start: ",item.start," end: ",item.end

	def RLScan(self):
		current_start_site = self.m
		l=[] # is the current interval
		for i in range(self.m-1,-1,-1): # start from the rightmost SNP
			for j in range(len(l)): # perform FGT against every SNP in the current active interval
				if FGT(np.concatenate((self.matrix[:,l[j]:l[j]+1],self.matrix[:,i:i+1]),axis=1)):
					self.CRL.append(block(current_start_site,i+2)) # position(i+1) is not compatible, so it ends at position(i+2)
					current_start_site=i+1 # i is the index of array, its real position is i+1
					l=[]
					break
			l.append(i)
		self.CRL.append(block(current_start_site,1))

	def displayCRL(self):
		for item in self.CRL:
			print "start: ",item.start," end: ",item.end

	def UberScan(self):
		s=1 # start position
		l=[]
		flagSNP = []
		count_interval_name = 1
		for i in range(self.m): # from the leftmost SNP
			for j in range(len(l)-1,-1,-1): # to get the nearst SNP sj that is incompatible with i+1
				#print "i: ",i," j: ",j," l[j]: ",l[j]
				if FGT(np.concatenate((self.matrix[:,l[j]:l[j]+1],self.matrix[:,i:i+1]),axis=1)):
					temp = block(s,i)
					flagSNP.append(s-1)
					flagSNP.append(i+1)
					temp.name = count_interval_name
					count_interval_name += 1
					self.CUber.append(temp) # position(i+1) is incompatible, so the end position is i
					s=l[j]+2 # l[j] + 1 is the incompatible SNP's position, so start site should be l[j] +2
					new=[] # retrieve the position is the rest of l (previous active interval)
					for x in range(j+1,len(l)): # retrieve with the original order
						new.append(l[x])
					l=new
					break
			l.append(i)
		temp = block(s,self.m)
		flagSNP.append(s-1)
		flagSNP.remove(0)
		self.flaggingSNP = list(set(flagSNP))
		temp.name = count_interval_name
		self.CUber.append(temp)

	def displayCUber(self):
		for item in self.CUber:
			item.display()

	def getCore(self):
		s=1
		core_name = 1
		for i in range(len(self.CLR)):# follow the definition of core
			j = len(self.CLR)-i-1
			LR_start = self.CLR[i].start
			RL_end = self.CRL[j].start
			#print "LR_start: ", LR_start, " RL_end ",RL_end
			temp = block(LR_start,RL_end)
			temp.name(core_name)
			self.core.append(temp) # this is the fomula of core range
			core_name += 1
			
	def displayCore(self):
		for item in self.core:
			item.display()
	
	def assignKGroup(self):
		self.kGroup = range(len(self.core))
		for core in self.core:
			self.kGroup[core.name-1] = []
			for interval in self.CUber:
				if (interval.start <= core.start) and (interval.end >= core.end):
					self.kGroup[core.name-1].append(interval)
		
	def displayKGroup(self):
		for i in range(len(self.kGroup)):
			print "group: ",i+1
			for interval in self.kGroup[i]:
				interval.display()
	
	def maximalKCover(self):
		self.diGraph = digraph(len(self.core)+2)
		#part0
		source = part()
		source.insertNode(node(-1))
		self.diGraph.insertPart(source,0)
		sink = part()
		sink.insertNode(node(-2))
		self.diGraph.insertPart(sink,len(self.core)+1)
		# insert part1
		part1=part()
		for interval in self.kGroup[0]:
			part1.insertNode(node(interval.name))
			part1.insertEdge(edge(-1,interval.name,0))
		self.diGraph.insertPart(part1,1)
		# insert last part
		last = part()
		for interval in self.kGroup[len(self.core)-1]:
			last.insertNode(node(-2))
			last.insertEdge(edge(interval.name,-2,0))
		self.diGraph.insertPart(last,len(self.core)+1)
		# insert other part
		for i in range(1,len(self.kGroup)):
			current_part = part()
			#print "self.kGroup[i] ",self.kGroup[i]
			for interval in self.kGroup[i]:
				#print "part",i+1,": interval name: ", interval.name
				current_part.insertNode(node(interval.name))
				for pre in self.kGroup[i-1]:
					weight = pre.end - interval.start + 1
					if weight >= 0:
						current_part.insertEdge(edge(pre.name,interval.name,weight))
			self.diGraph.insertPart(current_part,i+1)
		
	def displayPath(self):
		[self.max_dis,self.best_path] = longestPath(self.diGraph)
		findPath(self.best_path,-2)
		
	def displayDiGraph(self):
		self.diGraph.display()
		
	def findSusSNP(self):
		for pos in self.flaggingSNP:
			temp = haplotype(self.file)
			i = pos-1
			temp.matrix = np.delete(self.matrix,i,1)
			[temp.n,temp.m] = temp.matrix.shape
			temp.RLScan()
			temp.LRScan()
			temp.getCore()	
			if len(temp.core) < len(self.core):
				self.criticalSNP.append(pos)
					
class genotype:

	def __init__(self,file):
		self.file = file
		self.matrix = None
		# optimal compatible interval variables
		self.opticore = []
		self.optiCLR = [] # store the blocks retrieved from L to R
		self.optiCRL = [] # store the blocks retrieved from R to L
		self.optiCUber = [] # store all the maximal interval
		self.optikGroup = [] # store the group divided by core
		self.optidiGraph = None
		self.optibest_path = None
		self.optimax_dis = None
		self.optiflaggingSNP = []
		self.opticriticalSNP = []
		# pesitic compatible interval variables
		self.pesicore = []
		self.pesiCLR = [] # store the blocks retrieved from L to R
		self.pesiCRL = [] # store the blocks retrieved from R to L
		self.pesiCUber = [] # store all the maximal interval
		self.pesikGroup = [] # store the group divided by core
		self.pesidiGraph = None
		self.pesibest_path = None
		self.pesimax_dis = None
		self.pesiflaggingSNP = []
		self.pesicriticalSNP = []
	
	def read_file(self):
		array=[]
		for line in open(self.file).readlines():
			col=[]
			line=line.split()
			for item in line:
				col.append(item)
			array.append(col)
		self.matrix = np.array(array)
		[self.n,self.m] = self.matrix.shape

	def optiLRScan(self):	
		array = [-5]*self.n
		row = np.matrix(array).T
		in_out_matrix = row # initial
		breakSignal = 0 # initial
		current_start_site=1
		l=[]
		for i in range(self.m):
			in_out_matrix = np.append(in_out_matrix,row,1)
			for j in range(len(l)):
				[pair22,state] = Compatible(np.concatenate((self.matrix[:,l[j]:l[j]+1],self.matrix[:,i:i+1]),axis=1))
				if state == 2: # do nothing
					continue
				## below is most difficult..................................................
				for record in pair22:
					current_state = checkState(in_out_matrix[record,i],in_out_matrix[record,l[j]])
					if current_state == 3 or current_state == -5:
						if state == 1:
							in_out_matrix[record,i] = 1
							in_out_matrix[record,l[j]] = 0
						else:
							in_out_matrix[record,i] = state
							in_out_matrix[record,l[j]] = state
					elif in_out_matrix[record,i] == -5 or in_out_matrix[record,i] == 3: # 0/1
						if state == 3:
							in_out_matrix[record,i] = 3
						else:
							in_out_matrix[record,i] = abs(in_out_matrix[record,l[j]]-state)
					elif in_out_matrix[record,l[j]] == -5 or in_out_matrix[record,l[j]] == 3:
						if state == 3:
							in_out_matrix[record,l[j]] = 3
						else:
							in_out_matrix[record,l[j]] = abs(in_out_matrix[record,i]-state)
					else:
						if state != 3:
							if current_state != state:
								breakSignal = 1
				## above is most difficult..................................................	
				if state == 4 or breakSignal == 1: # break			
					self.optiCLR.append(block(current_start_site,i)) # because the SNP(i+1) is incompatible, so SNP(i+1) is not belong to that interval, so the end position is i
					current_start_site=i+1
					l=[]
					[x,y] = in_out_matrix.shape
					for ii in range(x):
						for jj in range(y):
							in_out_matrix[ii,jj] = -5
					breakSignal = 0 # initial
					break
			l.append(i)
		self.optiCLR.append(block(current_start_site,self.m))
		
	def pesiLRScan(self):	
		current_start_site = 1
		l=[] # is the current interval
		# l will grow like this : l=[0] -> l=[0,1] -> l=[0,1,2] till the current interval is closed
		for i in range(self.m): # start from the leftmost SNP
			for j in range(len(l)): # perform FGT against every SNP in the current active interval
				[pair22,state] = Compatible(np.concatenate((self.matrix[:,l[j]:l[j]+1],self.matrix[:,i:i+1]),axis=1))
				if state != 2:
					self.pesiCLR.append(block(current_start_site,i)) # because the SNP(i+1) is incompatible, so SNP(i+1) is not belong to that interval, so the end position is i
					current_start_site=i+1 # i is the index of array, its real position is i+1
					l=[]
					break
			l.append(i)
		self.pesiCLR.append(block(current_start_site,self.m))
		
	def pesiRLScan(self):
		current_start_site = self.m
		l=[] # is the current interval
		for i in range(self.m-1,-1,-1): # start from the rightmost SNP
			for j in range(len(l)): # perform FGT against every SNP in the current active interval
				[pair22,state] = Compatible(np.concatenate((self.matrix[:,l[j]:l[j]+1],self.matrix[:,i:i+1]),axis=1))
				if state != 2:
					self.pesiCRL.append(block(current_start_site,i+2)) # position(i+1) is not compatible, so it ends at position(i+2)
					current_start_site=i+1 # i is the index of array, its real position is i+1
					l=[]
					break
			l.append(i)
		self.pesiCRL.append(block(current_start_site,1))
		
	def optiRLScan(self):	
		array = [-5]*self.n
		row = np.matrix(array).T
		in_out_matrix = row # initial
		for i in range(self.m):
			in_out_matrix = np.append(in_out_matrix,row,1)
		breakSignal = 0 # initial
		current_start_site=self.m
		l=[]
		for i in range(self.m-1,-1,-1):
			for j in range(len(l)):
				[pair22,state] = Compatible(np.concatenate((self.matrix[:,l[j]:l[j]+1],self.matrix[:,i:i+1]),axis=1))
				if state == 2: # do nothing
					continue
				## below is most difficult..................................................
				for record in pair22:
					current_state = checkState(in_out_matrix[record,i],in_out_matrix[record,l[j]])
					if current_state == 3 or current_state == -5:
						if state == 1:
							in_out_matrix[record,i] = 1
							in_out_matrix[record,l[j]] = 0
						else:
							in_out_matrix[record,i] = state
							in_out_matrix[record,l[j]] = state
					elif in_out_matrix[record,i] == -5 or in_out_matrix[record,i] == 3: # 0/1
						if state == 3:
							in_out_matrix[record,i] = 3
						else:
							in_out_matrix[record,i] = abs(in_out_matrix[record,l[j]]-state)
					elif in_out_matrix[record,l[j]] == -5 or in_out_matrix[record,l[j]] == 3:
						if state == 3:
							in_out_matrix[record,l[j]] = 3
						else:
							in_out_matrix[record,l[j]] = abs(in_out_matrix[record,i]-state)
					else:
						if state != 3:
							if current_state != state:
								breakSignal = 1
				## above is most difficult..................................................	
				if state == 4 or breakSignal == 1: # break			
					self.optiCRL.append(block(current_start_site,i+2)) # because the SNP(i+1) is incompatible, so SNP(i+1) is not belong to that interval, so the end position is i
					current_start_site=i+1
					l=[]
					[x,y] = in_out_matrix.shape
					for ii in range(x):
						for jj in range(y):
							in_out_matrix[ii,jj] = -5
					breakSignal = 0 # initial
					break
			l.append(i)
		self.optiCRL.append(block(current_start_site,1))
	
	def optiUberScan(self):
		interval_name_count = 1 # give their name
		array = [-5]*self.n
		row = np.matrix(array).T
		in_out_matrix = row # initial
		breakSignal = 0 # initial
		current_start_site=1
		l=[]
		flagSNP=[]
		i = 0
		while i<self.m:
			in_out_matrix = np.append(in_out_matrix,row,1)
			for j in range(len(l)-1,-1,-1):
				[pair22,state] = Compatible(np.concatenate((self.matrix[:,l[j]:l[j]+1],self.matrix[:,i:i+1]),axis=1))
				if state == 2: # do nothing
					continue
				## below is most difficult..................................................
				for record in pair22:
					current_state = checkState(in_out_matrix[record,i],in_out_matrix[record,l[j]])
					if current_state == 3 or current_state == -5:
						if state == 1:
							in_out_matrix[record,i] = 1
							in_out_matrix[record,l[j]] = 0
						else:
							in_out_matrix[record,i] = state
							in_out_matrix[record,l[j]] = state
					elif in_out_matrix[record,i] == -5 or in_out_matrix[record,i] == 3: # 0/1
						if state == 3:
							in_out_matrix[record,i] = 3
						else:
							in_out_matrix[record,i] = abs(in_out_matrix[record,l[j]]-state)
					elif in_out_matrix[record,l[j]] == -5 or in_out_matrix[record,l[j]] == 3:
						if state == 3:
							in_out_matrix[record,l[j]] = 3
						else:
							in_out_matrix[record,l[j]] = abs(in_out_matrix[record,i]-state)
					else:
						if state != 3:
							if current_state != state:
								breakSignal = 1
				## above is most difficult..................................................	
				if state == 4 or breakSignal == 1: # break	
					temp = block(current_start_site,i)
					flagSNP.append(current_start_site-1)
					flagSNP.append(i+1)
					temp.name(interval_name_count)
					self.optiCUber.append(temp) # because the SNP(i+1) is incompatible, so SNP(i+1) is not belong to that interval, so the end position is i
					current_start_site=l[j]+2
					l=[]
					[x,y] = in_out_matrix.shape
					for ii in range(x):
						for jj in range(y):
							in_out_matrix[ii,jj] = -5
					breakSignal = 0 # initial
					i = current_start_site-1
					interval_name_count += 1
					break
			l.append(i)		
			i += 1
		temp = block(current_start_site,self.m)
		flagSNP.append(current_start_site-1)
		flagSNP.remove(0)
		self.optiflaggingSNP = list(set(flagSNP))
		temp.name(interval_name_count)
		self.optiCUber.append(temp)
	
	def pesiUberScan(self):
		s=1 # start position
		l=[]
		flagSNP = []
		count_interval_name = 1
		for i in range(self.m): # from the leftmost SNP
			for j in range(len(l)-1,-1,-1): # to get the nearst SNP sj that is incompatible with i+1
				[pair22,state] = Compatible(np.concatenate((self.matrix[:,l[j]:l[j]+1],self.matrix[:,i:i+1]),axis=1))
				if state != 2:
					temp = block(s,i)
					flagSNP.append(s-1)
					flagSNP.append(i+1)
					temp.name = count_interval_name
					count_interval_name += 1
					self.pesiCUber.append(temp) # position(i+1) is incompatible, so the end position is i
					s=l[j]+2 # l[j] + 1 is the incompatible SNP's position, so start site should be l[j] +2
					new=[] # retrieve the position is the rest of l (previous active interval)
					for x in range(j+1,len(l)): # retrieve with the original order
						new.append(l[x])
					l=new
					break
			l.append(i)
		temp = block(s,self.m)
		flagSNP.append(s-1)
		flagSNP.remove(0)
		self.pesiflaggingSNP = list(set(flagSNP))
		temp.name = count_interval_name
		self.pesiCUber.append(temp)
	
	def displayoptiCRL(self):
		for item in self.optiCRL:
			print "start: ",item.start," end: ",item.end
	
	def displaypesiCRL(self):
		print "pesiCRL"
		for item in self.pesiCRL:
			print "start: ",item.start," end: ",item.end

	def displayoptiCLR(self):
		for item in self.optiCLR:
			print "start: ",item.start," end: ",item.end

	def displaypesiCLR(self):
		print "pesiCLR"
		for item in self.pesiCLR:
			print "start: ",item.start," end: ",item.end
		
	def displayoptiCUber(self):
		for item in self.optiCUber:
			item.display()
	
	def displaypesiCUber(self):
		for item in self.pesiCUber:
			item.display()
	
	def getoptiCore(self):
		s=1
		core_name = 1
		for i in range(len(self.optiCLR)):# follow the definition of core
			j = len(self.optiCLR)-i-1
			LR_start = self.optiCLR[i].start
			RL_end = self.optiCRL[j].start
			temp = block(LR_start,RL_end)
			temp.name(core_name)
			self.opticore.append(temp) # this is the fomula of core range
			core_name += 1
			
	def getpesiCore(self):
		s=1
		core_name = 1
		for i in range(len(self.pesiCLR)):# follow the definition of core
			j = len(self.pesiCLR)-i-1
			LR_start = self.pesiCLR[i].start
			RL_end = self.pesiCRL[j].start
			temp = block(LR_start,RL_end)
			temp.name(core_name)
			self.pesicore.append(temp) # this is the fomula of core range
			core_name += 1
	
	def displayoptiCore(self):
		for item in self.opticore:
			item.display()
	
	def displaypesiCore(self):
		for item in self.pesicore:
			item.display()
	
	def assignpesiKGroup(self):
		self.pesikGroup = range(len(self.pesicore))
		for core in self.pesicore:
			self.pesikGroup[core.name-1] = []
			for interval in self.pesiCUber:
				if (interval.start <= core.start) and (interval.end >= core.end):
					self.pesikGroup[core.name-1].append(interval)

	def assignoptiKGroup(self):
		self.optikGroup = range(len(self.opticore))
		for core in self.opticore:
			self.optikGroup[core.name-1] = []
			for interval in self.optiCUber:
				if (interval.start <= core.start) and (interval.end >= core.end):
					self.optikGroup[core.name-1].append(interval)
					
	def displaypesiKGroup(self):
		#print "asd"
		for i in range(len(self.pesikGroup)):
			print "group: ",i+1
			for interval in self.pesikGroup[i]:
				interval.display()
	
	def displayoptiKGroup(self):
		for i in range(len(self.optikGroup)):
			print "group: ",i+1
			for interval in self.optikGroup[i]:
				interval.display()
	
	def maximaloptiKCover(self):
		self.optidiGraph = digraph(len(self.opticore)+2)
		#part0
		source = part()
		source.insertNode(node(-1))
		self.optidiGraph.insertPart(source,0)
		sink = part()
		sink.insertNode(node(-2))
		self.optidiGraph.insertPart(sink,len(self.opticore)+1)
		# insert part1
		part1=part()
		for interval in self.optikGroup[0]:
			part1.insertNode(node(interval.name))
			part1.insertEdge(edge(-1,interval.name,0))
		self.optidiGraph.insertPart(part1,1)
		# insert last part
		last = part()
		for interval in self.optikGroup[len(self.opticore)-1]:
			last.insertNode(node(-2))
			last.insertEdge(edge(interval.name,-2,0))
		self.optidiGraph.insertPart(last,len(self.opticore)+1)
		# insert other part
		for i in range(1,len(self.optikGroup)):
			current_part = part()
			for interval in self.optikGroup[i]:
				current_part.insertNode(node(interval.name))
				for pre in self.optikGroup[i-1]:
					weight = pre.end - interval.start + 1
					if weight >= 0:
						current_part.insertEdge(edge(pre.name,interval.name,weight))
			self.optidiGraph.insertPart(current_part,i+1)
		
	def maximalpesiKCover(self):
		self.pesidiGraph = digraph(len(self.pesicore)+2)
		#part0
		source = part()
		source.insertNode(node(-1))
		self.pesidiGraph.insertPart(source,0)
		sink = part()
		sink.insertNode(node(-2))
		self.pesidiGraph.insertPart(sink,len(self.pesicore)+1)
		# insert part1
		part1=part()
		for interval in self.pesikGroup[0]:
			part1.insertNode(node(interval.name))
			part1.insertEdge(edge(-1,interval.name,0))
		self.pesidiGraph.insertPart(part1,1)
		# insert last part
		last = part()
		for interval in self.pesikGroup[len(self.pesicore)-1]:
			last.insertNode(node(-2))
			last.insertEdge(edge(interval.name,-2,0))
		self.pesidiGraph.insertPart(last,len(self.pesicore)+1)
		# insert other part
		for i in range(1,len(self.pesikGroup)):
			current_part = part()
			for interval in self.pesikGroup[i]:
				current_part.insertNode(node(interval.name))
				for pre in self.pesikGroup[i-1]:
					weight = pre.end - interval.start + 1
					if weight >= 0:
						current_part.insertEdge(edge(pre.name,interval.name,weight))
			self.pesidiGraph.insertPart(current_part,i+1)
	
	def displayoptiPath(self):
		[self.optimax_dis,self.optibest_path] = longestPath(self.optidiGraph)
		findPath(self.optibest_path,-2)
		
	def displaypesiPath(self):
		[self.pesimax_dis,self.pesibest_path] = longestPath(self.pesidiGraph)
		findPath(self.pesibest_path,-2)
		
	def displayoptiDiGraph(self):
		self.optidiGraph.display()
		
	def displaypesiDiGraph(self):
		self.pesidiGraph.display()
		
	def findoptiSusSNP(self):
		for pos in self.optiflaggingSNP:
			temp = genotype(self.file)
			i = pos-1
			temp.matrix = np.delete(self.matrix,i,1)
			[temp.n,temp.m] = temp.matrix.shape
			temp.optiRLScan()
			temp.optiLRScan()
			temp.getoptiCore()
			if len(temp.opticore) < len(self.opticore):
				self.opticriticalSNP.append(pos)
			
	def findpesiSusSNP(self):
		for pos in self.pesiflaggingSNP:
			temp = genotype(self.file)
			i = pos-1
			temp.matrix = np.delete(self.matrix,i,1)
			[temp.n,temp.m] = temp.matrix.shape
			temp.pesiRLScan()
			temp.pesiLRScan()
			temp.getpesiCore()
			if len(temp.pesicore) < len(self.pesicore):
				self.pesicriticalSNP.append(pos)
				
