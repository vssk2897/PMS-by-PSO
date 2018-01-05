import random,itertools,time,math
from random import randint
# The Driver function to this module
def pso(sequences,l,d) :
	cp=list(sequences)
	print "\t\t----PSO  Initialising----"
	numIteration=200
	pop=30
	gBest=[0 for i in range(len(cp)+1)]
	# Randomly initialising the initial Population
	initialPopulation=matrixGenerate(cp,pop,l)
	X=list(initialPopulation)
	vel=initParticle(X)
	for count in range(numIteration) :
		for i in range(len(X)):
			for j in range(len(X[i])) :
				X[i][j]=int(math.ceil(X[i][j]))
		print "Iteration No :",count
		cp=list(sequences)
		# Refer paper provided in readme for the exact desciption of these constatnts
		c1=1.8
		c2=1.9
		w=0.5
		pBest=[0 for j in range(len(cp)+1)]
		for i in range(len(initialPopulation)) :
			pBest[0 : len(cp)]=list(X[i])
			children=childrenGenerate(cp,l,list(X[i]))
			seq=clever(children,d,l)
			possiblePairs=bestPairs(seq,d,l)
			if purity(possiblePairs)==True :
				comb=bestCombination(possiblePairs)
				locPop=list(initialPopulation[i])
				for it in range(len(comb)) :
					matrix=alignmentMatrix(comb[it])
					cons=consensus(matrix)
					sco=score(matrix)
					if sco > pBest[len(cp)] :
						pBest[len(cp)]=sco
						if it>=0 and it<=4 :
							for m in range(len(locPop)):
  								 pBest[m]= pBest[m]-(it+1)
  						if it>=5 and it<=9 :
							for m in range(len(locPop)):
  								 pBest[m]= pBest[m]+(it-5+1)		 
  						print "pBest",pBest		 		 
			if gBest[len(cp)] <= pBest[len(cp)] :
				gBest=list(pBest)
				print "gBest :",gBest	
			for j in range(len(X[i])):
				rand1=round(random.random(),3)
				rand2=round(random.random(),3)
				vel[i][j]= w*vel[i][j]+c1*rand1*(pBest[j]-X[i][j])+c2*rand2*(gBest[j]-X[i][j])
				X[i][j]=X[i][j]+vel[i][j]		
		
	#print '\n',gBest
	print "\t\t----PSO Terminated----"
	for i in range(len(X)):
			for j in range(len(X[i])) :
				X[i][j]=int(math.ceil(X[i][j]))	
	print X				
	return X				
def purity(possiblePairs) :
	for i in possiblePairs:
		for j in i :
			for k in j :
				if len(k)==0:
					return False
	return True							
	
# This function calculates the best combination from the posible pairs
def bestCombination(possiblePairs) :
	closePairs=list(possiblePairs)
	comb=[[[None] for j in range(len(possiblePairs[i])) ] for i in range(len(possiblePairs))]
	safe=[None for i in range(len(closePairs))]
	try :
		start=time.clock()
		while time.clock()-start <= 0.01 :
			for i in range(len(possiblePairs)):
				closePairs=list(possiblePairs)
				comb[i][len(closePairs[i])-1]=list(itertools.chain.from_iterable(closePairs[i][len(closePairs[i])-1]))
				comb[i][len(closePairs[i])-2]=list(itertools.chain.from_iterable(closePairs[i][len(closePairs[i])-2]))
				for j in xrange(len(closePairs[i])-3,-1,-1):
					count=0
					flag=False
					while flag==False and count<=10 :
						if len(closePairs[i][j]) ==0:
							count=count+1
							pass
						else :	
							k=random.randrange(0,len(closePairs[i][j]),1)
							if occurredBefore(comb,i,list(closePairs[i][j][k])) == False :
								comb[i][j]=list(closePairs[i][j][k])
								flag=True	
							else :
								count=count+1
					if flag==False :
						safe[i]=False
					else :
						safe[i]=True
	except KeyboardInterrupt :
		pass
							
	return comb
	
def occurredBefore(comb,i,key) :
		if testDim(comb[i]) !=0 :
			for	m in xrange(len(comb[i])-1,-1,-1) :
				if comb[i][m] == key :
						return True		
		return False			
# generates the random matrix	
def matrixGenerate(sequences,pop,l):
	matrix=[[random.randrange(0,len(sequences[0])-l+1,1) for i in range(len(sequences))] for j in range(pop)]
	return matrix
# Initialises the particle positions
def initParticle(X) :
	return [[0.0 for j in range(len(X[i]))] for i in range(len(X))]	
# Calculates the Score for the given consensus matrix
def score(matrix) :
	value=0
	transpose=zip(*matrix)
	for i in range(len(transpose)) :
		value=value+max(transpose[i])
	return value	
# Generates the possible combination of children for the given sequences and initial population
def childrenGenerate(sequences,l,localPopulation) :
	children=[]
	dupli=list(sequences)
	original=selectMotifCandidate(dupli,l,localPopulation)
	for k in range(5) :
		dupl=list(sequences)
		for i in range(len(localPopulation)) :
			if len(dupl[i])!=0:
				dupl[i][localPopulation[i]-5:localPopulation[i]+1]=shift(dupl[i][localPopulation[i]-5:localPopulation[i]+1],k+1)
		children.append(selectMotifCandidate(dupl,l,localPopulation))
	for k in range(5) :
		dupl=list(sequences)
		for i in range(len(localPopulation)) :
			dupl[i][localPopulation[i]:localPopulation[i]+6]=shift(dupl[i][localPopulation[i]:localPopulation[i]+6],-(k+1))
		children.append(selectMotifCandidate(dupl,l,localPopulation))		
		
	children.append(original)
	return children
	
# selects the possible list of motif candidates
def selectMotifCandidate(sequences,l,localPopulation) :
	motifCandidate=[]
	for i in range(len(localPopulation)) :
		motifCandidate.append(sequences[i][localPopulation[i]-1:localPopulation[i]+l-1])
	return motifCandidate	


# Selects the best pairs of sequences
def bestPairs(seq,d,l) :
	dupl=list(seq)
	for i in dupl :
		for j in i :
			for k in j :
				spec=shift(k,-1)
				if len(spec)!=l:
					if spec[0] != d :
						j.remove(k)
					
	final=dupl								
	for i in final :
		for j in i :
			for k in j :
				if len(k)>0 and k[len(k)-1]==d:
					f=k.pop()
					
	return final
	
						
def shift(seq, n):
	if len(seq) ==0:
		return seq
	else :		
		n = n % len(seq)
		return seq[n:] + seq[:n]


           
def clever(children,d,l):
	seq=list(children)
	dupl=list(children)
	key=[]
	for i in range(len(dupl)) :
		alignedSites=list(dupl[i])
		flag=True
		for j in range(len(alignedSites)):
			if  len(alignedSites[j]) == 0 :
				flag=False
				break
		if flag==True :			
			matrix=list(alignmentMatrix(alignedSites))
			cons=list(consensus(matrix))
			key=list(dupl[i][len(dupl[i])-1])
			for j in range(len(dupl[i])) :
				del seq[i][j][:]
				for k in range(len(alignedSites)) :
					try :
						ch=hamm(alignedSites[k],cons)-d
						if ch >=0 and ch <=5 :
							generateHammSeq(list(alignedSites[k]),len(alignedSites[k])-1,ch,cons,seq[i][j])
							if testDim(seq[i][k]) > 2 :
								seq[i][k] = list(itertools.chain.from_iterable(seq[i][k]))
							for m in range(len(seq[i][k])) :
								if len(seq[i][k][m]) == l :
									seq[i][k][m].append(hamm(seq[i][k][m],cons))
					except IndexError:
						pass					
			for j in range(len(seq[i]))	:
				try :
					if testDim(seq[i][j]) == 0 :
						ch=hamm(key,cons)-d
						generateHammSeq(key,len(key)-1,ch,cons,seq[i][j])
				except IndexError :
					pass
		else :
			pass								
	return seq			
				
def testDim(testlist, dim=0):
   """tests if testlist is a list and how many dimensions it has
   returns -1 if it is no list at all, 0 if list is empty 
   and otherwise the dimensions of it"""
   if isinstance(testlist, list):
      if testlist == []:
          return dim
      dim = dim + 1
      dim = testDim(testlist[0], dim)
      return dim
   else:
      if dim == 0:
          return -1
      else:
          return dim
# Calculates the consensus of a matrix
def consensus(matrix) :
	consensusString=[]
	transpose=zip(*matrix)
	for i in range(len(transpose)) :
		c=transpose[i].index(max(transpose[i]))
		if c==0 :
			consensusString.append('A')
		if c==1 : 
			consensusString.append('C')
		if c==2 :
			consensusString.append('G')
		if c==3 :
			consensusString.append('T')
	return consensusString
# Calculates the Hamming Distance
def hamm(str1, str2):
	diffs = 0
	for ch1, ch2 in zip(str1, str2):
		if ch1 != ch2:
			diffs += 1
	return diffs
# Gives the Alignment of the matrix
def alignmentMatrix(alignedSites) :
	alignment=[[0 for i in range(4)] for j in range(len(alignedSites[0]))]
	columns=zip(*alignedSites)
	temp =[]
	for i in columns :
		temp.append(list(i))
	count=-1
	for i in temp :
		count=count+1
		for j in i :
			if j == 'A' :
				alignment[count][0]=alignment[count][0]+1
			if j == 'C' :
				alignment[count][1]=alignment[count][1]+1
			if j == 'G' :
				alignment[count][2]=alignment[count][2]+1
			if j == 'T' :
				alignment[count][3]=alignment[count][3]+1
	alignment=zip(*alignment)
	for i in range(len(alignment)) :
		alignment[i]=list(alignment[i])						
	return alignment

# Generates the Hamming Sequences
def generateHammSeq(candidate,i,changesLeft,cons,seq):
	motifCandidate=list(candidate)
	if i >=0 :
		if changesLeft==0 :
			seq.append(motifCandidate)
		else :
			if motifCandidate[i] == 'A' :
				if cons[i]=='A' :
					generateHammSeq(motifCandidate,i-1,changesLeft,cons,seq)
				 	
				if cons[i]=='C' :
					motifCandidate[i]='C'
				 	generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
			 	
				if cons[i]=='G' :
					motifCandidate[i]='G'
				 	generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
				 	
				if cons[i]=='T' :
					motifCandidate[i]='T'
				 	generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
				 	
			if motifCandidate[i] == 'C' :
				if cons[i]=='A' :
					motifCandidate[i]='A'
				 	generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
				 	
				if cons[i]=='C' :
					generateHammSeq(motifCandidate,i-1,changesLeft,cons,seq)
				
				if cons[i]=='G' :
					motifCandidate[i]='G'
				 	generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
			 	
				if cons[i]=='T' :
					motifCandidate[i]='T'
				 	generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
			 	
			if motifCandidate[i] == 'G' :
				if cons[i]=='A' :
					motifCandidate[i]='A'
				 	generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
				 	
				if cons[i]=='C' :
					motifCandidate[i]='C'
					generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
					
				if cons[i]=='G' :
					generateHammSeq(motifCandidate,i-1,changesLeft,cons,seq)
					
				if cons[i]=='T' :
					motifCandidate[i]='T'
				 	generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
				 	
			if motifCandidate[i] == 'T' :
				if cons[i]=='A' :
					motifCandidate[i]='A'
				 	generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
			 	
				if cons[i]=='C' :
					motifCandidate[i]='C'
					generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
				
				if cons[i]=='G' :
					motifCandidate[i]='G'
					generateHammSeq(motifCandidate,i-1,changesLeft-1,cons,seq)
				
				if cons[i]=='T' :
					generateHammSeq(motifCandidate,i-1,changesLeft,cons,seq)
# Reads the input from the File 
def readInput(fileName) :
	with open(fileName) as f:
    		content = f.readlines()
	content = [x.strip() for x in content]
	k=[]
	for i in range(len(content)) :
		k.append(list(content[i]))
	return k					
					
