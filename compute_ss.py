#!/usr/bin/env python3

import sys
import argparse

#the function return the number of segregating sites in a sample of n sequences conteined in the "a" variable
def diff_letters_nseq(a):
	nseq=len(a)
	nchar=len(a[0])
	#print (nseq, nchar)
	#print (a[0][:])
	tot=0
	for i in range(nchar):
		s=[]
		for j in range(nseq):
			s.append(a[j][i])
		#print(s)
		if len(set(s))>1 and len(set(s))<=4: 
			tot=tot+1
		elif len(set(s))>4:
			print ("Error in comparing char function: base composition problems")
	return tot

#the function return the number of segregating sites private of each population and pairwise shared (a=list of lists of segsites, n=sample size from each pop)
#order: p1m2, p2m1, p1p2 and m1m2 for each pariwise comparison. ie: 4 pop= 6 pairwise comparison= 4 stat * 6 comp= 24 stats  
def priv_and_shared_s(a,n):
	#print(a)
	nseq=len(a)
	nchar=len(a[0])
	if (nseq % n) == 0:
		npop=int(nseq/n)
		#print(npop)
		ris=[0]*(4*int((npop*(npop-1)/2)))
	else:
		print("Incorrect sample size or number of populations")
		return None
	#print (nseq, nchar)
	#print (a[0][:])
	for i in range(nchar):
		s=[]
		for j in range(nseq):
			s.append(a[j][i])
		#print(s,len(set(s)))
		if len(set(s))>1 and len(set(s))<=4:
			c=0
			#for y in range(0,npop*n-(npop),n):
			for y in range(0,npop*n-(n+1),n):
				for z in range(y+n,npop*n-1,n):
					#print(range(y+n,npop*n-1,n))
					#print("between pop: "+str(y)+" and pop: "+str(z))
					l1=list(s[y:(y+n)])
					l2=list(s[z:(z+n)])
					if len(set(l1))>1 and len(set(l2))==1: ris[c]=ris[c]+1#polymorphic in p1 and monomorphic in p2 (Sx1)
					if len(set(l1))==1 and len(set(l2))>1: ris[c+1]=ris[c+1]+1#polymorphic in p2 and monomorphic in p1 (Sx2)
					if len(set(l1))>1 and len(set(l2))>1: ris[c+2]=ris[c+2]+1#polymorphic in both (Ss)
					if len(set(l1))==1 and len(set(l2))==1 and set(l1)!=set(l2): ris[c+3]=ris[c+3]+1#fixed for diff alleles (Sf)
					c=c+4
		elif len(set(s))>4:
			print ("Error in comparing char function: base composition problems")
	return ris


#this function return the number of differences between 2 sequences
def diff_letters(a,b):
    return sum ( a[i] != b[i] for i in range(len(a)) )

#This function initialize a vector of zeros of length z
def list_zeros(z):
	return([0]*z)

#This function return the frequency table for pi using 0:llw categories for intra pop variation and 0:llb categories for pairwise pop comparisons
def do_freq_table_pi(out,llw,llb,npop,nsim,debug):
	out1=[]			
	for j in range(len(out[0])):
		ftw=list_zeros(llw+1)
		ftb=list_zeros(llb+1)
		if j < npop:
			#do freq tab within pops
			for i in range(nsim):
				if out[i][j]<=llw:
					ftw[out[i][j]]=ftw[out[i][j]]+1
				else:
					ftw[llw]=ftw[llw]+1
			out1.append(ftw)
			if debug: out1.append("|")#debug
		else:
			#do freq tab between pops
			for i in range(nsim):
				if out[i][j]<=llb:
					ftb[out[i][j]]=ftb[out[i][j]]+1
				else:
					ftb[llb]=ftb[llb]+1
			out1.append(ftb)
			if debug: out1.append("|")#debug
	#print(out1)
	result = " ".join(" ".join(map(str,l)) for l in out1)
	return result

#This function return the frequency table for segsites using 0:llb categories for pairwise pop comparisons
def do_freq_table_s(out,llb,npop,nsim,debug):
	out1=[]			
	for j in range(len(out[0])):
		ftb=list_zeros(llb+1)
		#do freq tab between pops
		for i in range(nsim):
			#print (out[i][j])
			if out[i][j]<=llb:
				ftb[out[i][j]]=ftb[out[i][j]]+1
			else:
				ftb[llb]=ftb[llb]+1
		out1.append(ftb)
		if debug: out1.append("|")#debug
	#print(out1)
	result = " ".join(" ".join(map(str,l)) for l in out1)
	return result


#MAIN PROGRAM
#-np NUMBER OF POPULATIONS -nc NUMBER OF CHR FOR EACH POP -w NUMBER OF CATEGORIES OF THE WITHIN POP FREQ TABLE -b NUMBER OF CATEGORIES OF THE WHITHIN POP DIFFERENCES 
parser = argparse.ArgumentParser(description='Compute a summary of the number of difference between chr within and between populations simulated with Hudson\'s ms of similar simulators')
parser.add_argument("-np", "--nrpop", dest="npop", help="Number of populations", type=int, default=1)
parser.add_argument("-nc", "--nrchr", dest="nchr", help="Number of chr for each pop", type=int, default=2)
parser.add_argument("-w", "--within", dest="llw", help="nr categories of whithin pop freq table", type=int, default=20)
parser.add_argument("-b", "--between", dest="llb" ,help="nr categories of between pop freq table", type=int, default=40)
parser.add_argument("-s", "--segsitesPartition", dest="segpar", help="compute private and shared segregating sites instead of pairwise diff, with a minimum of two pops", action='store_true', default=False)
parser.add_argument("-d", "--debug", dest="debug", help="print delimiters between pops and pairwise comparison", action='store_true', default=False)

args = parser.parse_args()


#npop=args.npop#number of populations
#nchr=args.nchr#number of chr per pop
#llw=args.llw#maximum number of difference between 2 chr within each pop, higher values are counted in the last column. (NB there are n+1 categories because of the 0)
#llb=args.llb#maximum number of difference between 4 chr comparison between pop pairs, higher values are counted in the last column (NB: with 2 chr only even number of differences are possible)
if args.npop<1 or (args.npop==1 and args.segpar==True):
	print("Number of populations not correct or incompatible with the -s option")
	sys.exit()
ok=0
c=0
m=[]
out=[]
nsim=0
for line in sys.stdin:
	#print (line)
	if "segsites" in line:
		p=line.split(" ")
		p=int(p[1])
		if ok==0:
			del m[:]
			ok=1
			c=0
	elif ok==1:
		if p>0:
			c=c+1
			if c > 1 and c < (args.npop * args.nchr + 2):
				m.append(line.rstrip('\n'))
				if c == (args.npop * args.nchr + 1):
					ok=0
					nsim=nsim+1
					#print(m)
					d=[]
					#matrix loaded
					if args.segpar:
						#compute private and shared segsites
						#d.append(priv_and_shared_s(m,args.nchr))
						#print(priv_and_shared_s(m,args.nchr))
						out.append(priv_and_shared_s(m,args.nchr))
					else:
						if (args.npop>1):
							#compute segregating sites whithin pop and pairwise diff between pops
							for i in range(0,args.npop*args.nchr-1,args.nchr):
								tmp=[]
								for k in range(i,i+args.nchr,1):
									tmp.append(m[k])
								#print(tmp)
								d.append(diff_letters_nseq(tmp))
							#find initial index for all pairwise population comparison
							for i in range(0,args.npop*args.nchr-(args.npop+1),args.nchr):
								for j in range(i+args.nchr,args.npop*args.nchr-1,args.nchr):
									#compare all pairs of chr
									tot=0
									for ii in range(args.nchr):
										for jj in range(args.nchr):
											tot=tot+diff_letters(m[i+ii],m[j+jj])
									d.append(tot)
							out.append(d)
						else:
							d.append(p)
							out.append(d)
		else:
			if args.segpar:
				z=int(4*(args.npop*(args.npop-1)/2))
				out.append(list_zeros(z))
				ok=0
				nsim=nsim+1
			else:
				z=int(args.npop+(args.npop*(args.npop-1)/2))
				out.append(list_zeros(z))
				ok=0
				nsim=nsim+1
#print(out)
#do freq table
if args.segpar:
	print(do_freq_table_s(out,args.llb,args.npop,nsim,args.debug))
else:
	print(do_freq_table_pi(out,args.llw,args.llb,args.npop,nsim,args.debug))


