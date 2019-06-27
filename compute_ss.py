#!/usr/bin/env python3

import sys
import argparse

#the function return the number of segregating sites in a sample of n sequences
def diff_letters_nseq(a):
	nseq=len(a)
	nchar=len(a[0])
	tot=0
	for i in range(nchar):
		s=[]
		for j in range(nseq):
			s.append(a[j][i])
		if len(set(s))>1 and len(set(s))<=4: 
			tot=tot+1
		elif len(set(s))>4:
			print ("Error in comparing char function: base composition problems")
	return tot

#the function return the number of segregating sites private of each population and pairwise shared (a=list of lists of segsites, n=vector of sample sizes from each pop, tn=total number of chrs)
#order: p1m2, p2m1, p1p2 and m1m2 for each pariwise comparison. ie: 4 pop= 6 pairwise comparison= 4 stat * 6 comp= 24 stats  
def priv_and_shared_s(a,n,tn):
	#print(a)
	nseq=len(a)
	nchar=len(a[0])
	if (nseq % tn) == 0:
		npop=len(n)
		ris=[0]*(4*int((npop*(npop-1)/2)))
	else:
		print("Incorrect sample size or number of populations")
		return None
	pi=get_starting_points(n)
	for i in range(nchar):
		s=[]
		for j in range(nseq):
			s.append(a[j][i])
		if len(set(s))>1 and len(set(s))<=4:
			c=0
			for y in range(len(n)-1):
				for z in range(y+1,len(n),1):
					l1=list(s[pi[y]:(pi[y]+n[y])])
					l2=list(s[pi[z]:(pi[z]+n[z])])
					if len(set(l1))>1 and len(set(l2))==1: ris[c]=ris[c]+1#polymorphic in p1 and monomorphic in p2 (Sx1)
					if len(set(l1))==1 and len(set(l2))>1: ris[c+1]=ris[c+1]+1#polymorphic in p2 and monomorphic in p1 (Sx2)
					if len(set(l1))>1 and len(set(l2))>1: ris[c+2]=ris[c+2]+1#polymorphic in both (Ss)
					if len(set(l1))==1 and len(set(l2))==1 and set(l1)!=set(l2): ris[c+3]=ris[c+3]+1#fixed for diff alleles (Sf)
					c=c+4
		elif len(set(s))>4:
			print ("Error in comparing char function: base composition problems")
	return ris

#This function extract the starting points given a chromosome set vector
def get_starting_points(n):
	pi=[]
	for i in range(len(n)):
		if i==0:
			pi.append(0)
		else:
			pi.append(pi[i-1]+n[i-1])
	return(pi)

#This function compute the unfolded sfs in single or in all possible pairs of populations using the chr numbers defined in the n vector list. The sfs is unfolded by default (0: ancestral, 1: derived allele) but it could be folded using the fold option.
#(a=list of lists of segsites;n=vector count of chr per pop; tn=total number of chrs)
def sfs(a,n,tn,fold=False):
	nseq=len(a)
	nchar=len(a[0])
	if (nseq % tn) == 0:
		ris=init_sfs(n)
	else:
		print("Incorrect sample size or number of populations")
		return None
	pi=get_starting_points(n)
	for i in range(nchar):
		s=[]
		for j in range(nseq):
			s.append(a[j][i])
		if len(set(s))>1 and len(set(s))<=4:
			c=0
			#check the size of sfs 1D or 2D
			if len(n)==1:
				l1=list(s[pi[0]:(pi[0]+n[0])])
				c1=l1.count('1')
				if fold:
					d1=l1.count('0')
					if d1<c1: c1=d1
				ris[c1]=ris[c1]+1
			else:
				for y in range(len(n)-1):
					for z in range(y+1,len(n),1):
						l1=list(s[pi[y]:(pi[y]+n[y])])
						l2=list(s[pi[z]:(pi[z]+n[z])])
						c1=l1.count('1')
						c2=l2.count('1')
						if fold:
							d1=l1.count('0')
							d2=l2.count('0')
							if d1<c1: c1=d1
							if d2<c2: c2=d2
						e1=int(n[y])
						e2=int(n[z])
						ris[c][((e1+1)*(c2+1)-1)-(e1-c1)]=(ris[c][((e1+1)*(c2+1)-1)-(e1-c1)])+1
						c=c+1
		elif len(set(s))>4:
			print ("Error in comparing char function: base composition problems")
	return ris	

#this function initialize a single unfolded or folded sfs, or the ones that are possible between each pair of pop (n is the vector of chromosome sample size)
def init_sfs(n):
	out=[]
	if len(n)>1:
		for y in range(len(n)-1):
			for z in range(y+1,len(n),1):
				out.append([0]*((int(n[y])+1)*(int(n[z])+1)))
	else:
		out=list_zeros(n[0]+1)
	return out

#this function print a planar sfs given the output of the sfs function. It collapses the spectra produced in multiple simulations to a single combination of spectra.
def collapse_sfs(x,n,debug):
	#print(x)
	nsim=len(x)
	#print(nsim,nsp,n)
	out=[]
	#single simulation
	if nsim==1:
		if n==1: 
			#if debug: x[0].append("|")
			return " ".join(" ".join(map(str,l)) for l in x)
		else:
			#if debug: x[0].append("|")
			return " ".join(" ".join(map(str,l)) for l in x[0])
	else:
		#multiple simulations
		if n==1:#single pop case->1 1D spectrum
			for j in range(1,len(x),1):
				#print(j) 
				if j==1: out.append(([sum(z) for z in zip(x[0], x[1])]))
				else:  out[0]=[sum(z) for z in zip(out[0], x[j])]
		elif n==2:#two pop case -> 1 1D spectrum
			for j in range(1,len(x),1):
				if j==1: out.append(([sum(z) for z in zip(x[0][0], x[1][0])]))
				else:  out[0]=[sum(z) for z in zip(out[0], x[j][0])]
		else:#multi pop case-> multiple 2D spectra
			for i in range(len(x[0])):#iterate over spectra
				for j in range(1,len(x),1):#iterate over simulations
					if j==1: out.append([sum(z) for z in zip(x[0][i], x[j][i])])
					else: out[i]=([sum(z) for z in zip(out[i], x[j][i])])
					#print("aaa",out)
					#if debug: out.append("|")
		return " ".join(" ".join(map(str,l)) for l in out)


#this function return the number of differences between 2 sequences
def diff_letters(a,b):
    return sum ( a[i] != b[i] for i in range(len(a)) )

#This function initialize a vector of zeros of length z
def list_zeros(z):
	return([0]*z)

#This function checks that the vector of genotypes is composed by genotypes ("0 or 1" char)
def check_geno(x):
	for i in range(len(x)):
		if (x[i].count('0')==0) and (x[i].count('1')==0):
			print("Problem reading genotypes in individual:",i+1,"- Check the chromosome vector or snp symbols")
			sys.exit()

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
parser = argparse.ArgumentParser(description='Compute a summary of the number of difference between chr within and between populations simulated with Hudson\'s ms or similar coalescent simulators')
parser.add_argument("-np", "--nrpop", dest="npop", help="Number of populations", type=int, default=1)
parser.add_argument("-nc", "--nrchr", dest="nchr", help="Vector containing the number of chr for each pop (es: 2,2,2 for a 3 population comparison)", type=str, default=2)
parser.add_argument("-w", "--within", dest="llw", help="nr categories of whithin pop freq table", type=int, default=20)
parser.add_argument("-b", "--between", dest="llb" ,help="nr categories of between pop freq table", type=int, default=40)
parser.add_argument("-s", "--segSitesPartition", dest="segpar", help="compute private and shared segregating sites instead of pairwise diff, with a minimum of two pops", action='store_true', default=False)
parser.add_argument("-sfs", "--siteFrequencySpectrum", dest="sfs", help="compute the 1D or 2D unfolded sfs", action='store_true', default=False)
parser.add_argument("-folded", "--folded", dest="fold", help="fold the site frequency spectrum", action='store_true', default=False)
parser.add_argument("-d", "--debug", dest="debug", help="print delimiters between pops and pairwise comparison", action='store_true', default=False)

args = parser.parse_args()


#npop=args.npop#number of populations
#nchr=args.nchr#number of chr per pop
#llw=args.llw#maximum number of difference between 2 chr within each pop, higher values are counted in the last column. (NB there are n+1 categories because of the 0)
#llb=args.llb#maximum number of difference between 4 chr comparison between pop pairs, higher values are counted in the last column (NB: with 2 chr only even number of differences are possible)
if args.npop<1 or (args.npop==1 and args.segpar==True):
	print("Number of populations not correct or incompatible with the -s option")
	sys.exit()
#check chromosome vector
vc=args.nchr.split(",")
if len(vc)<1 or len(vc)!=args.npop:
	print("Wrong number of chromosomes. Formatting error")
	sys.exit()
else:
	vc = [int(i) for i in vc]
	totc=sum(vc)
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
			if c > 1 and c < (totc + 2):
				m.append(line.rstrip('\n'))
				if c == (totc + 1):
					ok=0
					nsim=nsim+1
					#print(m)
					check_geno(m)
					d=[]
					#matrix loaded
					if args.segpar:
						#compute private and shared segsites
						out.append(priv_and_shared_s(m,vc,totc))
					elif args.sfs:
						out.append(sfs(m,vc,totc,fold=args.fold))
					else:
						if (args.npop>1):
							#compute segregating sites whithin pop and pairwise diff between pops
							pi=get_starting_points(vc)
							for i in range(len(pi)):
								tmp=[]
								for k in range(pi[i],pi[i]+vc[i],1):
									tmp.append(m[k])
								d.append(diff_letters_nseq(tmp))
							#find initial index for all pairwise population comparison
							for i in range(len(vc)-1):
								for j in range(i+1,len(vc),1):
									#compare all pairs of chr
									tot=0
									for ii in range(vc[i]):
										for jj in range(vc[j]):
											tot=tot+diff_letters(m[pi[i]+ii],m[pi[j]+jj])
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
			elif args.sfs:
				z=init_sfs(vc)
				out.append(z)
				ok=0
				nsim=nsim+1
			else:
				z=int(args.npop+(args.npop*(args.npop-1)/2))
				out.append(list_zeros(z))
				ok=0
				nsim=nsim+1
#print(out)
#print(args)
#do freq table
if args.segpar:
	#print("fdss")
	print(do_freq_table_s(out,args.llb,args.npop,nsim,args.debug))
elif args.sfs:
	#print("sfs")
	print(collapse_sfs(out,args.npop,args.debug))
else:
	#print("pdiff")
	print(do_freq_table_pi(out,args.llw,args.llb,args.npop,nsim,args.debug))



