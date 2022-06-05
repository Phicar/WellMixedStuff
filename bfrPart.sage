from sage.symbolic.integration.integral import indefinite_integral
nn = 0
rr = 0
#using formulas
def rStir(n,k,r):
	return stirling_number2(r,k)*stirling_number2(n-r+1,k)*factorial(k-1)
def rBell(n,r):
	if n<r:
		return 0
	s = 0
	for k in range(1,r+1):
		s+=rStir(n,k,r)#stirling_number2(r,k)*stirling_number2(n-r+1,k)*factorial(k-1)
	return s
def bfBell(n,r):
	s = 0
	for p in SetPartitions(n):
		si = True
		for b in p:
			mi = min(b)
			ma = max(b)
			if ma<r or mi>r:
				si = False
		if si:
			s+=1
			#print(p)
	return s
def binomTransfBell(n):
	s =0
	for r in range(0,n+1):
		s+=binomial(n,r)*rBell(n,r)
	return s
#complementary Bell numbers
def compBell(n):
	s =0
	for k in range(0,n+1):
		s+=((-1)**k)*stirling_number2(n,k)
	return s
#checking prop2.1
def binoTra(n):
	s =0
	for a in range(0,n+2):
		for b in range(0,n+2):
			for c in range(0,n+2):
				d = n+1-a-b-c
				if 0<=d and d<n+1:
					t=factorial(n)/(factorial(a)*factorial(b)*factorial(c)*factorial(d))*(-1)**d
					t*=(2^(a))*bell_number(a)*compBell(b)*compBell(c)*bernoulli(d)
					s+=t
					#print("->",t)
	return s
	
#using SetPartitions of sage
def bfStirl(n,r,k):
        s = 0
        for p in SetPartitions(n):
		if len(p)!=k:
			continue
                si = True
                for b in p:
                        mi = min(b)
                        ma = max(b)
                        if ma<r or mi>r:
                                si = False
                if si:
                        s+=1
                        #print(p)
        return s
#Checking thm2.1
def checkThm21(n,r):
	for rr in range(1,r+1):
		print("--->",rr)
		for nn in range(rr,n+1):
			for k in range(1,nn+1):
				lhs = bfStirl(nn,rr,k)
				rhs = stirling_number2(rr,k)*stirling_number2(nn-rr+1,k)*factorial(k-1)
				if lhs!=rhs:
					print(nn,k,rr,lhs,rhs)

#Checking thm2.2
def checkThm22(s,r):
	ss = 0
	for k in range(1,r+1):
		for l in range(1,k+1):
			ss+=stirling_number2(r,k)*binomial(k-1,l-1)*((-1)^(k-l))*l^(s-r)
	return ss
#checking the determinental prop
def Henk(n,r):
	p=[[rBell(a+b,r) for b in range(0,n)] for a in range(0,n)]
	pp = matrix(p)
	print(pp)
	return pp.det()
#touchard polynomial
def Tou(n,x):
	s =0
	for k in range(0,n+1):
		s+=stirling_number2(n,k)*x^k
	return s
def GF(r,x):
	f(x)=(e^x/(e^x-1))*Tou(r,e^x-1)
	for k in range(0,r):
		f(x)=indefinite_integral(f(x),x)
		
	return f
#bruteforce for permutations
def xyPermutation(r,n,k,no):
	s = 0
	for p in Permutations(n):
		pp = p.to_cycles()
		if len(pp)!=k and no:
			continue
		#print(pp,len(pp))
		si = True
		for x in pp:
			if not si:
				break
			spp = set(x)
			ma = max(spp)
			mi = min(spp)
			if mi>r:
				si = False
				break
		if si:
			s+=1
			
	return s
def abStirling(a,b,k):
        s = 0
        for p in Permutations(a+b):
                pp = p.to_cycles()
                if len(pp)!=k:
                        continue
                #print(pp,len(pp))
                si = True
                for x in pp:
                        if not si:
                                break
                        spp = set(x)
                        ma = max(spp)
                        mi = min(spp)
                        if mi>a or ma<a+1:
                                si = False
                                break
                if si:
                        s+=1

        return s
def rwellStir1(n,r,k):
	s = 0
        for p in Permutations(n):
                pp = p.to_cycles()
                if len(pp)!=k:
                        continue
                #print(pp,len(pp))
                si = True
                for x in pp:
                        if not si:
                                break
                        spp = set(x)
                        ma = max(spp)
                        mi = min(spp)
                        if mi>r or ma<r:
                                si = False
                                break
                if si:
                        s+=1

        return s
@CachedFunction
def recurProp41(a,b,k):
	if k>min(a,b):
		return 0
	if k==0:
		if a==b and a==0:
			return 1
	if k==1:
		return factorial(a+b-1)
	s = (a+b-1)*recurProp41(a,b-1,k)
	for i in range(1,a+1):
		s+=binomial(a,i)*factorial(i)*recurProp41(a-i,b-1,k-1)
	return s
def prop42(n,k,r):
	s= 0
	for i in range(0,r):
		for j in range(0,n-r+1):
			s+=binomial(r-1,i)*binomial(n-r,j)*factorial(i+j)*recurProp41(r-1-i,n-r-j,k-1)
	return s	
def rWellFactorial(n,r):
	s = 0
	for k in range(0,n+1):
		s+=prop42(n,k,r)
	return s
def checkCor41(n,r):
	lhs = r*factorial(n-1)
	rhs = 0
	for k in range(0,r):
		rr = r-k
		nn = n-k
		ss = factorial(nn-2)*(rr*nn-1-rr*(rr-1))
		rhs+=binomial(r-1,k)*factorial(k)*ss
	return lhs #-rhs
#Statistics
def Denspart(n,p):
	posr = [1 for r in range(0,n+1)]
	for x in p:
		mi = min(x)
		ma = max(x)
		#print(p,x,mi,ma)
		for i in range(1,mi):
			posr[i] = 0
		for i in range(ma+1,n+1):
			posr[i]=0
	#print("jode")
	return sum(posr)-1
def Densperm(n,pr):
	global esta
	p = pr.to_cycles()
        posr = [1 for r in range(0,n+1)]
        for xx in p:
		x=set(xx)
                mi = min(x)
                ma = max(x)
                #print(p,x,mi,ma)
                for i in range(1,mi):
                        posr[i] = 0
                for i in range(ma+1,n+1):
                        posr[i]=0
        if sum(posr)==2:
		for j in range(0,len(posr)):
			esta[j]+=posr[j]
		#print(p,posr[1:],esta)
        return sum(posr)-1
#bruteforce Dens partition polynomial
def DensPolPar(n,x):
	tal = [0 for i in range(0,n+1)]
	for p in SetPartitions(n):
		pp = Denspart(n,p)
		#print(p,pp,len(tal))
		tal[pp]+=1
	f(x)=0
	for i in range(0,n+1):
		f+=tal[i]*x^i
	return f
def DensPolPer(n,x):
	global esta
	esta = [0 for i in range(0,n+1)]
        tal = [0 for i in range(0,n+1)]
        for p in Permutations(n):
                pp = Densperm(n,p)
                #print(p,pp,len(tal))
                tal[pp]+=1
        f(x)=0
        for i in range(0,n+1):
                f+=tal[i]*x^i
	print(esta,n)
	esta2 = [0 for i in range(0,n+2)]
	for k in range(1,n+1):
		esta2[k]= esta[k-1]*(k-1)+(n+1-k)*esta[k]
	esta2[n+1]=0
	print(esta2)
        return f

	
def coefZeroPar(n):
	s = 0
	for r in range(1,n):
		for i in range(1,n-r+1):
			s+=stirling_number2(n-r,i)*(stirling_number2(r-1,i)*factorial(i)+stirling_number2(r-1,i+1)*factorial(i+1))
	return bell_number(n)-1-s
