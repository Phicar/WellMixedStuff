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

