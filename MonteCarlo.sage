def tal(n,m):
	s = 0
	p = Permutations(n)
	for k in range(0,m):
		x = p.random_element()
		xx = x.to_cycles()
		perm = [1 for i in range(1,n+1)]
		for c in xx:
			mc = set(c)
			mi = min(mc)
			ma = max(mc)
			for r in range(1,mi):
				perm[r-1]=0
			for r in range(ma+1,n+1):
				perm[r-1]=0
		st = sum(perm)
		s+=st
	return s/m	
