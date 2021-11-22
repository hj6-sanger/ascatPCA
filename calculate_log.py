import math
f=open('ascat_sample.txt')
fline = f.readline()
fline = f.readline()
while fline :
	fline = fline.strip()
	fso = fline.split('\t')
	#normal
	a=[]
	#tumor
	b=[]
	#chr_pos
	pos=[]

	
	#normal_count_file
	k=open(fso[0] + '.count')
	kline = k.readline()
	kline = k.readline()
	tot1=0
	tot2=0
	while kline :
		kline = kline.strip()
		kso = kline.split('\t')
		
		a.append(float(kso[6]))
		pos.append(kso[0]+'@'+kso[1])
		kline = k.readline()
	k.close()
	tot1=sum(a)
	
	#tumor_count_file
	k=open(fso[1] + '.count')
	kline = k.readline()
	kline = k.readline()
	while kline :
		kline = kline.strip()
		kso = kline.split('\t')
		b.append(float(kso[6]))	
		kline = k.readline()
	k.close()
	tot2=sum(b)
	g=open(fso[1] +'.logR','w')
	i=0
	#calculation of logR for each position
	while i < len(a) :
		if a[i] == 0 or b[i] == 0 :

			g.write(pos[i] + '\t' + 'NA' + '\n')
		else :
			a1=a[i]/tot1
			b1=b[i]/tot2
			ans = math.log2(b1/a1)
			ans = round(ans,4)
			g.write(pos[i] + '\t' + str(ans) + '\n')
		i=i+1
	g.close()
	fline = f.readline()
f.close()	
