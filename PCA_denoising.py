import rpy2.robjects as robjects
from rpy2.robjects import r
import sys
position=[]
snp_id=[]
#To store chromosome & position
f=open('/lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/ascat/SnpPositions.tsv')
fline = f.readline()
while fline :
	fline = fline.strip()
	fso = fline.split('\t')
	snp_id.append(fso[0])
	fso[1] = fso[1].replace('chr','')
	position.append(fso[1]+'@'+fso[2])
	fline = f.readline()
f.close()

#read metadata
#for each sample, control matrix should be generated using samples from differnt donors
#control matrix should be made up of samples of the same gender. 

f=open('pca_metadata.txt')
fline = f.readline()
fline = f.readline()
patient=[]
sample=[]
gender=[]
while fline :
	fline = fline.strip()
	fso = fline.split('\t')
	sample.append(fso[0])
	patient.append(fso[1])
	gender.append(fso[2])
	fline = f.readline()
f.close()


dt_patient = dict(zip(sample,patient))
dt_gender = dict(zip(sample,gender))


#PCA denoising code using RPY2 (input : case logR & control logR matrix)

i=0
while i < len(sample) :
	case_NA=set()
	print (sample[i])
	f=open('../step1/' + sample[i] + '.logR')
	fline = f.readline()
	while fline :
		fline = fline.strip()
		fso = fline.split('\t')
		if fso[1] == 'NA' :
			case_NA.add(fso[0])
		fline = f.readline()
	f.close()
	
	a=[]
	b=[]

	f=open('pca_metadata.txt')
	fline = f.readline()
	fline = f.readline()

	control_sample=set()
	while fline :
		fline = fline.strip()
		fso = fline.split('\t')
		if dt_patient.get(sample[i]) != fso[1] and dt_gender.get(sample[i]) == fso[2] :
			
			k=open('../step1/' + fso[0] + '.logR')
			control_sample.add(fso[0])
			kline = k.readline()
			while kline :
				kline = kline.strip()
				kso = kline.split('\t')
				if kso[1] != 'NA' and kso[0] not in case_NA :
					a.append(fso[0] + '@' + kso[0])
					b.append(kso[1])
				kline = k.readline()
			k.close()
		fline = f.readline()
	f.close()

	dt_control=dict(zip(a,b))
	
	g=open(sample[i]+'.matrix.txt','w')
	g.write('sample')
	for control in control_sample :
		g.write('\t' + control)
	g.write('\n')
	
	j=0
	control_NA=set()
	while j < len(position) :
		num=[]
		for control in control_sample :
			tmp=control + '@' + position[j]
			if dt_control.get(tmp):
				num.append(dt_control.get(tmp))
		
		if len(num)==len(control_sample):
			g.write(position[j])
			jj=0
			while jj <len(num):
				#num[jj] = round(float(num[jj]),4)
				g.write('\t' + str(num[jj]))
				jj=jj+1
			g.write('\n')
		else :
			control_NA.add(position[j])
		j=j+1
	g.close()
	
	k=open('../step1/' + sample[i] + '.logR')
	g=open('v2_' + sample[i] + '.logR','w')
	kline = k.readline()
	while kline :
		kline = kline.strip()
		kso = kline.split('\t')
		if kso[1] != 'NA' and kso[0] not in control_NA :
			kso = kline.split('\t')
			kso[1] = round(float(kso[1]),4)

			g.write(kso[0] + '\t' + str(kso[1])  + '\n')
		kline = k.readline()
	k.close()
	g.close()
 
	# the number of PCs 
	# default : the number of control samples
	
	if len(sys.argv) >= 2:

		input_num = int(sys.argv[1])
	else :
		input_num = len(control_sample)


	r.assign('num',input_num)
	r.assign('case_file','v2_' + sample[i] + '.logR')
	r.assign('control_matrix',sample[i]+'.matrix.txt')
	r('case=read.table(case_file,row.names=1)')
	r('control=read.table(control_matrix,row.names=1,header=TRUE)')
	r('sv=svd(control)')
	r('t2=t(as.matrix(case)) %*% sv$u[,1:num]')
	r('t3=t2 %*% t(sv$u[,1:num])')
	r('t4=t(t3)')
	r('t5=case-t4')
	
	r.assign('denoised_logR','pre_' +sample[i]+'_denoised.tumour.LogR.txt')
	r('write.table(t5,denoised_logR,col.names = F,sep="\t")')
	
	k=open('pre_' +sample[i]+'_denoised.tumour.LogR.txt')
	kline = k.readline()

	denoise_pos=[]
	denoise_log=[]
	while kline :
		kline = kline.strip()
		kso = kline.split('\t')
		kso[0] = kso[0].replace('"','')
		denoise_pos.append(kso[0])
		denoise_log.append(kso[1])
		kline = k.readline()
	k.close()

	denoise = dict(zip(denoise_pos,denoise_log))


	g=open(sample[i]+'_denoised.tumour.LogR.txt','w')
	g.write('\t' + 'Chr' + '\t' + 'Position' + '\t' + sample[i] + '\n')
	j=0
	while j < len(snp_id) :
		position2 = position[j].split('@')
		g.write(snp_id[j] + '\t' +position2[0] + '\t' + position2[1])
		if denoise.get(position[j]) :	
			g.write('\t' + denoise.get(position[j]) + '\n')
		else :
			g.write('\t' + 'NA' + '\n')
		j=j+1
	g.close()
	i=i+1
	break

	
			
	
		
