# This script is under development. It uses cantera to convert a cantera chemistry file to a set of FDS inputs.
# CAUTION: At this point scaling factors for A and E need to be verified and ensure handling of A for third body and catalytic reactions is correct.
# The cantera file for gas can be either cti or yaml format. Change the file name to the file to be converted.
# The output is a stream of txt to the terminal. Use > to redirect to a file.
# The spec_list needs to have all the species present in the cti or yaml file. Major species are there, but not some species may be identified by different names.
# For example the ISOPROPYL RADICAL could be C3H7, iC3H7, i-C3H7, or something else.
# Long term plans are to allow FDS to read in polynomial coeffiicents on a SPEC line.  Then the spec_list here would not be needed.
import cantera as ct
import numpy as np

gas=ct.Solution('smooke.yaml')

numreac = len(gas.reactions())

spec_list=[]

spec_list.append(['AR','ARGON'])
spec_list.append(['aC3H4','ALLENE'])
spec_list.append(['aC3H5','ALLYL RADICAL'])
spec_list.append(['C','SOOT'])
spec_list.append(['C2H','ETHYNYL RADICAL'])
spec_list.append(['C2H2','ACETYLENE'])
spec_list.append(['C2H3','VINYL RADICAL'])
spec_list.append(['C2H3CHO','ACROLEIN'])
spec_list.append(['C2H4','ETHYLENE'])
spec_list.append(['C2H5','ETHYL RADICAL'])
spec_list.append(['C2H6','ETHANE'])
spec_list.append(['C3H2','CYCLOPROPENYLIDENE'])
spec_list.append(['C3H3','1-PROPENYL'])
spec_list.append(['C3H6','CYCLOPROPANE'])
spec_list.append(['C3H8','PROPANE'])
spec_list.append(['C4H2','BUTADIYNE'])
spec_list.append(['C4H81','1-BUTENE'])
spec_list.append(['CH','METHYLIDYNE'])
spec_list.append(['CH2','METHYLENE'])
spec_list.append(['CH2*','METHYLENEs'])
spec_list.append(['CH2CHO','CH2CHO'])
spec_list.append(['CH2CO','KETENE'])
spec_list.append(['CH2O','FORMALDEHYDE'])
spec_list.append(['CH2OH','HYDROXYMETHYL RADICAL'])
spec_list.append(['CH3','METHYL RADICAL'])
spec_list.append(['CH3CCH2','CYCLOPROPYL'])
spec_list.append(['CH3CHO','ETHANAL'])
spec_list.append(['CH3O','METHOXY RADICAL'])
spec_list.append(['CH3OH','METHANOL'])
spec_list.append(['CH4','METHANE'])
spec_list.append(['CO','CARBON MONOXIDE'])
spec_list.append(['CO2','CARBON DIOXIDE'])
spec_list.append(['H','HYDROGEN ATOM'])
spec_list.append(['H2','HYDROGEN'])
spec_list.append(['H2O','WATER VAPOR'])
spec_list.append(['H2O2','HYDROGEN PEROXIDE'])
spec_list.append(['HCCO','ETHYNYLOXY RADICAL'])
spec_list.append(['HCO','FORMYL RADICAL'])
spec_list.append(['HNO','NITROSYL HYDRIDE'])
spec_list.append(['HO2','HYDROPEROXY RADICAL'])
spec_list.append(['iC3H7','ISOPROPYL RADICAL'])
spec_list.append(['N','NITROGEN ATOM'])
spec_list.append(['N2','NITROGEN'])
spec_list.append(['N2O','NITROUS OXIDE'])
spec_list.append(['NA','SODIUM'])
spec_list.append(['NAH','SODIUM HYDRIDE'])
spec_list.append(['NAO','SODIUM OXIDE'])
spec_list.append(['NAO2','SODIUM SUPEROXIDE'])
spec_list.append(['NAOH','SODIUM HYDROXIDE'])
spec_list.append(['nC3H7','PROPYL RADICAL'])
spec_list.append(['NH','IMIDOGEN'])
spec_list.append(['NNH','NNH'])
spec_list.append(['NO','NITRIC OXIDE'])
spec_list.append(['O','OXYGEN ATOM'])
spec_list.append(['O2','OXYGEN'])
spec_list.append(['OH','HYDROXYL RADICAL'])
spec_list.append(['pC3H4','PROPYNE'])

dumpspec=[False]*len(spec_list)

#for i in range(len(spec_list)):
#	outstr="&SPEC ID='"+spec_list[i][1]+",FYI='"+spec_list[i][0]+"'/"
#	print(outstr)

A=[]
Ea=[]
b=[]
rlist=[]
plist=[]
three=[]
efflist=[]
for i in range(numreac):
	rlist.append(list(gas.reaction(i).reactants.items()))
	plist.append(list(gas.reaction(i).products.items()))
	try:
		rate=gas.reaction(i).rate.input_data['low-P-rate-constant']
		A.append(rate['A'])
		Ea.append(rate['Ea'])
		b.append(rate['b'])

	except (KeyError):
		rate=gas.reaction(i).rate.input_data['rate-constant']
		A.append(rate['A'])
		Ea.append(rate['Ea'])
		b.append(rate['b'])
	if gas.reaction(i).reaction_type=='reaction':
		three.append(False)
		efflist.append([])
	else:
		three.append(True)
		efflist.append(list(gas.reaction(i).efficiencies.items()))
		for j in range(len(efflist[i])):
			for k in range(len(spec_list)):
				if(efflist[i][j][0]==spec_list[k][0]):
					x=list(efflist[i][j])
					dumpspec[k]=True
					x[0]=spec_list[k][1]
					efflist[i][j]=tuple(x)
					break

for i in range(len(spec_list)):
	for j in range(len(rlist)):
		for k in range(len(rlist[j])):
			rlist[j][k]=list(rlist[j][k])
			if rlist[j][k][0]==spec_list[i][0]:
				dumpspec[i]=True
				rlist[j][k][0]=spec_list[i][1]
	for j in range(len(plist)):
		for k in range(len(plist[j])):
			plist[j][k]=list(plist[j][k])
			if plist[j][k][0]==spec_list[i][0]:
				dumpspec[i]=True
				plist[j][k][0]=spec_list[i][1]

for i in range(len(spec_list)):
	if dumpspec[i]:
		outstr="&SPEC ID='"+spec_list[i][1]+"',FYI='"+spec_list[i][0]+"'/"
		print(outstr)

for i in range(len(rlist)):
	print(f"&REAC ID='R{i+1}',")
	print("     REVERSE=T,")
	print("     RADIATIVE_FRACTION=0,")
	if (three[i]):
		print("     THIRD_BODY=T,")
		print("     A=",A[i]*1E6,",")
		if (len(efflist[i])>0):
			effs=str(np.array(efflist[i])[:,0])
			effn=str(np.array(efflist[i])[:,1])
			effs=effs.replace('[','').replace(']','').replace("' '","','").replace("\n",",")
			effn=effn.replace('[','').replace(']','').replace(" ",",").replace("'","").replace("\n",",")
			print("     THIRD_EFF_ID=",effs,",")
			print("     THIRD_EFF=",effn,",")
	else:
		print("     A=",A[i]*1E3,",")
	print("     E=",Ea[i]*0.001,",")
	if (b[i]!=0): print("     N_T=",b[i],",")
	rlist2=str(np.array(rlist[i])[:,0])
	plist2=str(np.array(plist[i])[:,0])
	rlist2=rlist2.replace('[','').replace(']','').replace("' '","','")
	plist2=plist2.replace('[','').replace(']','').replace("' '","','")
	print("     SPEC_ID_NU=",rlist2,",",plist2,",")
	rlist2=str(np.array(rlist[i])[:,1])
	plist2=str(np.array(plist[i])[:,1])
	rlist2=rlist2.replace('[','-').replace(']','').replace(" ",",-").replace("'","")
	plist2=plist2.replace('[','').replace(']','').replace(" ",",").replace("'","")
	print("     NU=",rlist2,",",plist2,"/")

