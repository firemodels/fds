#Converts a yaml file into FDS SPEC and REAC inputs.
#yaml file must have thermodynamic data inlcuding NASA7 or NASA9 polynomial and diameter and well-depth.
#run from command line as:
#python cantera2fds.py yamlfile.yaml background_species> fdsfile.fds
# where yamlfile.yaml is the name of the yaml file to process and background spcecies is the species name to set as BACKGROUND

import cantera as ct
import numpy as np
import sys
import math

gas = ct.Solution(sys.argv[1])
bg = sys.argv[2]

#gas = ct.Solution('./mechanisms/GRIMECH/grimech30.yaml')
#bg = 'N2'

k_b = 1.380649E-23
t_r = 298.15
r0 = 8.314472

n_species = len(gas.species())

for i in range(n_species):
   if (bg==gas.species(i).name):
      bg_index = i
      break

for i2 in range(n_species):
   if (i2==0):
      i=bg_index
   elif (i2<=bg_index):
      i=i2-1
   elif (i2>bg_index):
      i=i2
   name = gas.species(i).name
   element_list = list(gas.species(i).composition.keys())
   atoms_list = list(gas.species(i).composition.values())
   sigma = gas.species(i).transport.diameter*1E10
   epsok = gas.species(i).transport.well_depth / k_b
   #sigma = 2.75*1E10
   #epsok = 100.0 / k_b
   n_elem = len(atoms_list)
   formula = ''
   for j in range(n_elem):
      formula += element_list[j]+str(atoms_list[j])

   gas_list = list(gas.species(i).thermo.input_data.values())
   gasinit=gas.species(i).name+':1'
   gas.TPX=293.15,101325,gasinit

   Pr = gas.cp_mass * gas.viscosity / gas.thermal_conductivity
   #Pr = 1
   poly = gas_list[0]
   temp_bands= gas_list[1]
   outstr="&SPEC ID='"+name+"',"
   if (i2==0):
      outstr+="BACKGROUND=T,"
   print(outstr)
   print("      PR_GAS=",int(Pr*1000)/1000,",")
   outstr="      FORMULA='"+formula+"',"
   print(outstr)
   print("      SIGMALJ=",int(sigma*1000)/1000,",")
   print("      EPSILONKLJ=",epsok,",")
   outstr="      POLYNOMIAL='"+poly+"',"
   print(outstr)
   band_l = len(temp_bands)
   outstr = ''
   for j in range(band_l):
      outstr+=str(temp_bands[j])+','
   print("      POLYNOMIAL_TEMP=",outstr)
   poly_len=9
   if (poly=='NASA7'):
      poly_len=7
   for j in range(band_l-1):
     if (j==0 and temp_bands[j]>t_r):
        t_r2 = t_r
        t_r = temp_bands[j]
        if (poly=='NASA7'):
           c_p = gas_list[2][j][0]+gas_list[2][j][1]*t_r+\
           gas_list[2][j][2]*t_r**2+gas_list[2][j][3]*t_r**3+gas_list[2][j][4]*t_r**4
           h_f = gas_list[2][j][0]*t_r+0.5*gas_list[2][j][1]*t_r**2+\
           1./3.*gas_list[2][j][2]*t_r**3+0.25*gas_list[2][j][3]*t_r**4+0.2*gas_list[2][j][4]*t_r**5+gas_list[2][j][5]
           h_f = h_f - c_p * (t_r - t_r2)
           h_f = int(h_f * r0 * 100)/100000.
           t_r = t_r2
        else:
           c_p = gas_list[2][j][0]/t_r**2+gas_list[2][j][1]/t_r+gas_list[2][j][2]+gas_list[2][j][3]*t_r+\
           gas_list[2][j][4]*t_r**2+gas_list[2][j][5]*t_r**3+gas_list[2][j][6]*t_r**4
           h_f = -gas_list[2][j][0]/t_r+gas_list[2][j][1]*math.log(t_r)+gas_list[2][j][2]*t_r+0.5*gas_list[2][j][3]*t_r**2+\
           1./3.*gas_list[2][j][4]*t_r**3+0.25*gas_list[2][j][5]*t_r**4+0.2*gas_list[2][j][6]*t_r**5+gas_list[2][j][7]
           h_f = h_f - c_p * (t_r - t_r2)
           h_f = int(h_f * r0 * 100)/100000.
           t_r = t_r2
     if (temp_bands[j]<t_r):
        if (poly=='NASA7'):
           h_f = gas_list[2][j][0]*t_r+0.5*gas_list[2][j][1]*t_r**2+\
           1./3.*gas_list[2][j][2]*t_r**3+0.25*gas_list[2][j][3]*t_r**4+0.2*gas_list[2][j][4]*t_r**5+gas_list[2][j][5]
           h_f = int(h_f * r0 * 100)/100000.
        else:
           h_f = -gas_list[2][j][0]/t_r+gas_list[2][j][1]*math.log(t_r)+gas_list[2][j][2]*t_r+0.5*gas_list[2][j][3]*t_r**2+\
           1./3.*gas_list[2][j][4]*t_r**3.+0.25*gas_list[2][j][5]*t_r**4+0.2*gas_list[2][j][6]*t_r**5+gas_list[2][j][7]
           h_f = int(h_f * r0 * 100)/100000.
     outstr = ''
     for k in range(poly_len):
        outstr+=str(gas_list[2][j][k])+','
     namestr='      POLYNOMIAL_COEFF(1:'+str(k+1)+','+str(j+1)+")="
     print(namestr,outstr)
   print("      ENTHALPY_OF_FORMATION=",h_f,"/")

numreac = len(gas.reactions())

A=[]
Ea=[]
b=[]
A_lowPr=[]
Ea_lowPr=[]
b_lowPr=[]
rlist=[]
plist=[]
reactType=[]
three=[]
A_Troe=[]
T1_Troe=[]
T2_Troe=[]
T3_Troe=[]
efflist=[]
explist=[]
for i in range(numreac):
    rlist.append(list(gas.reaction(i).reactants.items()))
    explist.append(0)
    for j in range(len(list(gas.reaction(i).reactants.items()))):
        explist[i]=explist[i]+list(gas.reaction(i).reactants.items())[j][1]
    plist.append(list(gas.reaction(i).products.items()))
    if gas.reaction(i).reaction_type=='Arrhenius':
        rate=gas.reaction(i).input_data['rate-constant']
        A.append(rate['A'])
        Ea.append(rate['Ea'])
        b.append(rate['b'])
        A_lowPr.append([])
        Ea_lowPr.append([])
        b_lowPr.append([])
        reactType.append('ARRHENIUS')
        three.append(False)
        efflist.append([])
        A_Troe.append([])
        T1_Troe.append([])
        T2_Troe.append([])
        T3_Troe.append([])
    elif gas.reaction(i).reaction_type=='three-body-Arrhenius':
        rate=gas.reaction(i).input_data['rate-constant']
        A.append(rate['A'])
        Ea.append(rate['Ea'])
        b.append(rate['b'])
        A_lowPr.append([])
        Ea_lowPr.append([])
        b_lowPr.append([])
        reactType.append('THREE-BODY-ARRHENIUS')
        three.append(True)
        efflist.append(list(gas.reaction(i).third_body.efficiencies.items()))
        A_Troe.append([])
        T1_Troe.append([])
        T2_Troe.append([])
        T3_Troe.append([])
    elif gas.reaction(i).reaction_type=='falloff-Lindemann':
        rate=gas.reaction(i).input_data['high-P-rate-constant']
        A.append(rate['A'])
        Ea.append(rate['Ea'])
        b.append(rate['b'])
        rate=gas.reaction(i).input_data['low-P-rate-constant']
        A_lowPr.append(rate['A'])
        Ea_lowPr.append(rate['Ea'])
        b_lowPr.append(rate['b'])
        reactType.append('FALLOFF-LINDEMANN')
        three.append(True)
        efflist.append(list(gas.reaction(i).third_body.efficiencies.items()))
        A_Troe.append([])
        T1_Troe.append([])
        T2_Troe.append([])
        T3_Troe.append([])
    elif gas.reaction(i).reaction_type=='falloff-Troe':
        rate=gas.reaction(i).input_data['high-P-rate-constant']
        A.append(rate['A'])
        Ea.append(rate['Ea'])
        b.append(rate['b'])
        rate=gas.reaction(i).input_data['low-P-rate-constant']
        A_lowPr.append(rate['A'])
        Ea_lowPr.append(rate['Ea'])
        b_lowPr.append(rate['b'])
        reactType.append('FALLOFF-TROE')
        three.append(True)
        efflist.append(list(gas.reaction(i).third_body.efficiencies.items()))
        rate=gas.reaction(i).input_data['Troe']
        A_Troe.append(rate['A'])
        T1_Troe.append(rate['T1'])
        if rate['T1'] in rate:
            T2_Troe.append(rate['T2'])
        else:
            T2_Troe.append(-2E20)
        T3_Troe.append(rate['T3'])

for i in range(len(rlist)):

    for ri in range(len(rlist[i])):
        spName1 = rlist[i][ri][0]
        nu1 = rlist[i][ri][1]
        for pi in  range(len(plist[i])):
            spName2 = plist[i][pi][0]
            nu2 = plist[i][pi][1]
            if (spName1 == spName2):
                print(" Reaction with same element both side:", i)

    print(f"&REAC ID='R{i+1}',")
    print("     REACTYPE='",reactType[i],"',", sep='')
    if (gas.reaction(i).reversible):
        print("     REVERSE=T,")
    print("     RADIATIVE_FRACTION=0,")
    if reactType[i]=='ARRHENIUS':
        print("     A=","{:.5e}".format((A[i]*10000*1000**(explist[i]-1))/1E4),",") #Convert (kmol/m3)^(1-efflist) to mol/cm3^(1-efflist)
    elif reactType[i]=='THREE-BODY-ARRHENIUS':
        print("     A=","{:.5e}".format((A[i]*10000*1000**(explist[i]))/1E4),",")  #Convert (kmol/m3)^(-efflist) to mol/cm3^(-efflist)
    elif reactType[i]=='FALLOFF-LINDEMANN' or reactType[i]=='FALLOFF-TROE':
        print("     A=","{:.5e}".format((A[i]*10000*1000**(explist[i]-1))/1E4),",")  #Convert (kmol/m3)^(1-efflist) to mol/cm3^(1-efflist)

    if (three[i]):
        if (len(efflist[i])>0):
            effs=str(np.array(efflist[i])[:,0])
            effn=str(np.array(efflist[i])[:,1])
            effs=effs.replace('[','').replace(']','').replace("' '","','").replace("\n",",")
            effn=effn.replace('[','').replace(']','').replace(" ",",").replace("'","").replace("\n",",")
            print("     THIRD_EFF_ID=",effs,",")
            print("     THIRD_EFF=",effn,",")
    print("     E=",int(Ea[i]*1E6)/1E9,",")         #kJ/mol to J/mol
    print("     N_T=",b[i],",")
    if reactType[i]=='FALLOFF-LINDEMANN' or reactType[i]=='FALLOFF-TROE':
       print("     A_LOW_PR=","{:.5e}".format((A_lowPr[i]*10000*1000**(explist[i]))/1E4),",")
       print("     E_LOW_PR=",int(Ea_lowPr[i]*1E6)/1E9,",")
       print("     N_T_LOW_PR=",b_lowPr[i],",")
    if reactType[i]=='FALLOFF-TROE':
       print("     A_TROE=",A_Troe[i],",")
       print("     T1_TROE=",T1_Troe[i],",")
       if T2_Troe[i] > -1E20:
           print("     T2_TROE=",T2_Troe[i],",")
       print("     T3_TROE=",T3_Troe[i],",")
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



#gas.species(1).name
#gas.species(1).composition
#{'element':val}
#list(gas.species().composition)
#len elements
#gas.speces().thermo.coeffs
#gas.species().thermo.min_temp  max_temp
#gas.species().transport.diameter * 1E10
#gas.species().transport.well_depth / k_b


#x=list(gas.species().thermo.input_data)
#x[0]=poly
#x[1]=temp range
#x[2][0:N] = poly data for range

