import math

def cp_w(tmp):
   cp = -0.394796083E+05/tmp**2+0.575573102E+03/tmp+0.931782653E+00+0.722271300E-02*tmp-0.734256000E-05*tmp**2+0.495504000E-08*tmp**3-0.133693000E-11*tmp**4
   cp = cp * 8314.472/18.01528

   cp=2000

   return cp

def h_w(tmp):
   h = 0.394796083E+05/tmp+0.575573102E+03*math.log(tmp)+0.931782653E+00*tmp+0.722271300E-02/2*tmp**2-0.734256000E-05/3*tmp**3+0.495504000E-08/4*tmp**4-0.133693000E-11/5*tmp**5-0.330397431E+05
   h = h * 8314.472/18.01528

   h = 2000 * tmp

   return h


y1 = 0.9
y2 = 0
yr = 0
yw = 0.1

cp1 = 1000
cp2 = 1000
cpr = 1000
cpw = 4184
rho10 = 500
rho20 = 500
rhor0 = 200
rhow0 = 1000
hor1 = 1000*1000
hor2 = 1000*1000
horw = 2500 * 1000
ref_t1 = 315 + 273.15
ref_t2 = 430 + 273.15
ref_tw = 100 + 273.15
ref_r1 = 0.0056
ref_r2 = 0.0075
ref_rw = 0.0016
rgas = 8.314472
refrate = 10/60
m12f = 0.4
m2rf = 0.15
hoc = 14988.3674E3

cpg = 179.4778852

e = math.exp(1)

e1 = e * ref_r1 * rgas * ref_t1**2/refrate
a1 = e * ref_r1 * math.exp(e1/(rgas*ref_t1))

e2 = e * ref_r2 * rgas * ref_t2**2/refrate
a2 = e * ref_r2 * math.exp(e2/(rgas*ref_t2))

ew = e * ref_rw * rgas * ref_tw**2/refrate
aw = e * ref_rw * math.exp(ew/(rgas*ref_tw))

rho_s = y1/rho10 + yw/rhow0
rho_s = 1/rho_s
m0 = rho_s * 1E-6
m1 = y1 * m0
m2 = 0
mr = 0
rho1 = y1 * rho_s
rho2 = 0
rho3 = 0
rhow = yw * rho_s
mw = yw * m0
msum = m1 + m2 + mr + mw
dtdt = 5/60
dt = 1

t = 0
Temp = 20
twrite = 12
mcc =0
dsc = 0
mlr1 = 0
mlr2 = 0
mlrr = 0
mlrw = 0
TempMax = 800

h_w_ref = h_w(ref_tw)

f = open('tga_analysis_exact.csv','w')

f.write('s,C,g/g,g/g,g/g,g/g,g/g,1/s,1/s,1/s,1/s,1/s,W/g,W/g\n')
f.write('Time,Temp,Total Mass,component 1 Mass,water Mass,component 2 Mass,residue Mass,Total MLR, component 1 MLR, water MLR, component 2 MLR, residue MLR,MCC,DSC\n')

outstr=[str('{:.4e}'.format(t)),str('{:.4e}'.format(Temp)),str('{:.4e}'.format(msum/m0)),str('{:.4e}'.format(m1/m0)),str('{:.4e}'.format(mw/m0)),str('{:.4e}'.format(m2/m0)),str('{:.4e}'.format(mr/m0)),str('{:.4e}'.format(0)),str('{:.4e}'.format(0)),str('{:.4e}'.format(0)),str('{:.4e}'.format(0)),str('{:.4e}'.format(0)),str('{:.4e}'.format(0)),str('{:.4e}'.format(0))]
outstr2 = ','.join(outstr)+'\n'
f.write(outstr2)

while True:
   TempK = Temp + 273.15
   m1p = m1
   m2p = m2
   mrp = mr
   mwp = mw

   rr1 = rho1 * a1 * math.exp(-e1/(rgas*TempK))
   rr2 = rho2 * a2 * math.exp(-e2/(rgas*TempK))
   rrw = rhow * aw * math.exp(-ew/(rgas*TempK))

   drho1 = rr1 * dt
   drho2 = rr2 * dt
   drhow = rrw * dt

   rho1n = max(0,rho1 - drho1)
   rho2n = max(0,rho2 - drho2)
   rhown = max(0,rhow - drhow)

   m1 = m1p * rho1n / (rho1+1E-60)
   m2 = m2p * rho2n / (rho2+1E-60) + m12f * (m1p - m1)
   mf1 = (1-m12f) * (m1p - m1)
   mr = mr + m2rf * m2p * (1- rho2n / (rho2+1E-60))
   mf2 = (1-m2rf) * m2p * (1- rho2n / (rho2+1E-60))
   mf = mf1 + mf2
   mw = mw * rhown / (rhow+1E-60)
   msum = m1 + m2 + mr + mw

   y1 = m1/(msum + 1E-60)
   y2 = m2/(msum + 1E-60)
   yr = mr/(msum + 1E-60)
   yw = mw/(msum + 1E-60)

   rho_s = y1/(rho10+1E-60) + y2/(rho20+1E-60) + yr/(rhor0+1E-60) + yw/(rhow0+1E-60)

   if (rho_s > 0):
      rho_s = 1/rho_s

   rho1 = y1 * rho_s
   rho2 = y2 * rho_s
   rhor = yr * rho_s
   rhow = yw * rho_s
   if (rho1<=0): rho1n = 0
   if (rho2<=0): rho2n = 0
   if (rhor<=0): rhorn = 0
   if (rhow<=0): rhown = 0
   mlr1 = (m1p - m1 )/(dt*(m0+1E-60))
   mlr2 = (m2p - m2 )/(dt*(m0+1E-60))
   mlrr = (mrp - mr )/(dt*(m0+1E-60))
   mlrw = (mwp - mw )/(dt*(m0+1E-60))
   mlrtot = mlr1 + mlr2 + mlrr + mlrw

   mcc = hoc * mf / dt * 0.001 / m0

   t = t + dt
   Temp = Temp + dtdt*dt

   dsc = dtdt * dt * (y1 * cp1 + y2 * cp2 + yr * cpr + yw * cpw) * msum /(m0 * 1000)
   dsc = dsc + (m1p - m1) * (hor1 - (TempK - ref_t1) * (cpg - cp1) * (1 - m12f)) / (m0 * 1000)
   dsc = dsc + (mr - mrp)/m2rf * (hor2 - (TempK - ref_t2) * (cpg - cp2) * (1 - m2rf)) / (m0 * 1000)
   h_w_TempK = h_w(TempK)
   dsc = dsc + (mwp - mw) * (horw + h_w_TempK - h_w_ref - (TempK - ref_tw) * cpw) / (m0 * 1000)

   if (t >= twrite):
      outstr=[str('{:.4e}'.format(t)),str('{:.4e}'.format(Temp)),str('{:.4e}'.format(msum/m0)),str('{:.4e}'.format(m1/m0)),str('{:.4e}'.format(mw/m0)),str('{:.4e}'.format(m2/m0)),str('{:.4e}'.format(mr/m0)),str('{:.4e}'.format(mlrtot)),str('{:.4e}'.format(mlr1)),str('{:.4e}'.format(mlr2)),str('{:.4e}'.format(mlrw)),str('{:.4e}'.format(mlrr)),str('{:.4e}'.format(mcc)),str('{:.4e}'.format(dsc))]
      outstr2 = ','.join(outstr)+'\n'
      f.write(outstr2)
      twrite = twrite+12
   if (TempK >= TempMax):
      break



