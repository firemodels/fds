#!/usr/bin/python
#McDermott
#2013-08-20
# last modified: 07-29-2015
#
# This version of the mms uses Stokes form of the Navier-Stokes equations to be more consistent
# with how FDS discretizes the momentum equation.
#
# L. Shunn, F. Ham, P. Moin, Verification of variable-density flow solvers using manufactured
#    solutions, J. Comput. Phys. 231 (2012) 3801--3827.
#
# Problem 3 (modified)

from __future__ import division # make floating point division default, e.g., 1/2=0.5
from mpmath import *
from sympy import *

init_printing(use_unicode=True)

x,y,z,t,r,r0,r1,k,w,uf,vf,f,g,D,p,mu,H = symbols('x y z t r r0 r1 k w uf vf f g D p mu H')

z = ( 1 + sin(pi*k*(x-uf*t))*sin(pi*k*(y-vf*t))*cos(pi*w*t) )/( (1+r0/r1) + \
    (1-r0/r1)*sin(pi*k*(x-uf*t))*sin(pi*k*(y-vf*t))*cos(pi*w*t) ) # mixture fraction
r = 1/( z/r1 + (1-z)/r0 ) # density
u = uf + (r1-r0)/r*(-w/(4*k))*cos(pi*k*(x-uf*t))*sin(pi*k*(y-vf*t))*sin(pi*w*t) # velocity x
v = vf + (r1-r0)/r*(-w/(4*k))*sin(pi*k*(x-uf*t))*cos(pi*k*(y-vf*t))*sin(pi*w*t) # velocity y
#p = r*u*v/2 # pressure (Shunn 3)
K  = (u*u + v*v)/2     # kinetic energy
#H = p/r + K # Bernoulli integral (FDS pseudo-pressure)
H = (u-uf)*(v-vf)/2 # pressure solution has been modified
p = r*(H-K)

print 'Finished with definitions ...'

f = r*u # mass flux x
g = r*v # mass flux y
print 'Finished with mass flux ...'

omz = diff(v,x)-diff(u,y) # vorticity z
print 'Finished with vorticity term ...'

Dif = D*(diff(diff(z,x),x)+diff(diff(z,y),y)) # diffusion term in scalar equation
print 'Finished with diffusion term ...'

Div = diff(u,x)+diff(v,y) # divergence from velocity solution
print 'Finished with divergence ...'

T12 = mu*(diff(u,y)+diff(v,x)) # stress tensor components
print 'Finished with tensor T12 ...'
T11 = 2*mu*(diff(u,x)-Div/3)
print 'Finished with tensor T11 ...'
T22 = 2*mu*(diff(v,y)-Div/3)
print 'Finished with tensor T22 ...'

Q_r = simplify( diff(r,t)+diff(f,x)+diff(g,y) ) # continuity
Q_z = diff(r*z,t)+diff(f*z,x)+diff(g*z,y)-Dif # mixture fraction
#Q_u = diff(r*u,t)+diff(r*u*u,x)+diff(r*u*v,y)+diff(p,x)-(diff(T11,x)+diff(T12,y)) # u momentum
#Q_v = diff(r*v,t)+diff(r*v*u,x)+diff(r*v*v,y)+diff(p,y)-(diff(T12,x)+diff(T22,y)) # v momentum
Q_u = diff(u,t)-v*omz-p*diff(1/r,x)+diff(H,x)-(diff(T11,x)+diff(T12,y))/r # u momentum (Stokes form)
Q_v = diff(v,t)+u*omz-p*diff(1/r,y)+diff(H,y)-(diff(T12,x)+diff(T22,y))/r # v momentum (Stokes form)
print 'Finished with source terms ...'

Div_EOS = simplify( (1/r1 -1/r0)*(Dif + Q_z) ) # divergence from EOS

DD = simplify(Div-Div_EOS) # check EOS divergence formula, DD is zero

# print output to screen
print 'Div = ' + str(Div)
print 'Div_EOS = ' + str(Div_EOS)
print 'Div - Div_EOS = ' + str(DD)
print 'Q_r = ' + str(Q_r)
print 'Q_z = ' + str(Q_z)
print 'Q_u = ' + str(Q_u)
print 'Q_v = ' + str(Q_v)

# also write to file
foo=open('shunn3_stokes_mms.out','w')

foo.write('Div = ' + str(Div) + '\n')
foo.write('\n')
foo.write('Div_EOS = ' + str(Div_EOS) + '\n')
foo.write('\n')
foo.write('DD = ' + str(DD) + '\n')
foo.write('\n')
foo.write('Q_r = ' + str(Q_r) + '\n')
foo.write('\n')
foo.write('Q_z = ' + str(Q_z) + '\n')
foo.write('\n')
foo.write('Q_u = ' + str(Q_u) + '\n')
foo.write('\n')
foo.write('Q_v = ' + str(Q_v) + '\n')

foo.close()


