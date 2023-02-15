#!/usr/bin/python3
# -*- coding: utf-8 -*-
import argparse
import math

from scipy.optimize import brentq
from scipy.integrate import ode

from numpy import real


def getREff(r1,r2):
	return r1*r2/(r1+r2)
	
def getM(r,rho):
	return 4.0*math.pow(r,3.0)*rho*math.pi / 3.0
	
def getMEff(m1,m2):
	return m1*m2/(m1+m2)

def getRhoMat(Y,rEff,nu):
	return 2.0*Y*math.sqrt(rEff) / ( 3.0*(1.0-nu*nu) )

def getA(beta,rhoMat,mEff):
	return 2.0*beta / ( 3.0*math.pow( rhoMat/mEff , 2.0/5.0 ) )

def getBeta(vStar,v):
	return math.pow(vStar,2.0) / math.pow(v,2.0/10.0)
	
def getEpsOfVStar(vStar):
	orderNum=2
	orderDen=5
	
	aList=[1.0,0.501086]
	
	bList=[1.0,0.501086,1.15345,0.577977,0.532178]
	
	num=0.0
	for order in range(orderNum):
		#print "order num={0}".format(str(order))
		num=num + aList[order]*math.pow( vStar , float(order) )
		
	denom=0.0
	for order in range(orderDen):
		#print "order denom={0}".format(str(order))
		denom=denom + bList[order]*math.pow( vStar  , float(order) )
	
	return num/denom
		

def getVStar(eps):
	def epsZero(vStar):
		return getEpsOfVStar(vStar)-eps
	return brentq(epsZero,0.0,100.0,xtol=1.0e-16)


def getAPhysical(eps, r1, r2, rho, Y, nu, v):
	vStar=getVStar(eps)
	beta=getBeta(vStar,v)
	
	A=getA( beta , getRhoMat(Y,getREff(r1,r2),nu), getMEff( getM(r1,rho) , getM(r2,rho) ) )
	return A

def getTElasticApprox(rMin, rho, Y, nu,v):
	tau=3.281*rMin*math.pow( math.sqrt(2.0)*math.pi*rho*( 1.0-nu*nu ) , 2.0/5.0 )\
		/ ( math.pow( Y , 2.0/5.0 )*math.pow( v , 1.0/5.0) )
	return tau

def getT(mEff,rho,A,v,dt):	
	def derivs(t,y):
		if y[0]<0.0:
			return [ y[1] , 0.0 ]
		return [ y[1] , \
				min( 0.0 , (-rho*math.pow( y[0] , 3.0/2.0 ) - 3.0*( A*rho*math.sqrt(y[0])*y[1] )/2.0 )/mEff ) ]
				
	def jac(t,y):
		if y[0]<0.0:
			return [ [0.0,1.0] , [0.0 , 0.0] ]
		return	[
						[0.0 , 1.0],
						[	( -3.0*rho*( math.sqrt( y[0] ) + A*rho*y[1] / ( 2.0*math.sqrt( y[0] ) ) )/2.0 )/mEff ,
							( -3.0*rho*A*math.sqrt(y[0])/2.0 )/mEff
						]
					]
					
	def inContact(xi,xiDot):
		if xi<0:
			return False
		fN=( -rho*math.pow( xi , 3.0/2.0 ) - 3.0*( A*rho*math.sqrt( xi )*xiDot )/2.0 )/mEff
		if ( fN<0.0 ):
			return True
		return False

	y0=[1.0e-12 , v]
	t0=0.0
	
	#solver=ode(derivs,jac).set_integrator('vode', method='adams', with_jacobian=True)
	solver=ode(derivs,jac).set_integrator('dopri5')
	solver.set_initial_value(y0,t0)


	fOut=open('deformation.dat','w')

	while ( solver.successful() and inContact(solver.y[0],solver.y[1]) ):
		solver.integrate(solver.t+dT)
		outStr="{0}\t{1}\t{2}\n".format(str(solver.t),str(solver.y[0]),str( solver.y[1] ))
		fOut.write(outStr)
		#print solver.t, solver.y
		
	fOut.close()
	return (solver.t,-solver.y[1]/v)

	
def getTPhysical(A, r1, r2, rho, Y, nu, v,dT):
	mEff=getMEff(getM(r1,rho),getM(r2,rho))
	rhoMat=getRhoMat(Y,getREff(r1,r2),nu)
	return getT(mEff,rhoMat,A,v,dT)
	
			


if __name__=="__main__":

	parser = argparse.ArgumentParser(prog='getA.py',
											description="Compute the dissipative constant \"A\" for a "\
															"central collision of viscoelastic spheres according "\
															"to a [1/4] pade approximant (see PRE 84 021302 (2011)). "\
															"Additionally, a rough (theoretical) estimate "\
															"for the corresponding contact duration is returned "\
															"(\"tApprox\") by computing the duration of an elastic contact "\
															"for r=min(r1,r2). \"tCNum\" is the numerically obtained "\
															"contact duration for the desired setup. \"epsNum\" "\
															"is the numerically obtained coefficient of restitution "\
															"for the desired setup. It should, up to numerical errors, "\
															"reproduce the input value specified by -e/--eps. If not, "\
															"numerical integration failed. Please let me now... "\
															"For checking purposes, a file \"deformation.dat\" is "\
															"created. The first colloumn contains the time, "\
															"the second the deformation and the third the "\
															"relative velocity. All quantities are measured in "\
															"SI-units. For the moment don't use particle radii smaller "\
															"than 1.0e-4 meters or eps<0.1 for the numerical stuff. "\
															"The computation of \"A\" should always be reliable or fail totally.",
											formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	materialGroup = parser.add_argument_group('material parameters')
	materialGroup.add_argument('--Young',nargs='?', type=float, help='set Youngs\'s modulus',default=None)
	materialGroup.add_argument('--rho',nargs='?', type=float, help='set mass density',default=None)
	materialGroup.add_argument('--nu',nargs='?', type=float, help='set Poisson\'s ratio ',default=None)
	materialGroup.add_argument('-e','--eps',nargs='?', type=float, 
										help='set desired coefficient of restitution ',
										default=None)
										
	setupGroup = parser.add_argument_group('setup')
	setupGroup.add_argument('--r1',nargs='?', type=float, help='set radius of particle 1',default=None)
	setupGroup.add_argument('--r2',nargs='?', type=float, help='set radius of particle 2',default=None)
	setupGroup.add_argument('-v',nargs='?', type=float, help='set impact velocity',default=None)



	args=parser.parse_args()


	eps=args.eps
	r1=args.r1
	r2=args.r2
	rho=args.rho
	Y=args.Young
	nu=args.nu
	v=args.v
	
	print("\ninput parameters:\n-----------------\n")
	print("eps={0}\nr1={1}\nr2={2}\nrho={3}\nY={4}\nnu={5}\nv={6}\n"\
			.format(str(eps),str(r1),str(r2),str(rho),str(Y),str(nu),str(v)))
	print("-----------------\n")
	
	
	A=getAPhysical(eps,r1,r2,rho,Y,nu,v)
	print("A={0}".format(str(A)))
	
	tMinApprox=getTElasticApprox(min([r1,r2]), rho, Y, nu,v)
	print("tApprox={0}\n".format(str(tMinApprox)))

	dT=tMinApprox/1000.0
	[tCNum,epsNum]=getTPhysical(A,r1,r2,rho,Y,nu,v,dT)
	
	print("tCNum={0}".format(str(tCNum)))
	print("epsNum={0}\n".format(str(epsNum)))
