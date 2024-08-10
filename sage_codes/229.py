import numpy as np
from sage.all import *

import pickle

import sys

sys.setrecursionlimit(10000)

def sg229productg1g2(tuple1, tuple2):
    x1,y1,z1,c1,cp1,t1,m1,p1 = tuple1
    x2,y2,z2,c2,cp2,t2,m2,p2 = tuple2
    return np.array([[x1+x2-c1*x2-cp1*x2-m1*x2+2*c1*m1*x2-2*p1*x2+2*c1*p1*x2+2*cp1*p1*x2+2*m1*p1*x2-4*c1*m1*p1*x2+c1*y2-cp1*y2+m1*y2-2*c1*m1*y2-2*c1*p1*y2+2*cp1*p1*y2-2*m1*p1*y2+4*c1*m1*p1*y2-c1*z2+cp1*z2+2*c1*p1*z2-2*cp1*p1*z2,c1*x2-2*c1*cp1*x2+m1*x2-2*c1*m1*x2-2*cp1*m1*x2+4*c1*cp1*m1*x2-2*c1*p1*x2+4*c1*cp1*p1*x2-2*m1*p1*x2+4*c1*m1*p1*x2+4*cp1*m1*p1*x2-8*c1*cp1*m1*p1*x2+y1+y2-c1*y2-2*cp1*y2+2*c1*cp1*y2-m1*y2+2*c1*m1*y2+2*cp1*m1*y2-4*c1*cp1*m1*y2-2*p1*y2+2*c1*p1*y2+4*cp1*p1*y2-4*c1*cp1*p1*y2+2*m1*p1*y2-4*c1*m1*p1*y2-4*cp1*m1*p1*y2+8*c1*cp1*m1*p1*y2-c1*z2+2*c1*cp1*z2+2*c1*p1*z2-4*c1*cp1*p1*z2,cp1*x2-2*c1*cp1*x2-2*cp1*m1*x2+4*c1*cp1*m1*x2-2*cp1*p1*x2+4*c1*cp1*p1*x2+4*cp1*m1*p1*x2-8*c1*cp1*m1*p1*x2-cp1*y2+2*c1*cp1*y2+2*cp1*m1*y2-4*c1*cp1*m1*y2+2*cp1*p1*y2-4*c1*cp1*p1*y2-4*cp1*m1*p1*y2+8*c1*cp1*m1*p1*y2+z1+z2-2*c1*z2-cp1*z2+2*c1*cp1*z2-2*p1*z2+4*c1*p1*z2+2*cp1*p1*z2-4*c1*cp1*p1*z2,c1+c2+cp2*m1-2*c2*cp2*m1-2*c1*(c2+cp2*m1-2*c2*cp2*m1),cp1+cp2-2*cp1*cp2,0,m1+m2-2*m1*m2,p1+p2-2*p1*p2],[x1+x2-c1*x2-cp1*x2-m1*x2+2*c1*m1*x2-2*p1*x2+2*c1*p1*x2+2*cp1*p1*x2+2*m1*p1*x2-4*c1*m1*p1*x2+c1*y2-cp1*y2+m1*y2-2*c1*m1*y2-2*c1*p1*y2+2*cp1*p1*y2-2*m1*p1*y2+4*c1*m1*p1*y2-c1*z2+cp1*z2+2*c1*p1*z2-2*cp1*p1*z2,c1*x2-2*c1*cp1*x2+m1*x2-2*c1*m1*x2-2*cp1*m1*x2+4*c1*cp1*m1*x2-2*c1*p1*x2+4*c1*cp1*p1*x2-2*m1*p1*x2+4*c1*m1*p1*x2+4*cp1*m1*p1*x2-8*c1*cp1*m1*p1*x2+y1+y2-c1*y2-2*cp1*y2+2*c1*cp1*y2-m1*y2+2*c1*m1*y2+2*cp1*m1*y2-4*c1*cp1*m1*y2-2*p1*y2+2*c1*p1*y2+4*cp1*p1*y2-4*c1*cp1*p1*y2+2*m1*p1*y2-4*c1*m1*p1*y2-4*cp1*m1*p1*y2+8*c1*cp1*m1*p1*y2-c1*z2+2*c1*cp1*z2+2*c1*p1*z2-4*c1*cp1*p1*z2,cp1*x2-2*c1*cp1*x2-2*cp1*m1*x2+4*c1*cp1*m1*x2-2*cp1*p1*x2+4*c1*cp1*p1*x2+4*cp1*m1*p1*x2-8*c1*cp1*m1*p1*x2-cp1*y2+2*c1*cp1*y2+2*cp1*m1*y2-4*c1*cp1*m1*y2+2*cp1*p1*y2-4*c1*cp1*p1*y2-4*cp1*m1*p1*y2+8*c1*cp1*m1*p1*y2+z1+z2-2*c1*z2-cp1*z2+2*c1*cp1*z2-2*p1*z2+4*c1*p1*z2+2*cp1*p1*z2-4*c1*cp1*p1*z2,c1+c2+cp2*m1-2*c2*cp2*m1-2*c1*(c2+cp2*m1-2*c2*cp2*m1),cp1+cp2-2*cp1*cp2,1+m1,m1+m2-2*m1*m2,p1+p2-2*p1*p2],[x1+x2-c1*x2-cp1*x2-m1*x2+2*c1*m1*x2-2*p1*x2+2*c1*p1*x2+2*cp1*p1*x2+2*m1*p1*x2-4*c1*m1*p1*x2+c1*y2-cp1*y2+m1*y2-2*c1*m1*y2-2*c1*p1*y2+2*cp1*p1*y2-2*m1*p1*y2+4*c1*m1*p1*y2-c1*z2+cp1*z2+2*c1*p1*z2-2*cp1*p1*z2,c1*x2-2*c1*cp1*x2+m1*x2-2*c1*m1*x2-2*cp1*m1*x2+4*c1*cp1*m1*x2-2*c1*p1*x2+4*c1*cp1*p1*x2-2*m1*p1*x2+4*c1*m1*p1*x2+4*cp1*m1*p1*x2-8*c1*cp1*m1*p1*x2+y1+y2-c1*y2-2*cp1*y2+2*c1*cp1*y2-m1*y2+2*c1*m1*y2+2*cp1*m1*y2-4*c1*cp1*m1*y2-2*p1*y2+2*c1*p1*y2+4*cp1*p1*y2-4*c1*cp1*p1*y2+2*m1*p1*y2-4*c1*m1*p1*y2-4*cp1*m1*p1*y2+8*c1*cp1*m1*p1*y2-c1*z2+2*c1*cp1*z2+2*c1*p1*z2-4*c1*cp1*p1*z2,cp1*x2-2*c1*cp1*x2-2*cp1*m1*x2+4*c1*cp1*m1*x2-2*cp1*p1*x2+4*c1*cp1*p1*x2+4*cp1*m1*p1*x2-8*c1*cp1*m1*p1*x2-cp1*y2+2*c1*cp1*y2+2*cp1*m1*y2-4*c1*cp1*m1*y2+2*cp1*p1*y2-4*c1*cp1*p1*y2-4*cp1*m1*p1*y2+8*c1*cp1*m1*p1*y2+z1+z2-2*c1*z2-cp1*z2+2*c1*cp1*z2-2*p1*z2+4*c1*p1*z2+2*cp1*p1*z2-4*c1*cp1*p1*z2,c1+c2+cp2*m1-2*c2*cp2*m1-2*c1*(c2+cp2*m1-2*c2*cp2*m1),cp1+cp2-2*cp1*cp2,2-m1,m1+m2-2*m1*m2,p1+p2-2*p1*p2],[x1+c1*x2-cp1*x2-2*c1*m1*x2+2*cp1*m1*x2-2*c1*p1*x2+2*cp1*p1*x2+4*c1*m1*p1*x2-4*cp1*m1*p1*x2-c1*y2+cp1*y2+2*c1*m1*y2-2*cp1*m1*y2+2*c1*p1*y2-2*cp1*p1*y2-4*c1*m1*p1*y2+4*cp1*m1*p1*y2+z2-c1*z2-cp1*z2-2*p1*z2+2*c1*p1*z2+2*cp1*p1*z2,x2-c1*x2-2*cp1*x2+2*c1*cp1*x2-m1*x2+2*cp1*m1*x2-2*p1*x2+2*c1*p1*x2+4*cp1*p1*x2-4*c1*cp1*p1*x2+2*m1*p1*x2-4*cp1*m1*p1*x2+y1-c1*y2+2*c1*cp1*y2+m1*y2-2*cp1*m1*y2+2*c1*p1*y2-4*c1*cp1*p1*y2-2*m1*p1*y2+4*cp1*m1*p1*y2+c1*z2-2*c1*cp1*z2-2*c1*p1*z2+4*c1*cp1*p1*z2,-cp1*x2+2*c1*cp1*x2+m1*x2-2*c1*m1*x2+2*cp1*p1*x2-4*c1*cp1*p1*x2-2*m1*p1*x2+4*c1*m1*p1*x2+y2-2*c1*y2-cp1*y2+2*c1*cp1*y2-m1*y2+2*c1*m1*y2-2*p1*y2+4*c1*p1*y2+2*cp1*p1*y2-4*c1*cp1*p1*y2+2*m1*p1*y2-4*c1*m1*p1*y2+z1+cp1*z2-2*c1*cp1*z2-2*cp1*p1*z2+4*c1*cp1*p1*z2,c1+c2-2*c1*(c2+cp2*(1-m1)-2*c2*cp2*(1-m1))+cp2*(1-m1)-2*c2*cp2*(1-m1),c2+cp1+cp2*m1-2*c2*cp2*m1-2*cp1*(c2+cp2*m1-2*c2*cp2*m1),1,m1+m2-2*m1*m2,p1+p2-2*p1*p2],[x1+c1*x2-cp1*x2-2*c1*m1*x2+2*cp1*m1*x2-2*c1*p1*x2+2*cp1*p1*x2+4*c1*m1*p1*x2-4*cp1*m1*p1*x2-c1*y2+cp1*y2+2*c1*m1*y2-2*cp1*m1*y2+2*c1*p1*y2-2*cp1*p1*y2-4*c1*m1*p1*y2+4*cp1*m1*p1*y2+z2-c1*z2-cp1*z2-2*p1*z2+2*c1*p1*z2+2*cp1*p1*z2,x2-c1*x2-2*cp1*x2+2*c1*cp1*x2-m1*x2+2*cp1*m1*x2-2*p1*x2+2*c1*p1*x2+4*cp1*p1*x2-4*c1*cp1*p1*x2+2*m1*p1*x2-4*cp1*m1*p1*x2+y1-c1*y2+2*c1*cp1*y2+m1*y2-2*cp1*m1*y2+2*c1*p1*y2-4*c1*cp1*p1*y2-2*m1*p1*y2+4*cp1*m1*p1*y2+c1*z2-2*c1*cp1*z2-2*c1*p1*z2+4*c1*cp1*p1*z2,-cp1*x2+2*c1*cp1*x2+m1*x2-2*c1*m1*x2+2*cp1*p1*x2-4*c1*cp1*p1*x2-2*m1*p1*x2+4*c1*m1*p1*x2+y2-2*c1*y2-cp1*y2+2*c1*cp1*y2-m1*y2+2*c1*m1*y2-2*p1*y2+4*c1*p1*y2+2*cp1*p1*y2-4*c1*cp1*p1*y2+2*m1*p1*y2-4*c1*m1*p1*y2+z1+cp1*z2-2*c1*cp1*z2-2*cp1*p1*z2+4*c1*cp1*p1*z2,c1+c2-2*c1*(c2+cp2*(1-m1)-2*c2*cp2*(1-m1))+cp2*(1-m1)-2*c2*cp2*(1-m1),c2+cp1+cp2*m1-2*c2*cp2*m1-2*cp1*(c2+cp2*m1-2*c2*cp2*m1),2-2*m1,m1+m2-2*m1*m2,p1+p2-2*p1*p2],[x1+c1*x2-cp1*x2-2*c1*m1*x2+2*cp1*m1*x2-2*c1*p1*x2+2*cp1*p1*x2+4*c1*m1*p1*x2-4*cp1*m1*p1*x2-c1*y2+cp1*y2+2*c1*m1*y2-2*cp1*m1*y2+2*c1*p1*y2-2*cp1*p1*y2-4*c1*m1*p1*y2+4*cp1*m1*p1*y2+z2-c1*z2-cp1*z2-2*p1*z2+2*c1*p1*z2+2*cp1*p1*z2,x2-c1*x2-2*cp1*x2+2*c1*cp1*x2-m1*x2+2*cp1*m1*x2-2*p1*x2+2*c1*p1*x2+4*cp1*p1*x2-4*c1*cp1*p1*x2+2*m1*p1*x2-4*cp1*m1*p1*x2+y1-c1*y2+2*c1*cp1*y2+m1*y2-2*cp1*m1*y2+2*c1*p1*y2-4*c1*cp1*p1*y2-2*m1*p1*y2+4*cp1*m1*p1*y2+c1*z2-2*c1*cp1*z2-2*c1*p1*z2+4*c1*cp1*p1*z2,-cp1*x2+2*c1*cp1*x2+m1*x2-2*c1*m1*x2+2*cp1*p1*x2-4*c1*cp1*p1*x2-2*m1*p1*x2+4*c1*m1*p1*x2+y2-2*c1*y2-cp1*y2+2*c1*cp1*y2-m1*y2+2*c1*m1*y2-2*p1*y2+4*c1*p1*y2+2*cp1*p1*y2-4*c1*cp1*p1*y2+2*m1*p1*y2-4*c1*m1*p1*y2+z1+cp1*z2-2*c1*cp1*z2-2*cp1*p1*z2+4*c1*cp1*p1*z2,c1+c2-2*c1*(c2+cp2*(1-m1)-2*c2*cp2*(1-m1))+cp2*(1-m1)-2*c2*cp2*(1-m1),c2+cp1+cp2*m1-2*c2*cp2*m1-2*cp1*(c2+cp2*m1-2*c2*cp2*m1),2*m1,m1+m2-2*m1*m2,p1+p2-2*p1*p2],[x1-c1*x2+cp1*x2+m1*x2-2*cp1*m1*x2+2*c1*p1*x2-2*cp1*p1*x2-2*m1*p1*x2+4*cp1*m1*p1*x2+y2-c1*y2-cp1*y2-m1*y2+2*cp1*m1*y2-2*p1*y2+2*c1*p1*y2+2*cp1*p1*y2+2*m1*p1*y2-4*cp1*m1*p1*y2+c1*z2-cp1*z2-2*c1*p1*z2+2*cp1*p1*z2,-c1*x2+2*c1*cp1*x2+2*c1*m1*x2-4*c1*cp1*m1*x2+2*c1*p1*x2-4*c1*cp1*p1*x2-4*c1*m1*p1*x2+8*c1*cp1*m1*p1*x2+y1+c1*y2-2*c1*cp1*y2-2*c1*m1*y2+4*c1*cp1*m1*y2-2*c1*p1*y2+4*c1*cp1*p1*y2+4*c1*m1*p1*y2-8*c1*cp1*m1*p1*y2+z2-c1*z2-2*cp1*z2+2*c1*cp1*z2-2*p1*z2+2*c1*p1*z2+4*cp1*p1*z2-4*c1*cp1*p1*z2,x2-2*c1*x2-cp1*x2+2*c1*cp1*x2-m1*x2+2*c1*m1*x2+2*cp1*m1*x2-4*c1*cp1*m1*x2-2*p1*x2+4*c1*p1*x2+2*cp1*p1*x2-4*c1*cp1*p1*x2+2*m1*p1*x2-4*c1*m1*p1*x2-4*cp1*m1*p1*x2+8*c1*cp1*m1*p1*x2+cp1*y2-2*c1*cp1*y2+m1*y2-2*c1*m1*y2-2*cp1*m1*y2+4*c1*cp1*m1*y2-2*cp1*p1*y2+4*c1*cp1*p1*y2-2*m1*p1*y2+4*c1*m1*p1*y2+4*cp1*m1*p1*y2-8*c1*cp1*m1*p1*y2+z1-cp1*z2+2*c1*cp1*z2+2*cp1*p1*z2-4*c1*cp1*p1*z2,c1+cp2-2*c1*cp2,c2+cp1-2*cp1*(c2+cp2*(1-m1)-2*c2*cp2*(1-m1))+cp2*(1-m1)-2*c2*cp2*(1-m1),2,m1+m2-2*m1*m2,p1+p2-2*p1*p2],[x1-c1*x2+cp1*x2+m1*x2-2*cp1*m1*x2+2*c1*p1*x2-2*cp1*p1*x2-2*m1*p1*x2+4*cp1*m1*p1*x2+y2-c1*y2-cp1*y2-m1*y2+2*cp1*m1*y2-2*p1*y2+2*c1*p1*y2+2*cp1*p1*y2+2*m1*p1*y2-4*cp1*m1*p1*y2+c1*z2-cp1*z2-2*c1*p1*z2+2*cp1*p1*z2,-c1*x2+2*c1*cp1*x2+2*c1*m1*x2-4*c1*cp1*m1*x2+2*c1*p1*x2-4*c1*cp1*p1*x2-4*c1*m1*p1*x2+8*c1*cp1*m1*p1*x2+y1+c1*y2-2*c1*cp1*y2-2*c1*m1*y2+4*c1*cp1*m1*y2-2*c1*p1*y2+4*c1*cp1*p1*y2+4*c1*m1*p1*y2-8*c1*cp1*m1*p1*y2+z2-c1*z2-2*cp1*z2+2*c1*cp1*z2-2*p1*z2+2*c1*p1*z2+4*cp1*p1*z2-4*c1*cp1*p1*z2,x2-2*c1*x2-cp1*x2+2*c1*cp1*x2-m1*x2+2*c1*m1*x2+2*cp1*m1*x2-4*c1*cp1*m1*x2-2*p1*x2+4*c1*p1*x2+2*cp1*p1*x2-4*c1*cp1*p1*x2+2*m1*p1*x2-4*c1*m1*p1*x2-4*cp1*m1*p1*x2+8*c1*cp1*m1*p1*x2+cp1*y2-2*c1*cp1*y2+m1*y2-2*c1*m1*y2-2*cp1*m1*y2+4*c1*cp1*m1*y2-2*cp1*p1*y2+4*c1*cp1*p1*y2-2*m1*p1*y2+4*c1*m1*p1*y2+4*cp1*m1*p1*y2-8*c1*cp1*m1*p1*y2+z1-cp1*z2+2*c1*cp1*z2+2*cp1*p1*z2-4*c1*cp1*p1*z2,c1+cp2-2*c1*cp2,c2+cp1-2*cp1*(c2+cp2*(1-m1)-2*c2*cp2*(1-m1))+cp2*(1-m1)-2*c2*cp2*(1-m1),m1,m1+m2-2*m1*m2,p1+p2-2*p1*p2],[x1-c1*x2+cp1*x2+m1*x2-2*cp1*m1*x2+2*c1*p1*x2-2*cp1*p1*x2-2*m1*p1*x2+4*cp1*m1*p1*x2+y2-c1*y2-cp1*y2-m1*y2+2*cp1*m1*y2-2*p1*y2+2*c1*p1*y2+2*cp1*p1*y2+2*m1*p1*y2-4*cp1*m1*p1*y2+c1*z2-cp1*z2-2*c1*p1*z2+2*cp1*p1*z2,-c1*x2+2*c1*cp1*x2+2*c1*m1*x2-4*c1*cp1*m1*x2+2*c1*p1*x2-4*c1*cp1*p1*x2-4*c1*m1*p1*x2+8*c1*cp1*m1*p1*x2+y1+c1*y2-2*c1*cp1*y2-2*c1*m1*y2+4*c1*cp1*m1*y2-2*c1*p1*y2+4*c1*cp1*p1*y2+4*c1*m1*p1*y2-8*c1*cp1*m1*p1*y2+z2-c1*z2-2*cp1*z2+2*c1*cp1*z2-2*p1*z2+2*c1*p1*z2+4*cp1*p1*z2-4*c1*cp1*p1*z2,x2-2*c1*x2-cp1*x2+2*c1*cp1*x2-m1*x2+2*c1*m1*x2+2*cp1*m1*x2-4*c1*cp1*m1*x2-2*p1*x2+4*c1*p1*x2+2*cp1*p1*x2-4*c1*cp1*p1*x2+2*m1*p1*x2-4*c1*m1*p1*x2-4*cp1*m1*p1*x2+8*c1*cp1*m1*p1*x2+cp1*y2-2*c1*cp1*y2+m1*y2-2*c1*m1*y2-2*cp1*m1*y2+4*c1*cp1*m1*y2-2*cp1*p1*y2+4*c1*cp1*p1*y2-2*m1*p1*y2+4*c1*m1*p1*y2+4*cp1*m1*p1*y2-8*c1*cp1*m1*p1*y2+z1-cp1*z2+2*c1*cp1*z2+2*cp1*p1*z2-4*c1*cp1*p1*z2,c1+cp2-2*c1*cp2,c2+cp1-2*cp1*(c2+cp2*(1-m1)-2*c2*cp2*(1-m1))+cp2*(1-m1)-2*c2*cp2*(1-m1),1-m1,m1+m2-2*m1*m2,p1+p2-2*p1*p2]][3*t1+t2])


def RandXYZ(XYZtuple3):
    x1,y1,z1,x2,y2,z2,x3,y3,z3 = XYZtuple3
    
    #return np.array([1,x1,y1,z1,x2,y2,z2,x3,y3,z3,x1*y1,x1*z1,x1*x2,x1*y2,x1*z2,x1*x3,x1*y3,x1*z3,y1*z1,x2*y1,y1*y2,y1*z2,x3*y1,y1*y3,y1*z3,x2*z1,y2*z1,z1*z2,x3*z1,y3*z1,z1*z3,x2*y2,x2*z2,x2*x3,x2*y3,x2*z3,y2*z2,x3*y2,y2*y3,y2*z3,x3*z2,y3*z2,z2*z3,x3*y3,x3*z3,y3*z3])
    return np.array([1,x1,y1,z1,x2,y2,z2,x3,y3,z3,1/2*x1*(1+x1),1/2*y1*(1+y1),1/2*z1*(1+z1),1/2*x2*(1+x2),1/2*y2*(1+y2),1/2*z2*(1+z2),1/2*x3*(1+x3),1/2*y3*(1+y3),1/2*z3*(1+z3),x1*y1,x1*z1,x1*x2,x1*y2,x1*z2,x1*x3,x1*y3,x1*z3,y1*z1,x2*y1,y1*y2,y1*z2,x3*y1,y1*y3,y1*z3,x2*z1,y2*z1,z1*z2,x3*z1,y3*z1,z1*z3,x2*y2,x2*z2,x2*x3,x2*y3,x2*z3,y2*z2,x3*y2,y2*y3,y2*z3,x3*z2,y3*z2,z2*z3,x3*y3,x3*z3,y3*z3,x1*y1*z1,x1*x2*y1,x1*y1*y2,x1*y1*z2,x1*x3*y1,x1*y1*y3,x1*y1*z3,x1*x2*z1,x1*y2*z1,x1*z1*z2,x1*x3*z1,x1*y3*z1,x1*z1*z3,x1*x2*y2,x1*x2*z2,x1*x2*x3,x1*x2*y3,x1*x2*z3,x1*y2*z2,x1*x3*y2,x1*y2*y3,x1*y2*z3,x1*x3*z2,x1*y3*z2,x1*z2*z3,x1*x3*y3,x1*x3*z3,x1*y3*z3,x2*y1*z1,y1*y2*z1,y1*z1*z2,x3*y1*z1,y1*y3*z1,y1*z1*z3,x2*y1*y2,x2*y1*z2,x2*x3*y1,x2*y1*y3,x2*y1*z3,y1*y2*z2,x3*y1*y2,y1*y2*y3,y1*y2*z3,x3*y1*z2,y1*y3*z2,y1*z2*z3,x3*y1*y3,x3*y1*z3,y1*y3*z3,x2*y2*z1,x2*z1*z2,x2*x3*z1,x2*y3*z1,x2*z1*z3,y2*z1*z2,x3*y2*z1,y2*y3*z1,y2*z1*z3,x3*z1*z2,y3*z1*z2,z1*z2*z3,x3*y3*z1,x3*z1*z3,y3*z1*z3,x2*y2*z2,x2*x3*y2,x2*y2*y3,x2*y2*z3,x2*x3*z2,x2*y3*z2,x2*z2*z3,x2*x3*y3,x2*x3*z3,x2*y3*z3,x3*y2*z2,y2*y3*z2,y2*z2*z3,x3*y2*y3,x3*y2*z3,y2*y3*z3,x3*y3*z2,x3*z2*z3,y3*z2*z3,x3*y3*z3,1/2*x1**2*(1+x1),1/2*x1*(1+x1)*y1,1/2*x1*(1+x1)*z1,1/2*x1*(1+x1)*x2,1/2*x1*(1+x1)*y2,1/2*x1*(1+x1)*z2,1/2*x1*(1+x1)*x3,1/2*x1*(1+x1)*y3,1/2*x1*(1+x1)*z3,1/2*x1*y1*(1+y1),1/2*y1**2*(1+y1),1/2*y1*(1+y1)*z1,1/2*x2*y1*(1+y1),1/2*y1*(1+y1)*y2,1/2*y1*(1+y1)*z2,1/2*x3*y1*(1+y1),1/2*y1*(1+y1)*y3,1/2*y1*(1+y1)*z3,1/2*x1*z1*(1+z1),1/2*y1*z1*(1+z1),1/2*z1**2*(1+z1),1/2*x2*z1*(1+z1),1/2*y2*z1*(1+z1),1/2*z1*(1+z1)*z2,1/2*x3*z1*(1+z1),1/2*y3*z1*(1+z1),1/2*z1*(1+z1)*z3,1/2*x1*x2*(1+x2),1/2*x2*(1+x2)*y1,1/2*x2*(1+x2)*z1,1/2*x2**2*(1+x2),1/2*x2*(1+x2)*y2,1/2*x2*(1+x2)*z2,1/2*x2*(1+x2)*x3,1/2*x2*(1+x2)*y3,1/2*x2*(1+x2)*z3,1/2*x1*y2*(1+y2),1/2*y1*y2*(1+y2),1/2*y2*(1+y2)*z1,1/2*x2*y2*(1+y2),1/2*y2**2*(1+y2),1/2*y2*(1+y2)*z2,1/2*x3*y2*(1+y2),1/2*y2*(1+y2)*y3,1/2*y2*(1+y2)*z3,1/2*x1*z2*(1+z2),1/2*y1*z2*(1+z2),1/2*z1*z2*(1+z2),1/2*x2*z2*(1+z2),1/2*y2*z2*(1+z2),1/2*z2**2*(1+z2),1/2*x3*z2*(1+z2),1/2*y3*z2*(1+z2),1/2*z2*(1+z2)*z3,1/2*x1*x3*(1+x3),1/2*x3*(1+x3)*y1,1/2*x3*(1+x3)*z1,1/2*x2*x3*(1+x3),1/2*x3*(1+x3)*y2,1/2*x3*(1+x3)*z2,1/2*x3**2*(1+x3),1/2*x3*(1+x3)*y3,1/2*x3*(1+x3)*z3,1/2*x1*y3*(1+y3),1/2*y1*y3*(1+y3),1/2*y3*(1+y3)*z1,1/2*x2*y3*(1+y3),1/2*y2*y3*(1+y3),1/2*y3*(1+y3)*z2,1/2*x3*y3*(1+y3),1/2*y3**2*(1+y3),1/2*y3*(1+y3)*z3,1/2*x1*z3*(1+z3),1/2*y1*z3*(1+z3),1/2*z1*z3*(1+z3),1/2*x2*z3*(1+z3),1/2*y2*z3*(1+z3),1/2*z2*z3*(1+z3),1/2*x3*z3*(1+z3),1/2*y3*z3*(1+z3),1/2*z3**2*(1+z3)])
    
    
    
def RandPT(PTtuple3):
    c1,c2,c3,cp1,cp2,cp3,t1,t2,t3,m1,m2,m3,p1,p2,p3 = PTtuple3
    singt = np.array([1,p1,p2,m1,m1*p1,m1*p2,m2,m2*p1,m2*p2,m1*m2,m1*m2*p1,m1*m2*p2,cp1,cp1*p1,cp1*p2,cp1*m1,cp1*m1*p1,cp1*m1*p2,cp1*m2,cp1*m2*p1,cp1*m2*p2,cp1*m1*m2,cp1*m1*m2*p1,cp1*m1*m2*p2,cp2,cp2*p1,cp2*p2,cp2*m1,cp2*m1*p1,cp2*m1*p2,cp2*m2,cp2*m2*p1,cp2*m2*p2,cp2*m1*m2,cp2*m1*m2*p1,cp2*m1*m2*p2,cp1*cp2,cp1*cp2*p1,cp1*cp2*p2,cp1*cp2*m1,cp1*cp2*m1*p1,cp1*cp2*m1*p2,cp1*cp2*m2,cp1*cp2*m2*p1,cp1*cp2*m2*p2,cp1*cp2*m1*m2,cp1*cp2*m1*m2*p1,cp1*cp2*m1*m2*p2,c1,c1*p1,c1*p2,c1*m1,c1*m1*p1,c1*m1*p2,c1*m2,c1*m2*p1,c1*m2*p2,c1*m1*m2,c1*m1*m2*p1,c1*m1*m2*p2,c1*cp1,c1*cp1*p1,c1*cp1*p2,c1*cp1*m1,c1*cp1*m1*p1,c1*cp1*m1*p2,c1*cp1*m2,c1*cp1*m2*p1,c1*cp1*m2*p2,c1*cp1*m1*m2,c1*cp1*m1*m2*p1,c1*cp1*m1*m2*p2,c1*cp2,c1*cp2*p1,c1*cp2*p2,c1*cp2*m1,c1*cp2*m1*p1,c1*cp2*m1*p2,c1*cp2*m2,c1*cp2*m2*p1,c1*cp2*m2*p2,c1*cp2*m1*m2,c1*cp2*m1*m2*p1,c1*cp2*m1*m2*p2,c1*cp1*cp2,c1*cp1*cp2*p1,c1*cp1*cp2*p2,c1*cp1*cp2*m1,c1*cp1*cp2*m1*p1,c1*cp1*cp2*m1*p2,c1*cp1*cp2*m2,c1*cp1*cp2*m2*p1,c1*cp1*cp2*m2*p2,c1*cp1*cp2*m1*m2,c1*cp1*cp2*m1*m2*p1,c1*cp1*cp2*m1*m2*p2,c2,c2*p1,c2*p2,c2*m1,c2*m1*p1,c2*m1*p2,c2*m2,c2*m2*p1,c2*m2*p2,c2*m1*m2,c2*m1*m2*p1,c2*m1*m2*p2,c2*cp1,c2*cp1*p1,c2*cp1*p2,c2*cp1*m1,c2*cp1*m1*p1,c2*cp1*m1*p2,c2*cp1*m2,c2*cp1*m2*p1,c2*cp1*m2*p2,c2*cp1*m1*m2,c2*cp1*m1*m2*p1,c2*cp1*m1*m2*p2,c2*cp2,c2*cp2*p1,c2*cp2*p2,c2*cp2*m1,c2*cp2*m1*p1,c2*cp2*m1*p2,c2*cp2*m2,c2*cp2*m2*p1,c2*cp2*m2*p2,c2*cp2*m1*m2,c2*cp2*m1*m2*p1,c2*cp2*m1*m2*p2,c2*cp1*cp2,c2*cp1*cp2*p1,c2*cp1*cp2*p2,c2*cp1*cp2*m1,c2*cp1*cp2*m1*p1,c2*cp1*cp2*m1*p2,c2*cp1*cp2*m2,c2*cp1*cp2*m2*p1,c2*cp1*cp2*m2*p2,c2*cp1*cp2*m1*m2,c2*cp1*cp2*m1*m2*p1,c2*cp1*cp2*m1*m2*p2,c1*c2,c1*c2*p1,c1*c2*p2,c1*c2*m1,c1*c2*m1*p1,c1*c2*m1*p2,c1*c2*m2,c1*c2*m2*p1,c1*c2*m2*p2,c1*c2*m1*m2,c1*c2*m1*m2*p1,c1*c2*m1*m2*p2,c1*c2*cp1,c1*c2*cp1*p1,c1*c2*cp1*p2,c1*c2*cp1*m1,c1*c2*cp1*m1*p1,c1*c2*cp1*m1*p2,c1*c2*cp1*m2,c1*c2*cp1*m2*p1,c1*c2*cp1*m2*p2,c1*c2*cp1*m1*m2,c1*c2*cp1*m1*m2*p1,c1*c2*cp1*m1*m2*p2,c1*c2*cp2,c1*c2*cp2*p1,c1*c2*cp2*p2,c1*c2*cp2*m1,c1*c2*cp2*m1*p1,c1*c2*cp2*m1*p2,c1*c2*cp2*m2,c1*c2*cp2*m2*p1,c1*c2*cp2*m2*p2,c1*c2*cp2*m1*m2,c1*c2*cp2*m1*m2*p1,c1*c2*cp2*m1*m2*p2,c1*c2*cp1*cp2,c1*c2*cp1*cp2*p1,c1*c2*cp1*cp2*p2,c1*c2*cp1*cp2*m1,c1*c2*cp1*cp2*m1*p1,c1*c2*cp1*cp2*m1*p2,c1*c2*cp1*cp2*m2,c1*c2*cp1*cp2*m2*p1,c1*c2*cp1*cp2*m2*p2,c1*c2*cp1*cp2*m1*m2,c1*c2*cp1*cp2*m1*m2*p1,c1*c2*cp1*cp2*m1*m2*p2])
    return np.concatenate([((t1==0)*1)*singt,((t1==1)*1)*singt,((t1==2)*1)*singt])
    #return np.concatenate([((t1==0 and t2==0)*1)*singt,((t1==0 and t2==1)*1)*singt,((t1==0 and t2==2)*1)*singt,((t1==1 and t2==0)*1)*singt,((t1==1 and t2==1)*1)*singt,((t1==1 and t2==2)*1)*singt,((t1==2 and t2==0)*1)*singt,((t1==2 and t2==1)*1)*singt,((t1==2 and t2==2)*1)*singt])
    
def ToBeCycMat(tuple1, tuple2, tuple3):
    x1,y1,z1,c1,cp1,t1,m1,p1 = tuple1
    x2,y2,z2,c2,cp2,t2,m2,p2 = tuple2
    x3,y3,z3,c3,cp3,t3,m3,p3 = tuple3
    return np.reshape(np.outer(RandXYZ(np.array([x1,y1,z1,x2,y2,z2,x3,y3,z3])),RandPT(np.array([c1,c2,c3,cp1,cp2,cp3,t1,t2,t3,m1,m2,m3,p1,p2,p3]))),(-1,))



def diff3Cocycle(XYZtuple4, PTtuple4, PC3tuple4):
    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4 = XYZtuple4
    c1,cp1,m1,p1,c2,cp2,m2,p2,c3,cp3,m3,p3,c4,cp4,m4,p4 = PTtuple4
    t1,t2,t3,t4 = PC3tuple4
    
    return ToBeCycMat(np.array([x2,y2,z2,c2,cp2,t2,m2,p2]),np.array([x3,y3,z3,c3,cp3,t3,m3,p3]),np.array([x4,y4,z4,c4,cp4,t4,m4,p4]))+\
                      ToBeCycMat(np.array([x1,y1,z1,c1,cp1,t1,m1,p1]),np.array([x2,y2,z2,c2,cp2,t2,m2,p2]),sg229productg1g2(np.array([x3,y3,z3,c3,cp3,t3,m3,p3]),np.array([x4,y4,z4,c4,cp4,t4,m4,p4])))+\
                      ToBeCycMat(np.array([x1,y1,z1,c1,cp1,t1,m1,p1]),sg229productg1g2(np.array([x2,y2,z2,c2,cp2,t2,m2,p2]),np.array([x3,y3,z3,c3,cp3,t3,m3,p3])),np.array([x4,y4,z4,c4,cp4,t4,m4,p4]))+\
                      ToBeCycMat(sg229productg1g2(np.array([x1,y1,z1,c1,cp1,t1,m1,p1]),np.array([x2,y2,z2,c2,cp2,t2,m2,p2])),np.array([x3,y3,z3,c3,cp3,t3,m3,p3]),np.array([x4,y4,z4,c4,cp4,t4,m4,p4]))+\
                      ToBeCycMat(np.array([x1,y1,z1,c1,cp1,t1,m1,p1]),np.array([x2,y2,z2,c2,cp2,t2,m2,p2]),np.array([x3,y3,z3,c3,cp3,t3,m3,p3]))

def randdiffGen():
    return diff3Cocycle(np.random.choice(a=[0,1,2,3], size=(12,)), np.random.choice(a=[0,1], size=(16,)), np.random.choice(a=[0,1,2], size=(4,)))


def Cnonstandard(tuple1,tuple2,tuple3):
    x1,y1,z1,c1,cp1,t1,m1 = tuple1
    x2,y2,z2,c2,cp2,t2,m2 = tuple2
    x3,y3,z3,c3,cp3,t3,m3 = tuple3
    return

def diff3CocycleCnonstandard(XYZtuple4, PTtuple4):
    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4 = XYZtuple4
    c1,cp1,t1,m1,c2,cp2,t2,m2,c3,cp3,t3,m3,c4,cp4,t4,m4 = PTtuple4
    return int(Cnonstandard(np.array([x2,y2,z2,c2,cp2,t2,m2]),np.array([x3,y3,z3,c3,cp3,t3,m3]),np.array([x4,y4,z4,c4,cp4,t4,m4]))+
                      Cnonstandard(np.array([x1,y1,z1,c1,cp1,t1,m1]),np.array([x2,y2,z2,c2,cp2,t2,m2]),sg229productg1g2(np.array([x3,y3,z3,c3,cp3,t3,m3]),np.array([x4,y4,z4,c4,cp4,t4,m4])))+
                      Cnonstandard(np.array([x1,y1,z1,c1,cp1,t1,m1]),sg229productg1g2(np.array([x2,y2,z2,c2,cp2,t2,m2]),np.array([x3,y3,z3,c3,cp3,t3,m3])),np.array([x4,y4,z4,c4,cp4,t4,m4]))+
                      Cnonstandard(sg229productg1g2(np.array([x1,y1,z1,c1,cp1,t1,m1]),np.array([x2,y2,z2,c2,cp2,t2,m2])),np.array([x3,y3,z3,c3,cp3,t3,m3]),np.array([x4,y4,z4,c4,cp4,t4,m4]))+
                      Cnonstandard(np.array([x1,y1,z1,c1,cp1,t1,m1]),np.array([x2,y2,z2,c2,cp2,t2,m2]),np.array([x3,y3,z3,c3,cp3,t3,m3])))%2

def randdiffCnonstandard():
    return diff3CocycleCnonstandard(np.random.choice(a=[0,1,2,3], size=(12,)), np.random.choice(a=[0,1], size=(16,)))

#To test the cocycle is a genuine one:
#print([randdiffCnonstandard() for i in range(1000)] == [0 for i in range(1000)])



MM = [randdiffGen() for _ in range(128000)]
MM = np.vstack(MM)


MMod2 = Matrix(MM, ring=GF(2))
ns=MMod2.right_kernel()



with open('229_basis', 'wb') as f:
    pickle.dump(ns.basis(), f)
    
    
#f = open('/Users/tony/Downloads/229_basis', 'rb')
#kers = pickle.load(f)
#kers.basis() #this will print the vectors; we then copy then to 133_basis which is to be read in Mathematica.











