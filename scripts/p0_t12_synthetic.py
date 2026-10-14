"""Synthetic local placement checks; never reads target CT."""
import numpy as np
from p0_t12_core import fit,coordinates,candidates,choose
def volume(height):
    d=np.arange(-8,9)[:,None,None]
    return np.rint(20+200*np.exp(-(d-height[None])**2/(2*.7**2))).astype(np.uint8)
def main():
    h=np.full((128,128),1.);r=fit(volume(h));assert r['plane']['accepted'];assert abs(r['plane']['coefficients'][0]-1)<.01
    assert not r['quadratic']['attempted']
    z=np.zeros((17,128,128),np.uint8);assert not fit(z)['plane']['accepted']
    peaks=[{'depth':-1.,'prominence':30,'height':100},{'depth':1.,'prominence':30,'height':100}];assert choose(peaks,0)[0] is None
    assert choose([{'depth':3.,'prominence':30,'height':100}],0)[0] is None
    u,v=np.meshgrid((np.arange(128)-63.5)/63.5,(np.arange(128)-63.5)/63.5);r=fit(volume(.5+.2*u-.2*v));assert r['plane']['accepted']
    curved=fit(volume(.4*u*u+.8*u*v+.4*v*v))
    assert curved['quadratic']['attempted'] and curved['quadratic']['accepted'],curved
    assert curved['quadratic']['fit_residual_holdout']<=.8*curved['plane']['fit_residual_holdout']
    c=coordinates([100,100,100],[1,0,0],[0,0,1],[0,1,0],[1,0,0],1)
    assert np.allclose(c[1,8],100+9.6/7.91)
    assert np.allclose(np.linalg.norm(c[:,1:]-c[:,:-1],axis=0),9.6/7.91)
    assert np.allclose(c[:,8,63,63]+c[:,8,64,64],2*np.array([100,100+9.6/7.91,100]))
    print('PASS: centered sheet, sloped plane, exercised accepted quadratic, ambiguous competitors, displacement rejection, blank data, signed depth spacing and immutable anchor')
if __name__=='__main__':main()
