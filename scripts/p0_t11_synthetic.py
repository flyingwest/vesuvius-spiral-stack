"""Synthetic orientation/axis checks; no target CT read."""
import numpy as np
import p0_t11_core as c
def main():
    z,y,x=np.meshgrid(*[np.arange(-16,17)]*3,indexing='ij')
    for tilt in [-.6,0,.6]:
        n=np.array([0.,np.sqrt(1-tilt*tilt),tilt]);v=(200*np.exp(-(n[1]*y+n[2]*z)**2/8)).astype(np.float32)
        c.cube=lambda P:v
        for sigma in [1.,2.]:
            f=c.frame([0,0,0],[1,0,0],[0,1,0],sigma)
            assert abs(f['n']@n)>.999
            assert np.allclose(np.cross(f['t_s'],f['t_z']),f['n'])
            assert f['image_y'][2]>0 and f['n'][1]>0
    c.cube=lambda P:np.zeros((33,33,33),np.float32)
    assert not c.frame([0,0,0],[1,0,0],[0,1,0],1.)['accepted']
    # Independent affine CT field checks XYZ ordering and centered offsets.
    c.chunk=lambda k,j,i:(np.indices((128,128,128))[2]+128*i).astype(np.uint8)
    a=c.render([64,64,64],[1,0,0],[0,0,1],[0,1,0])
    expected=np.rint(64+(np.arange(128)-63.5)*c.STEP).clip(0,255).astype(np.uint8)
    # Avoid boundary/wrap points of this synthetic uint8 coordinate field.
    assert np.array_equal(a[:,:,20:108],np.broadcast_to(expected[20:108],a[:,:,20:108].shape))
    print('PASS: six tilted sheets, blank-volume rejection, handedness/depth sign, affine axis sampling')
if __name__=='__main__':main()
