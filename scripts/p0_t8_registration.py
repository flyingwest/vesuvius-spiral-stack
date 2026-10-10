"""Registration checks without CT-versus-ink feature tests."""
from p0_t8_acquire import *
from PIL import Image
from scipy.ndimage import distance_transform_edt
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 ct=np.load(DATA/'ct_profiles_dyx.npy');g=np.load(DATA/'grid_geometry.npy')
 images={name:np.array(Image.open(DATA/f'{name}.png')) for name in ['inklabels','mask','ir']}
 assert all(a.shape[:2]==(7606,5249) for a in images.values())
 labels=images['inklabels']>0;mask=images['mask']>0;ir=images['ir'];ir=ir.mean(2) if ir.ndim==3 else ir
 crop=(slice(2784,4832),slice(1600,3648));lab=labels[crop];m=mask[crop];photo=ir[crop]
 assert set(np.unique(images['inklabels'])).issubset({0,1,255});assert set(np.unique(images['mask'])).issubset({0,1,255})
 d=distance_transform_edt(lab);dn=distance_transform_edt(~lab);dm=distance_transform_edt(m)
 usable=(dm>16)&((d>2)|(dn>2));norm=np.linalg.norm(g[:,:,3:],axis=2)
 result={'global_shape':[7606,5249],'crop_shape':list(lab.shape),'same_xy_registration':'official shared surface-volume/inklabels/ir/mask coordinates, no flips/rotations/crop offsets fitted','label_positive_pixels':int(labels.sum()),'positive_outside_mask':int((labels&~mask).sum()),'crop_positive_pixels':int((lab&m).sum()),'crop_mask_pixels':int(m.sum()),'usable_grid_samples':int(usable[8::16,8::16].sum()),'normal_norm_range_on_usable_grid':[float(norm[usable[8::16,8::16]].min()),float(norm[usable[8::16,8::16]].max())],'ir_median_positive':float(np.median(photo[lab&m])),'ir_median_negative':float(np.median(photo[~lab&m])),'ct_nonzero_fraction_inside_mask':float(np.mean(ct[:,m]>0)),'registered_depth':32,'depth_spacing_um':3.24,'physical_ink_depth':'not independently established; fixed nominal surface, not chosen using labels','subpixel_registration_error':'not measured; two-pixel class-boundary exclusion does not guarantee complete removal of original alignment error'}
 assert result['usable_grid_samples']>1000
 assert np.all(np.isfinite(g)) and np.all(norm[usable[8::16,8::16]]>0.99)
 fig,ax=plt.subplots(2,3,figsize=(15,10))
 ax[0,0].imshow(photo,cmap='gray');ax[0,0].set_title('Registered IR (not CT input)')
 ax[0,1].imshow(photo,cmap='gray');ax[0,1].contour(lab,levels=[.5],colors='red',linewidths=.35);ax[0,1].set_title('Official ink-label boundary on IR')
 ax[0,2].imshow(ct[3],cmap='gray',vmin=0,vmax=65535);ax[0,2].set_title('CT layer32: same XY, no transform')
 ax[1,0].imshow(m,cmap='gray');ax[1,0].imshow(np.ma.masked_where(~lab,lab),cmap='autumn',alpha=.5);ax[1,0].set_title('Papyrus validity mask + ink labels')
 ax[1,1].imshow(photo,cmap='gray');ax[1,1].axhline(1024,color='cyan');ax[1,1].set_title('Top: direction only / Bottom: confirmation')
 ax[1,2].imshow(g[:,:,4],cmap='coolwarm');ax[1,2].set_title('Original surface normal Y at grid')
 for a in ax.ravel():a.axis('off')
 fig.tight_layout();fig.savefig(OUT/'registration_diagnostic.png',dpi=160);plt.close(fig)
 np.savez_compressed(OUT/'registration_arrays.npz',labels=lab,mask=m,ir=photo,usable=usable)
 dump(OUT/'registration.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
