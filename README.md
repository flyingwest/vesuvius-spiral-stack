# Vesuvius Direct Spiral-Stack

Experimental direct reconstruction method for Herculaneum CT volumes.

## P0 hypothesis

For each transverse CT slice:

1. identify the papyrus spiral,
2. parameterize the spiral by arc length,
3. sample a grayscale signal along the spiral,
4. straighten it into one 1-D row,
5. stack consecutive rows along the scroll axis.

The result is a 2-D raster I(s,z).

## P0 success criterion

P0 PASS requires a spatially coherent 2-D reconstruction from real Vesuvius CT data.

Readable letters are NOT required for P0.

Systematic stretch, shear, or curvature is acceptable.

## P0 failure criterion

STOP or reconsider the approach if reliable winding correspondence cannot be maintained
between neighboring transverse CT slices without effectively reconstructing the full
3-D surface manifold.

## Explicitly out of scope for P0

- OCR
- Greek language reconstruction
- synthetic CT
- physical surrogate scrolls
- learned distortion correction
- full X-ray simulation
- reproducing the existing Vesuvius mesh pipeline

## Idea B

Deferred until after P0:

Known flat document -> roll/deform -> CT or CT simulation -> known ground truth ->
learned distortion correction / reconstruction.
