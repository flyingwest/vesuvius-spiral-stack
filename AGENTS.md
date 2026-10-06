# Research instructions

This repository investigates Direct Spiral-Stack reconstruction for Vesuvius CT data.

The researchers are the user and ChatGPT. The agent implements experiments,
tests hypotheses, gathers measurements, and reports evidence.

## Frozen P0

Given real Vesuvius CT data:

    transverse CT slice
        -> trace papyrus spiral
        -> parameterize by arc length s
        -> sample grayscale / ink-related signal
        -> produce one 1-D row

Stack rows across z:

    I(s,z)

P0 PASS means a coherent 2-D reconstruction is produced.

Readable text is not required.

## Important

Do not silently replace this method with the existing Vesuvius 3-D mesh/surface
unwrapping approach.

Existing Villa segmentation and spiral-fitting code may be studied for comparison,
data access, coordinate conventions, and validation.

Do not assume raw CT darkness equals ink.

Do not use OCR, Greek-language priors, or language models to repair output during P0.

Keep experiments reproducible and preserve negative results.

Prefer small falsifiable experiments over large speculative implementations.

## Research discipline

For every experiment record:

- hypothesis
- exact input volume/slice range
- parameters
- output artifacts
- quantitative metrics
- qualitative observations
- conclusion
- next recommended experiment

Never claim success solely from a visually attractive image.
