# REN in Hutchinson Synthesis

This folder documents, and re-implements as a standalone module, how the
REN formula (see `../latex/ren.tex`) is used inside **Hutchinson
Synthesis** to shape the spatial distribution of sound events on the
sphere S^2:

> https://github.com/anthonydifuria/hutchinson-synthesis

## What's here

`ren_on_sphere.py` -- a clean, dependency-free extraction of the mechanism
implemented in `hutchinson_pure/main_v1_independent.py` of the
hutchinson-synthesis repository (functions `_campiona_ren_2d` and the REN
block inside `_aggiunge_derived`). It reproduces the same three steps:

1. **Placing centers.** One REN "flower" center per acoustic niche is
   sampled on the sphere and spread apart via repulsive relaxation
   (`sample_ren_centers`).
2. **Local boundary.** Around each center, the REN curve is evaluated in a
   local tangent-plane (gnomonic) projection (`ren_boundary`) -- the same
   base shape as the 1D formula in `../latex/ren.tex`, adapted to the
   projected geometry.
3. **Magnet + clamp.** Candidate event positions are pulled toward the
   nearest niche's REN lobes with a `force` parameter in [0, 1]
   (`apply_ren_magnet`); at `force = 1.0`, `clamp_inside_ren` guarantees
   the final position falls inside the lobes.

This is the mechanism behind the REN spatial-distribution audio and
visual examples in the `examples/` folder of that repository.

## Relationship to the base REN formula

The 1D formula in `../latex/ren.tex` defines a boundary curve over a full
angular range `[0, 2*pi)`. In Hutchinson Synthesis it is applied locally,
once per niche, to the tangent plane at each niche's center on the sphere
-- not globally over the sphere itself. This is why the local
implementation here uses `k * phi / 2` instead of `k * phi`: the angle in
the tangent-plane projection is halved to keep the intended number of
lobes after projection.

Note: the Hutchinson Synthesis codebase also contains a separate,
exploratory function (`ren_2d`) combining independent azimuth/elevation
lobe counts multiplicatively. It is not called anywhere in the generation
pipeline and is not the mechanism reproduced here.

## Usage

```bash
python ren_on_sphere.py
```

runs a minimal example: it samples 3 niche centers and pulls one candidate
point toward the nearest niche's lobes.

This module has no dependency on the rest of either repository -- it only
needs NumPy.
