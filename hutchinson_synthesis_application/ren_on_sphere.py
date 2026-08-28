"""
REN on the sphere -- spatial application used in Hutchinson Synthesis.

This module is a clean, standalone extraction of the mechanism that applies
the REN formula (see ../latex/ren.tex) to constrain the spatial distribution
of sound events on the sphere S^2, as implemented in:

    https://github.com/anthonydifuria/hutchinson-synthesis

It is derived from hutchinson_pure/main_v1_independent.py in that
repository (functions `_campiona_ren_2d` and the REN
block of `_aggiunge_derived`), simplified and decoupled from the rest of the
synthesiser's configuration/state so it can be read and reused on its own.

Pipeline
--------
1. `sample_ren_centers` places one REN "flower" center per acoustic niche on
   the sphere, spreading them apart with a repulsive relaxation.
2. `ren_boundary` evaluates the (locally projected) REN boundary around a
   given center, using the same base formula as the 1D curve:
       r(phi) = exp(-lambda * |sin(k * phi / 2)|^q)
   (the factor 2 in the angle adapts the lobe count to the tangent-plane
   projection used here -- see "REN in Hutchinson Synthesis" in
   ../latex/ren.tex).
3. `apply_ren_magnet` pulls a candidate (azimuth, elevation) position toward
   the nearest niche's REN lobes, with a `force` in [0, 1] controlling how
   strongly the constraint is enforced (0 = no effect, 1 = the point is
   guaranteed to end up inside the lobes, via `clamp_inside_ren`).

This file has no dependency on the rest of either repository -- it only
needs NumPy.
"""

import numpy as np


# ---------------------------------------------------------------------------
# 1. Placing REN centers on the sphere
# ---------------------------------------------------------------------------

def sample_ren_centers(n_centers, rng, relaxation_iterations=200, step=0.3):
    """Sample `n_centers` points on S^2, spread apart by mutual repulsion.

    Each point starts at a random position on the sphere and is iteratively
    pushed away from the others, projected back onto the tangent plane and
    renormalised, so the final centers are approximately equidistant.

    Returns a list of (azimuth, elevation) pairs, in radians.
    """
    points = []
    for _ in range(n_centers):
        u = rng.uniform(-1, 1)
        el = float(np.arcsin(u))
        az = float(rng.uniform(-np.pi, np.pi))
        points.append([el, az])

    for _ in range(relaxation_iterations):
        xyz = [_to_cartesian(el, az) for el, az in points]
        new_points = []
        for i in range(n_centers):
            xi, yi, zi = xyz[i]
            fx = fy = fz = 0.0
            for j in range(n_centers):
                if i == j:
                    continue
                dx, dy, dz = xi - xyz[j][0], yi - xyz[j][1], zi - xyz[j][2]
                dist2 = dx**2 + dy**2 + dz**2 + 1e-8
                fx += dx / dist2
                fy += dy / dist2
                fz += dz / dist2
            # Project the repulsive force onto the tangent plane
            dot = fx * xi + fy * yi + fz * zi
            fx -= dot * xi
            fy -= dot * yi
            fz -= dot * zi
            nx, ny, nz = xi + step * fx, yi + step * fy, zi + step * fz
            norm = np.sqrt(nx**2 + ny**2 + nz**2) or 1.0
            nx, ny, nz = nx / norm, ny / norm, nz / norm
            new_points.append([
                float(np.arcsin(np.clip(nz, -1.0, 1.0))),
                float(np.arctan2(ny, nx)),
            ])
        points = new_points
        step *= 0.98

    return [(float(np.clip(az, -np.pi, np.pi)),
             float(np.clip(el, -np.pi / 2, np.pi / 2)))
            for el, az in points]


def _to_cartesian(el, az):
    return (np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el))


# ---------------------------------------------------------------------------
# 2. REN boundary around a center (local gnomonic projection)
# ---------------------------------------------------------------------------

def _gnomonic_basis(center_az, center_el):
    """Local tangent-plane basis at (center_az, center_el)."""
    cos_el_c, sin_el_c = np.cos(center_el), np.sin(center_el)
    cos_az_c, sin_az_c = np.cos(center_az), np.sin(center_az)
    center = np.array([cos_el_c * cos_az_c, cos_el_c * sin_az_c, sin_el_c])
    e_az = np.array([-sin_az_c, cos_az_c, 0.0])
    e_el = np.array([-sin_el_c * cos_az_c, -sin_el_c * sin_az_c, cos_el_c])
    return center, e_az, e_el


def _project(az, el, center, e_az, e_el):
    """Project (az, el) onto the tangent plane at `center`.

    Returns the local polar coordinates (radius, angle) in that plane.
    """
    point = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    vec = point - center
    dx = float(np.dot(vec, e_az))
    dy = float(np.dot(vec, e_el))
    r = np.sqrt(dx * dx + dy * dy)
    phi = float(np.arctan2(dy, dx)) if r > 1e-8 else 0.0
    return r, phi


def _unproject(r, phi, center, e_az, e_el):
    """Map a local (radius, angle) back onto the sphere, as (az, el)."""
    offset = r * np.cos(phi) * e_az + r * np.sin(phi) * e_el
    point = center + offset
    norm = np.linalg.norm(point)
    if norm < 1e-10:
        return 0.0, 0.0
    point = point / norm
    return (float(np.arctan2(point[1], point[0])),
            float(np.arcsin(np.clip(point[2], -1.0, 1.0))))


def ren_boundary(phi, k=5, lam=2.0, q=2.0):
    """Local REN boundary value at angle `phi`, in the tangent plane.

    r(phi) = exp(-lambda * |sin(k * phi / 2)|^q)

    Same base shape as the 1D REN curve (eq. 1 in ../latex/ren.tex); the
    angle is halved to adapt the lobe count to this local projection.
    """
    return float(np.exp(-lam * np.abs(np.sin(k * phi / 2.0)) ** q))


# ---------------------------------------------------------------------------
# 3. Applying the REN constraint to a point ("magnet" + clamp)
# ---------------------------------------------------------------------------

def apply_ren_magnet(az, el, center_az, center_el, radius,
                      force=1.0, k=5, lam=2.0, q=2.0):
    """Pull (az, el) toward the REN lobes of the niche centered at
    (center_az, center_el), with a given `radius` and `force` in [0, 1].

    force = 0.0  -> the point is left untouched.
    force = 1.0  -> the point is pulled fully onto the target boundary
                     (combine with `clamp_inside_ren` for a hard guarantee).
    """
    center, e_az, e_el = _gnomonic_basis(center_az, center_el)
    r_curr, phi = _project(az, el, center, e_az, e_el)
    if r_curr < 1e-8:
        return az, el

    r_ren = radius * ren_boundary(phi, k=k, lam=lam, q=q)
    r_target = min(r_curr, r_ren)
    r_final = r_curr * (1.0 - force) + r_target * force
    return _unproject(r_final, phi, center, e_az, e_el)


def clamp_inside_ren(az, el, center_az, center_el, radius, k=5, lam=2.0, q=2.0):
    """Guarantee (az, el) falls strictly inside the REN lobes, projecting it
    onto the nearest lobe border if it currently falls outside."""
    center, e_az, e_el = _gnomonic_basis(center_az, center_el)
    r_v, phi_v = _project(az, el, center, e_az, e_el)
    r_ren = radius * ren_boundary(phi_v, k=k, lam=lam, q=q)
    if r_v <= r_ren + 1e-6:
        return az, el  # already inside

    lobe_angles = [2.0 * np.pi * n / k for n in range(k)]
    phi_near = min(lobe_angles, key=lambda p: abs(
        np.arctan2(np.sin(phi_v - p), np.cos(phi_v - p))))
    r_lobe = radius * ren_boundary(phi_near, k=k, lam=lam, q=q)
    r_final = min(r_v, r_lobe * 0.99)
    return _unproject(r_final, phi_near, center, e_az, e_el)


if __name__ == "__main__":
    # Minimal usage example: 3 niches, one candidate event pulled toward
    # the nearest niche's REN lobes.
    rng = np.random.default_rng(0)
    centers = sample_ren_centers(n_centers=3, rng=rng)
    print("Niche centers (az, el) in radians:", centers)

    az_c, el_c = centers[0]
    candidate_az, candidate_el = az_c + 0.3, el_c + 0.1
    pulled_az, pulled_el = apply_ren_magnet(
        candidate_az, candidate_el, az_c, el_c, radius=0.4, force=1.0
    )
    print(f"Candidate ({candidate_az:.3f}, {candidate_el:.3f}) "
          f"-> pulled to ({pulled_az:.3f}, {pulled_el:.3f})")
