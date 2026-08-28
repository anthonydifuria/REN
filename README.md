# REN

REN (Radial Exponential Normalised) is a minimal mathematical framework for generating complex radial shapes using a simple, normalized polar equation.

Core Idea

The REN formula is defined as:

r(\theta) = \exp\big(-\lambda \, |\sin(k\theta)|^q \big)

It produces structured, flower-like or star-like geometries with very few parameters, while remaining fully normalized (maximum radius is always 1).

Key Properties

* Normalized: the maximum radius is always 1
* Minimal: only 3 core parameters
* Expressive: capable of generating a wide range of organic and geometric shapes
* Symmetric by design: evenly distributed radial features

Parameters

* k — controls the number of radial peaks
* λ (lambda) — controls depth (how much the shape contracts)
* q — controls sharpness (smooth vs pointed features)

Extensions

REN can be extended by introducing an angular transformation:

r(\theta) = \exp\big(-\lambda \, |\sin(k\,\phi(\theta))|^q \big)

This allows:

* local clustering of features
* controlled asymmetry
* smooth redistribution of angular space

Purpose

REN is designed as a compact alternative to more complex parametric shape formulas, focusing on clarity, control, and visual richness with minimal structure.

Documentation

The full write-up (construction, properties, extensions, and the application of REN inside Hutchinson Synthesis) is in `latex/ren.pdf` (source: `latex/ren.tex`).

⸻

This repository contains implementations, experiments, and visual explorations of the REN framework.
