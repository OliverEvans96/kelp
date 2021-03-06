# *swdlyt*

It's a light model for seaweed!

Numerically solve the Radiative Transfer Equation in Fortran 90 in order to model the light field in a vertical-line kelp aquaculture environment.

This is the code used in my Master's Thesis in applied math.
The full thesis can be found at https://github.com/OliverEvans96/msthesis

Abstract:

> A mathematical model is developed to describe the light field in
> vertical line seaweed cultivation to determine the
> degree to which the seaweed shades itself and limits the
> amount of light available for photosynthesis.
> A probabilistic description of the spatial distribution of kelp
> is formulated using simplifying assumptions about frond geometry and orientation.
> An integro-partial differential equation called the  radiative transfer equation
> is used to describe the light field as a function of position and angle.
> A finite difference solution is implemented, providing robustness and accuracy
> at the cost of large CPU and memory requirements, and
> a less computationally intensive asymptotic approximation is explored for the case of low
> scattering.
> Conditions for applicability of the asymptotic approximation are discussed,
> and depth-dependent light availability is compared to the predictions of simpler light models.
> The 3D model of this thesis is found to predict significantly lower light levels than the simpler 1D models,
> especially in regions of high kelp density where a precise description of self-shading is most important.
