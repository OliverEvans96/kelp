# Fortran 3D Kelp Distribution Outline

Given properties at each depth:
- Length mean
- Length std
- Water current angle
- Water current velocity

Produce:
- Probability of kelp at each point on 3D grid

## Procedure

- Calculate $P_{2D}$
    - Calculate length distribution
        - PDF
            - numerical integration
        - CDF
            - `erf` (or numerical integration)
    - Calculate angle distribution
        - PDF
            - von Mises distribution
                - Modified Bessel function of the first kind of order 0 ($I_0$)
        - CDF
            - numerical integration
- Calculate $R_s(\theta_p, \theta_r)$
    - $\theta_p - \alpha < \theta_f + \theta_p + \alpha$
        - Calculate $\alpha$ from $f_s, f_r$
    - $L > L_{min}(\theta_p, r_p)$
        - $L_min(\theta_p, r_p)$
            - $r_f(\theta)$
                - $r_f'(\theta')$
                    - $S(\theta')$
- Integrate $P_{2D}$ over $R_s$
    - Polar $\leftrightarrow$ Cartesian conversion 
    - numerical integration over non-rectangular region