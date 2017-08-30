
def P_shade(rr,theta,vw=0,theta_w_deg=0,L_mean=3,L_std=1,fs=0.5,fr=2,stiff=1):

    # Convert theta_w to radians
    theta_w = theta_w_deg*np.pi/180

    # Frond half-span angle
    alpha = np.arctan((1+fs)/(2*fs*fr))

    # Try single integral
    def single_pdf(theta_f):

        theta_prime = theta - theta_f + np.pi/2
        S = np.sign(np.pi/2-theta_prime)

        # Minimum L for shading as a function of central frond angle
        L_star = rr*(np.sin(theta_prime)+2*S*fr/(1+fs)*np.cos(theta_prime))

        # Integrated normalized L distribution
        C_L = erf((L_star-L_mean)/(np.sqrt(2)*L_std)) / (1-erf(-L_mean/(np.sqrt(2)*L_std)))

        # Theta_f distribution
        if vw > 0:
            P_theta_f = vonmises.pdf(theta_f,vw/stiff,theta_w)
        else:
            P_theta_f = np.ones_like(theta_f)/(2*np.pi)

        return (1-C_L)*P_theta_f/2
    SI = quad(single_pdf,theta-alpha,theta+alpha)

    return SI[0]

# Spatial resolution
ds = 1e-1
# Angular resolution
#da = pi/6

# Define domain
xmin, xmax, dx = -1, 1, ds
ymin, ymax, dy = -1, 1, ds
zmin, zmax, dz =  0, 1, ds
#thmin, thmax, nth, = 0, 2*np.pi, da
#phmin, phmax, nph = 0, np.pi, da

nx = int((xmax - xmin) / dx)
ny = int((ymax - ymin) / dy)
nz = int((zmax - zmin) / dz)

# Allocate arrays
x = np.linspace(xmin, xmax, nx, endpoint=False)
y = np.linspace(ymin, ymax, ny, endpoint=False)
z = np.linspace(zmin, zmax, nz, endpoint=False)
#th = np.linspace(thmin, thmax, nth)
#ph = np.linspace(phmin, phmax, nph)

# Convert x and y to polar coordinates
xg, yg = np.meshgrid(x, y, indexing='ij')
rr = np.sqrt(xg**2 + yg**2)
th = np.arctan2(xg, yg)


#= Calculate P_shade for every point =#

P3D = np.zeros([nx, ny, nz])

for ii in range(nx):
    print("{}/{}".format(ii+1, nx))
    for jj in range(ny):
        for kk in range(nz):
            P3D[ii, jj, kk] = P_shade(
                rr=rr[ii,jj],
                theta=th[ii,jj],
                vw=values['vw'][kk],
                theta_w_deg=values['theta_w_deg'][kk],
                L_mean=values['L_mean'][kk],
                L_std=values['L_std'][kk]
            )

