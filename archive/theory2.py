import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Parameters
N = 2            # Number of slits
slit_width = 0.2  # Relative slit width (arbitrary units)
d = 1.0           # Relative slit spacing (same units)
llambda = 0.1      # Wavelength (same units as slit_width and d)

# Angular range in degrees
theta_deg = np.linspace(-35, 35, 10000)
theta_rad = np.deg2rad(theta_deg)

# Compute beta and alpha
beta  = np.pi * slit_width * np.sin(theta_rad) / llambda
alpha = np.pi * d / llambda * np.sin(theta_rad)

# Diffraction term: (sin β / β)^2 with handling at β=0
diffraction = np.ones_like(beta)
mask_beta = beta != 0
diffraction[mask_beta] = (np.sin(beta[mask_beta]) / beta[mask_beta])**2

# Interference term: (sin Nα / sin α)^2 with handling at α=0
interference = np.ones_like(alpha)
mask_alpha = alpha != 0
interference[mask_alpha] = (np.sin(N*alpha[mask_alpha]) / np.sin(alpha[mask_alpha]))**2

# Total normalized intensity (max = 1)
total_intensity = diffraction * interference / (N**2)
# rainbow colours
colours = plt.cm.rainbow(np.linspace(0, 1, 10))

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(theta_deg, diffraction,     label='Diffraction Envelope', color=colours[0],   linewidth=2, linestyle='--')
plt.plot(theta_deg, interference/N**2,label='Interference (scaled)', color=colours[7], linewidth=2, linestyle='-.')
plt.plot(theta_deg, total_intensity, label='Total Intensity',       color=colours[8],   linewidth=2)
plt.xlabel(r'$\theta$ (°)', fontsize=16)
plt.ylabel(r'$I/I_0$ (normalized)', fontsize=16)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()
