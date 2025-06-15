import numpy as np
import matplotlib.pyplot as plt

# === Parameters ===
N = 6                     # Number of slits
slit_width_ratio = 2.0    # a/λ
period_ratio = 5.0       # d/λ

# === sin(θ) range ===
sin_theta = np.linspace(-0.98, 0.98, 10000)  # sin(θ) from -1 to 1 (physical limit)

# === Compute β and α ===
beta = np.pi * slit_width_ratio * sin_theta
alpha = np.pi * period_ratio * sin_theta

# === Diffraction term: (sin β / β)^2 with β=0 handling ===
diffraction = np.ones_like(beta)
nonzero_beta = beta != 0
diffraction[nonzero_beta] = (np.sin(beta[nonzero_beta]) / beta[nonzero_beta]) ** 2

# === Interference term: (sin Nα / sin α)^2 with α=0 handling ===
interference = np.ones_like(alpha)
nonzero_alpha = alpha != 0
interference[nonzero_alpha] = (np.sin(N * alpha[nonzero_alpha]) / np.sin(alpha[nonzero_alpha])) ** 2

# === Total intensity ===
total_intensity = diffraction * interference / (N ** 2)

# rainbow colours
colours = plt.cm.rainbow(np.linspace(0, 1, 10))

# === Plotting ===
plt.figure(figsize=(12, 6))
plt.plot(sin_theta, diffraction,        label='Diffraction Envelope', linestyle='--', linewidth=2, color=colours[0])
plt.plot(sin_theta, interference / N**2, label='Interference (scaled)', linestyle='-.', linewidth=2, color=colours[7])
plt.plot(sin_theta, total_intensity,     label='Total Intensity', linewidth=2, color=colours[9])
plt.xlabel(r'$\sin(\theta)$', fontsize=16)
plt.ylabel(r'$I/I_0$ (normalized)', fontsize=16)
# plt.title(r'Diffraction Pattern vs. $\sin(\theta)$', fontsize=16)
plt.legend(fontsize=14)
plt.grid(False, alpha=0.3)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.show()
