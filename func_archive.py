import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from diffraction_utils import DiffractionLine

def sensitivity_vs_strain(periodicities, theta_inc_deg, obs_angle_deg):
    strain_percent_span = 10
    sensitivities = {}

    plt.figure(figsize=(8, 6))

    for periodicity in periodicities:
        # Generate strain values
        top_period = periodicity * (1 + strain_percent_span / 100)
        strain_points = np.linspace(periodicity, top_period, strain_percent_span+1)
        wavelengths = []

        for p in strain_points:
            diff_line = DiffractionLine(grating_period=p, 
                                        inc_angle=theta_inc_deg,
                                        diff_angle=obs_angle_deg)
            wavelength = diff_line.wavelength
            wavelengths.append(wavelength)

        # Convert to strain in decimal form
        strains = (strain_points - periodicity) / periodicity  

        # Perform linear regression to get sensitivity (slope)
        slope, intercept, _, _, _ = linregress(strains, wavelengths)
        sensitivities[periodicity] = slope / 100  # Convert to nm/strain(%)

        # Plot results
        plt.scatter(strains * 100, wavelengths, label=f"P = {periodicity} nm")
        plt.plot(strains * 100, slope * strains + intercept, linestyle="--")

    # Format the title with all sensitivities
    sensitivity_text = ", ".join([f"P={p}nm: {s:.2f} nm/%" for p, s in sensitivities.items()])
    plt.title(f"Sensitivities: {sensitivity_text}")

    plt.xlabel("Strain (%)")
    plt.ylabel("Wavelength (nm)")
    plt.legend()
    plt.grid(True)
    plt.show()

    
def angle_of_orders(wavelength, n_inc, n_t, m_values, theta_inc):
    """Plots the angle of one order at specific wavelength vs periodicity"""
    periods = np.linspace(200, 1500, 1000)  # 100 points for smoother curve
    plt.figure()  # Create a new figure for plotting

    # Loop through all diffraction orders
    for m in m_values:
        angles = []  # Reset angles list for each order
        
        for period in periods:
            angle = DiffractionLine(wavelength, period, theta_inc, m=m).diff_angle
            angles.append(angle)
        
        # Plot each curve with label for the diffraction order
        plt.plot(periods, angles, label=f'm = {m}')
    
    # Plotting the angles versus period
    plt.xlabel('Period (nm)')
    plt.ylabel('Angle (degrees)')
    title = f'Angle vs Period for Grating Equation\nWavelength: {wavelength} nm, Incidence Angle: {theta_inc}°'
    plt.title(title)
    plt.legend()  # Add a legend to distinguish diffraction orders
    plt.grid(True)

