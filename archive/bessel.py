# Plot first kind bessel function J(x) whwre I can define x
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jn

def plot_bessel_function(h_values, n_values, title='Bessel Functions of the First Kind'):
    """
    Plot Bessel functions of the first kind for given x values and orders.

    Parameters:
    - x_values: array-like, values at which to evaluate the Bessel functions.
    - n_values: list, orders of the Bessel functions to plot.
    - title: str, title of the plot.
    """
    plt.figure(figsize=(10, 6))
    


    for n in n_values:
        
        theta_m = np.arcsin(n*wl/P)
        print(f"Order: {n}, Theta_m: {(theta_m)} degrees and 1 + np.cos(theta_m): {1 + np.cos(theta_m)}")
        x_values = h_values * (np.pi / wl) * (np.cos(theta_m)) # Convert h to x for Bessel function
        y_values = jn(n, x_values)**2 
        plt.plot(h_values, y_values, label=f'J_{n}(x)')
    
    plt.title(title)
    plt.xlabel('h [nm]')
    plt.ylabel('J_n^2(x)')
    plt.axhline(0, color='black', lw=0.5, ls='--')
    plt.axvline(0, color='black', lw=0.5, ls='--')
    plt.grid()
    plt.legend()
    plt.show()

# Example usage
if __name__ == "__main__":
    wl = 532
    h_values = np.linspace(0, 532, 1000)  # Range of x values
    n_values = [0, 1, 2]           # Orders of Bessel functions to plot
    P = 2100
    plot_bessel_function(h_values, n_values)


