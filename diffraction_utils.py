import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.stats import linregress
from colour_utils import *




#----------------------------------
# -------- BASIC FUNCTIONS --------
#----------------------------------

from archive.wavelen2rgb import wavelen2rgb

def set_dark_mode(dark_mode=True):
    """Configures Matplotlib for dark mode styling if enabled."""
    if dark_mode:
        plt.rcParams.update({
            'axes.facecolor': '#2E2E2E',
            'axes.edgecolor': 'grey',
            'figure.facecolor': '#2E2E2E',
            'figure.edgecolor': 'white',
            'axes.labelcolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'grid.color': 'gray',
            'text.color': 'white',
            'legend.facecolor': '#2E2E2E',
            'legend.edgecolor': 'white'
        })

#----------------------------------
# ------------ CLASS  -------------
#----------------------------------

class DiffractionLine:
    
    def __init__(self, wavelength=None, grating_period=None, inc_angle=None, diff_angle=None, m=1):
        self.wavelength = wavelength  # λ (nm)
        self.grating_period = grating_period  # d (m)
        self.inc_angle = inc_angle  # θ_inc (degrees)
        self.diff_angle = diff_angle  # θ_obs (degrees)
        self.order = m  # Diffraction order
        
        # Automatically compute missing attributes
        if self.wavelength is None:
            self.wavelength = - (grating_period/m) * (np.sin(np.radians(inc_angle)) - np.sin(np.radians(diff_angle))) 
        if self.diff_angle is None:
            sine_term = (np.sin(np.radians(inc_angle)) - m * (wavelength / grating_period))
            if np.abs(sine_term) > 1:
                if m < 0:
                    self.diff_angle = 90
                else:
                    self.diff_angle = -90
            else:
                self.diff_angle = np.degrees(np.arcsin(sine_term))

    @property
    def x_proj(self):
        return - np.cos(self.diff_angle + 90)


    def __repr__(self):
        return (f"Grating(wavelength={self.wavelength}, grating_period={self.grating_period}, "
                f"diffraction_angle={self.diff_angle}, order={self.order})")

#---------------------------------
# -------- PROGRAMS --------
#---------------------------------

def rotatinator_view(period, inc_angle=0, 
                     angle_res=0.5, 
                     wl_range=[400, 720], 
                     angle_range=[0, 90], 
                     scale_height=1, 
                     m_values=[1, 2, 3, 4],
                     x_grid = False,
                     h_sep_line = False, 
                     dark_mode=True):
    
    # Calculate wavelength at lines in plot for each order
    angle_vals = np.arange(angle_range[0], angle_range[1]+1, angle_res)
    diffraction_data = []  # Store (order, angle, wavelength, color)

    for m in m_values:
        for diff_angle in angle_vals:
            diffraction_line = DiffractionLine(grating_period=period, inc_angle=inc_angle, m=m, diff_angle=diff_angle)
            wl = diffraction_line.wavelength
            
            if wl_range[0] <= wl <= wl_range[1]:
                color = rgb_from_wavelength(wl) if m != 0 else (0.5, 0.5, 0.5)
                diffraction_data.append((m, diff_angle, wl, color))
    
    # Initalize dark mode plot
    set_dark_mode(dark_mode)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 3*scale_height))
    plt.subplots_adjust(left=0.1, right=0.9, top=0.95*scale_height, bottom=0.15*(1/scale_height))
    
    height = 0.85

    for m, diff_angle, wl, color in diffraction_data:
        y_offset = -m  # Stagger orders vertically
        rect = plt.Rectangle((diff_angle - angle_res/2, y_offset - height/2), angle_res, height, color=color, alpha=1)
        ax.add_patch(rect)
    
    # Add mixing line below the last order
    # The mixing line will be drawn at a y-position below -max(m_values)
    mixing_y = -max(m_values) - 1.0   # y center of the mixing row
    mixing_height = height * 0.8      # mixing row height (adjustable)

    # For each diffraction angle, determine the color to plot:
    # - If there are two overlapping orders, mix them.
    # - If there is one order, use its color.
    # - If there are no orders at that angle, skip drawing.
    for diff_angle in angle_vals:
        overlapping = [d for d in diffraction_data if abs(d[1] - diff_angle) < 1e-6]
        if len(overlapping) == 0:
            # Leave blank if no orders are present at this angle.
            continue
        elif len(overlapping) == 1:
            mix_color = overlapping[0][3]
        elif len(overlapping) == 2:
            wl1 = overlapping[0][2]
            wl2 = overlapping[1][2]
            mix_color = mix_two_wavelengths(wl1, wl2)
        else:
            # Fallback: if more than two orders overlap, just take the first two.
            wl1 = overlapping[0][2]
            wl2 = overlapping[1][2]
            mix_color = mix_two_wavelengths(wl1, wl2)
        
        rect = plt.Rectangle((diff_angle - angle_res/2, mixing_y - mixing_height/2), angle_res, mixing_height, color=mix_color, alpha=1)
        ax.add_patch(rect)
    
    # Fontsize
    plt.rcParams.update({'font.size': 8})
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.xaxis.label.set_size(8)
    ax.yaxis.label.set_size(8)

    # Grid and formatting
    ax.set_xlim(angle_range[0], angle_range[1])
    # Extend the y-limit to include the mixing line.
    ax.set_ylim(mixing_y - mixing_height, -0.5)
    ax.set_yticks([-m for m in m_values])
    ax.set_yticklabels([f'{m}' for m in m_values])
    ax.set_xlabel("Diffraction Angle (degrees)")
    ax.set_ylabel("Order")
    title_line_1 = f"Diffraction Orders for Grating Period = {period:.0f} nm, Wavelength Range: {wl_range[0]} - {wl_range[1]} nm"
    title_line_2 = f"\nIncident Angle = {inc_angle:.1f}°, Angle Resolution = {angle_res}°"
    ax.set_title(title_line_1 + title_line_2)
    if x_grid == True:
        ax.grid(axis='x', linestyle='--', alpha=0.5, color='gray')
    
    # Horizontal separation lines
    if h_sep_line == True:
        for m in range(min(m_values) - 1, max(m_values) + 2):
            ax.axhline(-m + 0.5, color='gray', linewidth=2*scale_height)

    plt.show(block=True)



def objective_view(period, 
                inc_angle, 
                obs_angle, 
                wl_range=[400, 720], 
                wl_step=50, 
                na=0.25, 
                m_values=[-3, -2, -1, 0, 1, 2, 3],
                dark_mode=False):
    
    wavelengths = np.arange(wl_range[0], wl_range[1] + wl_step, wl_step)
    x_proj_data = []  # Store (x_proj, wavelength, color)

    for m in m_values:
        for wl in wavelengths:
            diffraction_line = DiffractionLine(wl, period, inc_angle, m=m)
            diff_angle = diffraction_line.diff_angle
            
            # Shift the angle by obs angle (centers the objective at 0 deg) and add 90 deg (for standard polar plot)
            theta_m_adj_rad = np.radians(diff_angle - obs_angle + 90)  
            theta_na_rad = np.pi/2 - np.arcsin(na)
            x_proj = np.sin(theta_na_rad)/np.tan(theta_m_adj_rad)  # x-projection to objective coordinate system

            if np.isfinite(x_proj):  # Avoid plotting invalid values
                if m == 0:
                    color = (0.5, 0.5, 0.5)  # Grey color for m=0
                else:
                    color = rgb_from_wavelength(wl)
                x_proj_data.append((x_proj, wl, color, diff_angle))

    # Initalize dark mode plot
    set_dark_mode(dark_mode)

    # Plotting
    fig, ax = plt.subplots(figsize=(8, 3))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.85, bottom=0.15)  # Adjust these values as needed
    
    for x_proj, wl, color, diff_angle in x_proj_data:
        ax.add_patch(plt.Rectangle((x_proj - 0.005, 0), 0.01, 1, color=color, alpha=0.6))
        # Add diff_angle on top x-axis here for every first rectangle in order. ->
        if wl%100 == 0:
            ax.text(x_proj, 0, f'{diff_angle:.0f}°', color='black', fontsize=8, ha='center', va='bottom')
    
    # Mark NA limits
    na_x_proj = np.cos((np.pi/2 - np.arcsin(na)))
    ax.axvline(na_x_proj, linewidth=4, color='black', linestyle='-', label=f'NA {na}')
    ax.axvline(-na_x_proj, linewidth=4, color='black', linestyle='-')
  
    ax.set_xlim(-1, 1)
    ax.set_yticks([])
    ax.set_xlabel("Projection onto objective x-coordinate / True angle of diffraction")
    title_line_1 = (f'Diffraction Orders for grating period = {period:.0f} nm, incident lighte angle = {inc_angle:.1f}°')
    title_line_2 = (f'Observance angle = {obs_angle:.1f}°, Objective NA = {na}')
    title_line_3 = (f'Wavelengths shown: [{wl_range}]')
    ax.set_title(title_line_1 + '\n' + title_line_2 + '\n' + title_line_3)
    ax.legend()
    
    plt.show(block=False)      

def polar_orders_overview(period, 
                inc_angle, 
                na, 
                obs_angle,
                m_values=[-4, -3, -2, -1, 0, 1, 2, 3, 4],
                wavelength_start=400, 
                wavelength_stop=680,
                wavelength_interval=100,
                dark_mode=False):
    
    wavelengths = np.arange(wavelength_start, wavelength_stop+wavelength_interval, wavelength_interval) 
    
    # Initalize dark mode plot
    set_dark_mode(dark_mode)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    plt.subplots_adjust(bottom=0)  # Adjust bottom space for two sliders

    ax.set_ylim(0, 1)  # Set radius limit
    ax.set_yticklabels([])  # Hide radial labels

    ax.set_xticks([0, np.pi/6, np.pi/3, np.pi/2, 2*np.pi/3, 5*np.pi/6, np.pi])
    ax.set_xticklabels(['-90°', '-60°', '-30°', '0°', '30°', '60°', '90°'])

    def update(val):
        ax.clear()  # Clear old plot
        ax.set_ylim(0, 1)
        ax.set_yticklabels([])
        ax.set_xticks([0, np.pi/6, np.pi/3, np.pi/2, 2*np.pi/3, 5*np.pi/6, np.pi])
        ax.set_xticklabels(['-90°', '-60°', '-30°', '0°', '30°', '60°', '90°'])

        new_period = slider_period.val  # Get current period value
        new_theta_inc = slider_theta.val  # Get current incident angle

        for wavelength in wavelengths:
            for m in m_values:
                sine_term = (np.sin(np.radians(new_theta_inc)) - m * (wavelength / new_period))
                
                if -1 <= sine_term <= 1:
                    theta_m_rad = np.arcsin(sine_term)
                    theta_polar = np.pi/2 - theta_m_rad

                    color = rgb_from_wavelength(wavelength) if m != 0 else (0, 0, 0)
                    ax.arrow(theta_polar, 0, 0, 0.8, head_width=0.02, head_length=0.02, fc=color, ec=color, alpha=0.7)
                    ax.text(theta_polar, 0.85, f'{m}', color=color, fontsize=8, ha='center', va='baseline')

        # Draw NA limit
        na_angle = np.arcsin(na) 
        for sign in [-1, 1]:
            theta_na = np.pi/2 - sign * na_angle + np.radians(obs_angle)
            ax.arrow(theta_na, 0, 0, 1.2, head_width=0, head_length=0, fc='grey', ec='grey', alpha=0.7, linewidth=3)
            ax.text(theta_na, 1, 'NA', color='black', fontsize=8, ha='left' if sign == -1 else 'right', va='top')

        new_theta_inc_rad = np.radians(new_theta_inc) + np.pi/2
        ax.arrow(new_theta_inc_rad, 0, 0, 1.2, head_width=0, head_length=0, fc='black', ec='black', alpha=1)

        ax.set_title(f'Diffraction Orders for period = {new_period:.0f} nm, θ_inc = {new_theta_inc:.1f}°\nObjective NA = {na}')
        
        # Update legend
        # legend_elements = [Patch(facecolor=wavelen2rgb(wl), edgecolor='black', label=f'{wl} nm') for wl in wavelengths]
        # ax.legend(handles=legend_elements, loc='upper right', fontsize=8, title="Wavelengths")
        
        fig.canvas.draw_idle()  # Redraw figure

    # Create sliders
    ax_slider_period = plt.axes([0.2, 0.12, 0.65, 0.03])  # Position: [left, bottom, width, height]
    slider_period = Slider(ax_slider_period, 'Period (nm)', 400, 2500, valinit=period, valstep=10)

    ax_slider_theta = plt.axes([0.2, 0.05, 0.65, 0.03])  # Second slider below period slider
    slider_theta = Slider(ax_slider_theta, 'θ_inc (°)', -90, 90, valinit=inc_angle, valstep=1)

    # Connect sliders to update function
    slider_period.on_changed(update)
    slider_theta.on_changed(update)

    update(None)  # Initial plot

    plt.show(block=True)

