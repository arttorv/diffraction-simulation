from matplotlib.widgets import Slider
from colour_utils import *

#----------------------------------
# -------- BASIC FUNCTIONS --------
#----------------------------------

def set_dark_mode(dark_mode=True):
    """Configures Matplotlib for dark mode styling if enabled."""
    if dark_mode:
        plt.rcParams.update({
            'axes.facecolor': 'black',
            'axes.edgecolor': 'grey',
            'figure.facecolor': 'black',
            'figure.edgecolor': 'white',
            'axes.labelcolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'grid.color': 'gray',
            'text.color': 'white',
            'legend.facecolor': 'black',
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
            self.wavelength = (grating_period/m) * (np.sin(np.radians(inc_angle)) + np.sin(np.radians(diff_angle))) 
        if self.diff_angle is None:
            sine_term = - (np.sin(np.radians(inc_angle)) - m * (wavelength / grating_period))
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

def rotatinator_view(period, 
                     inc_angle=0, 
                     angle_res=0.5, 
                     wl_range=[420, 700], 
                     angle_range=[0, 90], 
                     scale_height=1, 
                     m_values=[1, 2, 3, 4],
                     x_grid=False,
                     h_sep_line=False, 
                     dark_mode=True,
                     export_data=False):
    '''
    
    '''
    
    if export_data:
        with open(f'export/diffraction_data_{period}.csv', 'w') as f:
            pass

    angle_vals = np.arange(angle_range[0], angle_range[1] + 1e-6, angle_res)
    diffraction_data = []

    for m in m_values:
        for diff_angle in angle_vals:
            diffraction_line = DiffractionLine(grating_period=period, inc_angle=inc_angle, m=m, diff_angle=diff_angle)
            wl = diffraction_line.wavelength
            if wl_range[0] <= wl <= wl_range[1]:
                color = rgb_from_wavelength(wl) if m != 0 else (0.5, 0.5, 0.5)
                diffraction_data.append((m, diff_angle, wl, color))

    set_dark_mode(dark_mode)
    fig, ax = plt.subplots(figsize=(10, 4 * scale_height))
    plt.subplots_adjust(left=0.1, right=0.9, top=0.7 * scale_height, bottom=0.2 * (1 / scale_height))

    height = 0.5

    for m, diff_angle, wl, color in diffraction_data:
        y_offset = -m
        rect = plt.Rectangle((diff_angle - angle_res / 2, y_offset - height / 2), angle_res, angle_res * 0.5, color=color, alpha=1)
        ax.add_patch(rect)

    mixing_y = -max(m_values) - 1.5
    mixing_height = 0.8

    for diff_angle in angle_vals:
        overlapping = [d for d in diffraction_data if abs(d[1] - diff_angle) < 1e-6]
        if not overlapping:
            continue
        elif len(overlapping) == 1:
            mix_color = overlapping[0][3]
        else:
            wls = [d[2] for d in overlapping]
            mix_color = mix_wavelengths(wls[:3])  # Max 3 wavelengths for safety

        rect = plt.Rectangle((diff_angle - angle_res / 2, mixing_y - mixing_height / 2), angle_res, mixing_height, color=mix_color, alpha=1)
        ax.add_patch(rect)

        if export_data:
            with open(f'export/diffraction_data_{period}.csv', 'a') as f:
                f.write(f"{round(diff_angle,1)},{mix_color[0]},{mix_color[1]},{mix_color[2]}\n")

    plt.rcParams.update({'font.size': 10})
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.xaxis.label.set_size(8)
    ax.yaxis.label.set_size(8)

    ax.set_xlim(angle_range[0], angle_range[1])
    ax.set_ylim(mixing_y - mixing_height, -0.5)
    ax.set_yticks([-m for m in m_values])
    ax.set_yticklabels([f'{m}' for m in m_values])
    ax.set_xlabel("Diffraction Angle (degrees)", fontsize=12)
    ax.set_ylabel("Order", fontsize=12)

    title_line_1 = f"Diffraction Orders for Grating Period = {period:.0f} nm, Wavelength Range: {wl_range[0]} - {wl_range[1]} nm"
    title_line_2 = f"\nIncident Angle = {inc_angle:.1f}°, Angle Resolution = {angle_res}°"
    # ax.set_title(title_line_1 + title_line_2, fontsize=14)

    if x_grid:
        ax.grid(axis='x', linestyle='--', alpha=0.5, color='gray')

    if h_sep_line:
        for m in range(min(m_values) - 1, max(m_values) + 2):
            ax.axhline(-m + 0.5, color='gray', linewidth=2 * scale_height)

    plt.show(block=True)

def polar_orders_overview(period, 
                inc_angle, 
                na, 
                obs_angle,
                m_values=[-4, -3, -2, -1, 0, 1, 2, 3, 4],
                wavelength_start=400, 
                wavelength_stop=700,
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

def diffraction_CIE_view(period, 
                         inc_angle=0, 
                         angle_res=0.5, 
                         wl_range=[420, 700], 
                         angle_range=[0, 90], 
                         annotate_angles=False, 
                         m_values=[1, 2, 3, 4],
                         dark_mode=True):
    """
    Calculates the mixed diffraction orders (mixed_orders) for a range of diffraction angles and plots
    them on the CIE 1931 Chromaticity Diagram.
    
    Parameters
    ----------
    period : float
        Grating period (in the same units as the wavelengths, e.g. nm).
    inc_angle : float, optional
        Incident angle in degrees (default is 0).
    angle_res : float, optional
        Angular resolution for diffraction angles (in degrees; default is 0.5).
    wl_range : list, optional
        Wavelength range [min, max] in nm to consider (default is [420, 700]).
    angle_range : list, optional
        Diffraction angle range [start, stop] in degrees (default is [0, 90]).
    annotate_angles : bool, optional
        If True, annotate each plotted point with its diffraction angle (default is False).
    m_values : list, optional
        List of diffraction orders to consider (default is [1, 2, 3, 4]).
    dark_mode : bool, optional
        If True, configure the plot in dark mode (default is True).
    """
    # Enable dark mode styling if needed.
    set_dark_mode(dark_mode)
    
    # Generate the set of diffraction angles to evaluate.
    angle_vals = np.arange(angle_range[0], angle_range[1] + angle_res, angle_res)
    diffraction_data = []  # Each entry is a tuple: (m, diff_angle, wavelength)
    
    # For each diffraction order and angle, calculate the corresponding wavelength.
    for m in m_values:
        for diff_angle in angle_vals:
            d_line = DiffractionLine(grating_period=period, inc_angle=inc_angle, m=m, diff_angle=diff_angle)
            wl = d_line.wavelength
            if wl_range[0] <= wl <= wl_range[1]:
                diffraction_data.append((m, diff_angle, wl))
    
    # Prepare a list for the mixed orders (one per diffraction angle).
    mixed_points = []  # Each element: (x, y, sRGB_color, diff_angle)

    # Helper: Convert a wavelength to its (x, y, Y) coordinates.
    def wavelength_to_xy(wl):
        XYZ = get_XYZ_from_file(round(wl))
        if XYZ is None:
            return None
        return XYZ_to_xyY(*XYZ)  # (x, y, Y)

    # Helper: Mix two wavelengths in the chromaticity diagram by averaging their xy coordinates.
    def mix_two_wavelengths_xy(wl1, wl2):
        try:
            data = load_data()
        except Exception:
            return None
        wl1 = round(wl1)
        wl2 = round(wl2)
        row1 = data[data["wavelength"] == wl1]
        row2 = data[data["wavelength"] == wl2]
        if row1.empty or row2.empty:
            return None
        X1, Y1, Z1 = row1.iloc[0]["X"], row1.iloc[0]["Y"], row1.iloc[0]["Z"]
        X2, Y2, Z2 = row2.iloc[0]["X"], row2.iloc[0]["Y"], row2.iloc[0]["Z"]
        x1, y1, Y1_xy = XYZ_to_xyY(X1, Y1, Z1)
        x2, y2, Y2_xy = XYZ_to_xyY(X2, Y2, Z2)
        x_mid = (x1 + x2) / 2
        y_mid = (y1 + y2) / 2
        Y_mid = (Y1_xy + Y2_xy) / 2
        return (x_mid, y_mid, Y_mid)
    
    # Loop over each diffraction angle and compute the mixed order.
    for diff_angle in angle_vals:
        # Find all orders at this angle (using a tolerance to account for floating-point comparisons).
        overlapping = [d for d in diffraction_data if abs(d[1] - diff_angle) < 1e-6]
        if not overlapping:
            continue
        if len(overlapping) == 1:
            _, angle_val, wl = overlapping[0]
            xy = wavelength_to_xy(wl)
            sRGB = rgb_from_wavelength(wl)
        else:
            # Mix the first two overlapping orders.
            _, angle_val, wl1 = overlapping[0]
            _, _, wl2 = overlapping[1]
            xy = mix_two_wavelengths_xy(wl1, wl2)
            sRGB = mix_two_wavelengths(wl1, wl2)
        if xy is None or sRGB is None:
            continue
        mixed_points.append((xy[0], xy[1], sRGB, diff_angle))
    
    # Plot the results on the CIE 1931 Chromaticity Diagram.
    import pylab
    import colour.plotting as cl

    # Plot the background diagram (note: using show=False to allow further plotting).
    cl.plot_chromaticity_diagram_CIE1931(show=False)
    
    # Extract xy coordinates and colors.
    xs = [pt[0] for pt in mixed_points]
    ys = [pt[1] for pt in mixed_points]
    colors = [pt[2] for pt in mixed_points]
    
    # Plot each mixed point using triangle markers.
    pylab.scatter(xs, ys, c=colors, marker='^', s=40, edgecolors='black')
    
    # Optionally annotate with the diffraction angle.
    if annotate_angles:
        for pt in mixed_points:
            x, y, col, angle_val = pt
            pylab.annotate(f"{angle_val:.1f}°",
                           xy=(x, y),
                           xytext=(3, 3),
                           textcoords='offset points',
                           fontsize=10,
                           color='black')
    
    # Render the complete plot.
    cl.render(show=True, limits=(-0.1, 0.9, -0.1, 0.9), x_tighten=True, y_tighten=True)
