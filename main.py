from diffraction_utils import *

def main():
    """
    Main function to run different functionalities related to diffraction grating analysis.
    The function presents a menu of options to the user and executes the corresponding functionality based on the user's choice.
    The available options are:
    1 - Diffraction Orders Overview:
        # Generates a polar plot showing the diffraction orders within a specified numerical aperture.
    2 - Objective View Plot:
        # Simulates and plots the view through an objective lens observing at an angle, showing the captured diffraction orders.
    3 - Rotatinator View Plot:
        # Simulates and plots the diffraction orders as seen by rotating the observing angle, showing the captured color at different angles.
    """

    print("Choose a function to run:")
    print("1 - Diffraction Orders Overview")
    print("2 - Objective View Plot")
    print("3 - Rotatinator View Plot")
    print("4 - Diffraction CIE View Plot")
    print("0 - Quit")

    case = int(input("Enter the simulation number: "))

    # 1 - Diffraction Orders Overview
    if case == 1:
        period = 1700  # Period of the grating
        inc_angle = 0  # Incident angle of white light
        obs_angle = -0 # Observation angle
        na = 0.28  # Numerical aperture of objective
        m_values = [-4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
        wavelength_interval = 10
        dark_mode = False 
        polar_orders_overview(period, 
                                inc_angle, 
                                na, 
                                obs_angle, 
                                m_values=m_values, 
                                wavelength_interval=wavelength_interval, 
                                dark_mode=dark_mode)

    # 2 - Objective View Plot
    elif case == 2:
        period = 1900  # Period of the grating
        inc_angle = 0  # Incident angle in degree
        obs_angle = -17.5 # Observation angle
        na = 0.28  # Numerical aperture of objective 
        m_values = [0, 1, 2, 3, 4]
        wavelength_step = 4
        dark_mode = True
        objective_view(period, 
                        inc_angle, 
                        obs_angle, 
                        na=na, 
                        wl_step=wavelength_step, 
                        m_values=m_values, 
                        dark_mode=dark_mode)

    # 3 - Rotatinator View Plot
    elif case == 3:
        period = 2365 # Period of the grating
        inc_angle = 0  # Incident angle in degree
        m_values = [1, 2, 3, 4, 5] # Orders to plot
        angle_resolution = 1 # Colors are plotted at this resolution
        angle_range = [4, 84] # [start angle, stop angle]
        wl_range=[420, 680] # Wavelength range in nm
        scale_height = 1.2
        x_grid = False
        h_sep_line = False
        dark_mode = True
        export_data = True

        rotatinator_view(period=period, 
                            inc_angle=inc_angle, 
                            angle_res=angle_resolution,
                            angle_range=angle_range, 
                            wl_range=wl_range,
                            m_values=m_values,
                            scale_height=scale_height,
                            x_grid=x_grid,
                            h_sep_line=h_sep_line,
                            dark_mode=dark_mode,
                            export_data=export_data)
    
    # 4 - Diffraction CIE View Plot
    elif case == 4:
        period = 1900       # Grating period (nm)
        inc_angle = 0       # Incident angle in degrees
        angle_resolution = 1  # Angular resolution in degrees
        wl_range = [420, 700]   # Wavelength range in nm
        angle_range = [0, 90]   # Diffraction angle range in degrees
        annotate_angles = True  # Annotate each point with its diffraction angle
        m_values = [1, 2, 3, 4] # Diffraction orders to consider
        dark_mode = False

        diffraction_CIE_view(period=period, 
                             inc_angle=inc_angle, 
                             angle_res=angle_resolution, 
                             wl_range=wl_range, 
                             angle_range=angle_range, 
                             annotate_angles=annotate_angles, 
                             m_values=m_values, 
                             dark_mode=dark_mode)
        

    else:
        print("Invalid selection. Please choose again.")

if __name__ == "__main__":
    main()
