from diffraction_utils import *

def main():
    """
    Main function to run different functionalities related to diffraction grating analysis.
    The available options are:
    1 - Diffraction Orders Overview:
            Generates a polar plot showing the diffraction orders within a specified numerical aperture.
    2 - Rotatinator View Plot:
            Calculates and plots the colors of the diffraction orders as seen by rotating the 
            observing angle onto the grating
    3 - CIE View: 
            Calculates the expected colors for a range of diffraction angles and plots the color
            on the CIE diagram
    """

    print("Choose a function to run:")
    print("1 - Diffraction Orders Overview")
    print("2 - Rotatinator View Plot")
    print("3 - Diffraction CIE View Plot")
    print("0 - Quit")

    while True:

        try:
            case = int(input("Enter the function number: "))
        except ValueError:
            print("Invalid input. Please enter a number.")
            break

        # 1 - Diffraction Orders Overview
        if case == 1:
            period = 900  # Period of the grating
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
            break

        # 2 - Rotatinator View Plot
        elif case == 2:
            period = 2000 # Period of the grating
            inc_angle = 0  # Incident angle in degree
            m_values = [1, 2, 3, 4] # Orders to plot
            angle_resolution = 1 # Colors are plotted at this resolution
            angle_range = [0, 90] # [start angle, stop angle]
            wl_range=[420, 680] # Wavelength range in nm
            scale_height = 1.2
            x_grid = False
            h_sep_line = False
            dark_mode = False
            export_data = False # Exported data is found in "export"-folder

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
            break

        # 3 - Diffraction CIE View Plot
        elif case == 3:
            period = 1200       # Grating period (nm)
            inc_angle = 0       # Incident angle in degrees
            angle_resolution = 1  # Angular resolution in degrees
            wl_range = [360, 740]   # Wavelength range in nm
            angle_range = [0, 80]   # Diffraction angle range in degrees
            annotate_angles = True  # Annotate each point with its diffraction angle
            m_values = [1, 2] # Diffraction orders to consider
            dark_mode = False

            diffraction_CIE_view(period=period, 
                                inc_angle=inc_angle, 
                                angle_res=angle_resolution, 
                                wl_range=wl_range, 
                                angle_range=angle_range, 
                                annotate_angles=annotate_angles, 
                                m_values=m_values, 
                                dark_mode=dark_mode)
            break 

        elif case == 0:
            print("Quiting.")
            break

        else:
            print("Invalid selection. Please choose again.")
            

if __name__ == "__main__":
    main()
