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

    case = int(input("Enter the simulation number: "))

    # 1 - Diffraction Orders Overview
    if case == 1:
        period = 1900  # Period of the grating
        inc_angle = 0  # Incident angle of whit3e light
        obs_angle = -30
        na = 0.28  # Numerical aperture of objective
        m_values = [-4,-3,-2,-1,0,1,2,3,4,5]
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
        inc_angle = 0 # Incident angle in degree
        obs_angle = -17.5
        na = 0.28 # Numerical aperture of objective 
        m_values = [0, 1, 2, 3, 4]
        wavelength_step = 4
        dark_mode = False
        objective_view(period, 
                       inc_angle, 
                       obs_angle, 
                       na=na, 
                       wl_step=wavelength_step, 
                       m_values=m_values)

    # 3 - Rotatinator View Plot
    elif case == 3:
        period = 1900  # Period of the grating
        inc_angle = 0 # Incident angle in degree
        m_values = [1, 2, 3, 4]
        angle_resolution = 1
        dark_mode = True
        rotatinator_view(period=period, 
                         inc_angle=inc_angle, 
                         angle_res=angle_resolution, 
                         m_values=m_values,
                         dark_mode=dark_mode)

if __name__ == "__main__":
    main()