import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

'''
Using datasets from CIE 1931 available at https://cie.co.at/data-tables 
'''

#----------------------------------
# -------- BASIC UTILITIES --------
#----------------------------------

def load_data(file="CIE_cc_1931_2deg.csv"):
    """
    Load the CSV file that contains the XYZ values.
    The file is assumed to have no header and four columns:
    wavelength (nm), X, Y, and Z.
    
    Returns:
        DataFrame: A pandas DataFrame with columns ['wavelength', 'X', 'Y', 'Z'].
    """
    try:
        data = pd.read_csv(file, header=None, names=["wavelength", "X", "Y", "Z"])
        return data
    except FileNotFoundError:
        print(f"Error: '{file}' not found in the current folder.")
        return None

def get_XYZ_from_file(wavelength, data=None):
    """
    Extract the XYZ values for a given wavelength from the provided data.
    
    Args:
        wavelength (int or float): The wavelength (nm) to look up.
        data (DataFrame, optional): The pre-loaded data. If None, load the data.
    
    Returns:
        tuple: (X, Y, Z) values if found, or None otherwise.
    """
    if data is None:
        data = load_data()
        if data is None:
            return None
    
    wavelength = round(wavelength)
    row = data[data["wavelength"] == wavelength]
    if row.empty:
        print(f"Wavelength {wavelength} nm not found in the data.")
        return None

    X_val = row.iloc[0]["X"]
    Y_val = row.iloc[0]["Y"]
    Z_val = row.iloc[0]["Z"]
    return (X_val, Y_val, Z_val)

def XYZ_to_xyY(X, Y, Z):
    """
    Convert from XYZ to xyY.
    """
    total = X + Y + Z
    if total == 0:
        return 0, 0, Y
    x = X / total
    y = Y / total
    return x, y, Y

def xyY_to_XYZ(x, y, Y):
    """
    Convert from xyY to XYZ.
    """
    if y == 0:
        return 0, Y, 0
    X = (x / y) * Y
    Z = ((1 - x - y) / y) * Y
    return X, Y, Z

def XYZ_to_sRGB(X, Y, Z):
    """
    Convert XYZ to sRGB using the standard transformation for D65 white point.
    https://en.wikipedia.org/wiki/SRGB 
    """
    # Transformation matrix from XYZ to linear sRGB.
    M = np.array([[ 3.2406, -1.5372, -0.4986],
                  [-0.9689,  1.8758,  0.0415],
                  [ 0.0557, -0.2040,  1.0570]])
    # Convert XYZ to linear sRGB.
    RGB_linear = np.dot(M, np.array([X, Y, Z]))
    # Clip any negative values.
    RGB_linear = np.clip(RGB_linear, 0, None)

    # Gamma correction function.
    def gamma_correct(c):
        return 12.92 * c if c <= 0.0031308 else 1.055 * (c ** (1 / 2.4)) - 0.055

    rgb = np.array([gamma_correct(c) for c in RGB_linear])
    # Clip to the [0, 1] range.
    return np.clip(rgb, 0, 1)

def mix_two_wavelengths(wl1, wl2):
    """
    Given two wavelengths, this function:
      1. Loads the XYZ values from the CSV file.
      2. Finds the rows corresponding to the given wavelengths.
      3. Converts each to xyY.
      4. Computes the midpoint in the chromaticity diagram (averaging x, y and Y).
      5. Converts the averaged xyY back to XYZ.
      6. Converts the XYZ to sRGB.
      7. Prints both the XYZ and xyY values.
    """
    # Load CSV file without headers, assigning columns: wavelength, X, Y, Z.
    try:
        data = pd.read_csv("CIE_cc_1931_2deg.csv", header=None, names=["wavelength", "X", "Y", "Z"])
    except FileNotFoundError:
        print("Error: 'CIE_cc_1931_2deg.csv' not found in the current folder.")
        return None

    wl1 = round(wl1)
    wl2 = round(wl2)
    # Retrieve the rows corresponding to each wavelength.
    row1 = data[data["wavelength"] == wl1]
    row2 = data[data["wavelength"] == wl2]

    if row1.empty or row2.empty:
        print("Error: One of the wavelengths was not found in the data.")
        return None

    # Extract the XYZ values.
    X1, Y1, Z1 = row1.iloc[0]["X"], row1.iloc[0]["Y"], row1.iloc[0]["Z"]
    X2, Y2, Z2 = row2.iloc[0]["X"], row2.iloc[0]["Y"], row2.iloc[0]["Z"]

    # Convert each to xyY.
    x1, y1, Y1_xyY = XYZ_to_xyY(X1, Y1, Z1)
    x2, y2, Y2_xyY = XYZ_to_xyY(X2, Y2, Z2)

    # For equal intensities, take the midpoint in the chromaticity diagram.
    x_mid = (x1 + x2) / 2
    y_mid = (y1 + y2) / 2
    Y_mid = (Y1_xyY + Y2_xyY) / 2

    # Convert the averaged xyY back to XYZ.
    X_mid, Y_mid, Z_mid = xyY_to_XYZ(x_mid, y_mid, Y_mid)

    # Convert the mixed XYZ to sRGB.
    rgb = XYZ_to_sRGB(X_mid, Y_mid, Z_mid)
    return rgb

def rgb_from_wavelength(wavelength, data=None):
    """
    Retrieve the XYZ values for a given wavelength and convert them to sRGB.
    
    Args:
        wavelength (int or float): The wavelength (nm) to look up.
        data (DataFrame, optional): The pre-loaded data. If None, load the data.
    
    Returns:
        np.array or None: The resulting sRGB color as a numpy array, or None if not found.
    """
    XYZ = get_XYZ_from_file(wavelength, data)
    if XYZ is None:
        return None

    rgb = XYZ_to_sRGB(*XYZ)
    
    return rgb



#----------------------------------
# -------- PLOT FUNCTIONS --------
#----------------------------------

def plot_color(rgb, title="Mixed Color"):
    """
    Plot a rectangle filled with the given sRGB color.
    """
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, color=rgb))
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.axis('off')
    plt.title(title)
    plt.show()

def display_color_from_wavelength(wavelength):
    """
    Given a wavelength, locate the corresponding XYZ values from the CSV file,
    convert these to sRGB, print the XYZ values and display a rectangle filled
    with the resulting color.
    """
    # Load CSV file (no header) and assign column names.
    try:
        data = pd.read_csv("CIE_cc_1931_2deg.csv", header=None, names=["wavelength", "X", "Y", "Z"])
    except FileNotFoundError:
        print("Error: 'CIE_cc_1931_2deg.csv' not found in the current folder.")
        return

    # Find the row corresponding to the given wavelength.
    row = data[data["wavelength"] == wavelength]
    if row.empty:
        print(f"Wavelength {wavelength} nm not found in the data.")
        return

    # Extract the XYZ values.
    X_val = row.iloc[0]["X"]
    Y_val = row.iloc[0]["Y"]
    Z_val = row.iloc[0]["Z"]

    print("XYZ values:", (X_val, Y_val, Z_val))

    # Convert the XYZ values to sRGB.
    rgb = XYZ_to_sRGB(X_val, Y_val, Z_val)
    print("sRGB values:", rgb)

    # Plot a rectangle with the resulting sRGB color.
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, color=rgb))
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.axis('off')
    plt.title(f"Wavelength: {wavelength} nm")
    plt.show()

# Example usage:
if __name__ == "__main__":
    # Replace with your desired wavelengths, e.g., 450 nm and 700 nm.
    color_a = 600
    color_b = 423
    mixed_rgb = mix_two_wavelengths(color_a, color_b)
    mono_rgb = rgb_from_wavelength(580)
    if mixed_rgb is not None:
        plot_color(mono_rgb, title=f"Mixed Color from {color_a} nm & {color_b} nm")