%matplotlib qt 

import numpy as np
import matplotlib.pyplot as plt
from mp_api.client import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.diffraction.tem import TEMCalculator
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from matplotlib.offsetbox import AnnotationBbox, TextArea

def get_interactive_diffraction_pattern(material_id, zone_axis, api_key=None, voltage=200, max_range=20):
    """
    Retrieves a structure from Materials Project and generates an interactive TEM diffraction pattern
    for a specific zone axis, limited to a maximum range from the central beam.
    Allows mouse hover to display intensity and d-spacing information.

    Args:
        material_id (str): The Materials Project ID (e.g., 'mp-149')
        zone_axis (list): The zone axis as a list of three integers [h, k, l]
        api_key (str, optional): Your Materials Project API key. If None, uses environment variable.
        voltage (float, optional): Electron beam energy in keV. Default is 200 keV.
        max_range (float, optional): Maximum range from center in nm⁻¹. Default is 20 nm⁻¹.

    Returns:
        tuple: (fig, ax, pattern) matplotlib figure, axis objects, and diffraction pattern data
    """
    # Convert zone axis to tuple of integers
    zone_axis = tuple(map(int, zone_axis))

    # Get the structure from Materials Project
    with MPRester(api_key) as mpr:
        structure = mpr.get_structure_by_material_id(material_id)
        print(f"Retrieved structure: {structure.formula}")

    # Create TEM calculator
    tem_calculator = TEMCalculator(voltage=voltage, beam_direction=zone_axis)

    # Generate the diffraction pattern
    pattern = tem_calculator.get_pattern(structure)

    # Create figure for the diffraction pattern
    fig, ax = plt.subplots(figsize=(10, 10))

    # Extract pattern data from the DataFrame
    pattern_data = pattern.to_dict('records')

    # Filter based on the actual distance
    filtered_pattern_data = []
    for dot in pattern_data:
        x, y = dot['Position']
        distance_from_center = np.sqrt(x**2 + y**2)
        if distance_from_center <= max_range:
            filtered_pattern_data.append(dot)
    
    print(f"Total spots in pattern: {len(pattern_data)}")
    print(f"Spots within {max_range} nm⁻¹ range: {len(filtered_pattern_data)}")
    
    pattern_data = filtered_pattern_data
    
    # Scale factor for spot sizes based on intensity
    intensities = [dot['Intensity (norm)'] for dot in pattern_data]
    max_intensity = max(intensities) if intensities else 1
    intensity_scale = 200 / max_intensity

    # Dictionary to store the spots with their data (for hover functionality)
    spots_dict = {}

    # Plot the diffraction spots
    for dot in pattern_data:
        # Skip the direct beam for now (it will be added separately)
        hkl = dot['(hkl)']
        if hkl == '(0, 0, 0)':
            continue

        # Extract position and data
        x, y = dot['Position']
        intensity = dot['Intensity (norm)']
        
        # Calculate d-spacing from the reciprocal lattice vector magnitude
        # d = 1/|g| where |g| is the magnitude of the reciprocal lattice vector
        g_magnitude = np.sqrt(x**2 + y**2)  # nm^-1
        d_spacing = 1.0/g_magnitude if g_magnitude > 0 else float('inf')  # in nm
        
        # Convert d-spacing from nm to Ångstroms (1 nm = 10 Å)
        d_spacing_angstrom = d_spacing * 10.0
        
        # Plot the spot
        spot = ax.scatter(x, y, s=intensity*intensity_scale, c='black', picker=5)
        
        # Store the spot data for hover functionality
        spots_dict[spot] = {
            'hkl': hkl,
            'intensity': intensity,
            'd_spacing': d_spacing_angstrom,
            'position': (x, y)
        }

        # Label the spots with their Miller indices
        hkl_str = ''.join(c for c in str(hkl) if c.isdigit() or c == '-')
        ax.text(x, y+0.05*max_range/20, hkl_str, ha='center', fontsize=8)

    # Add a central spot (direct beam)
    central_spot = ax.scatter(0, 0, s=200, c='black', marker='o', picker=5)
    ax.text(0, 0.02*max_range/20, '000', ha='center', fontsize=8)
    spots_dict[central_spot] = {
        'hkl': '(0, 0, 0)',
        'intensity': 1.0,  # Normalized intensity for central beam
        'd_spacing': float('inf'),  # Infinite d-spacing for central beam
        'position': (0, 0)
    }

    # Set plot limits to the specified max_range
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)

    # Add labels and title
    ax.set_xlabel('Distance (nm⁻¹)')
    ax.set_ylabel('Distance (nm⁻¹)')
    ax.set_title(f'TEM Diffraction Pattern for {structure.formula}\nZone Axis: [{zone_axis[0]} {zone_axis[1]} {zone_axis[2]}] (Limited to {max_range} nm⁻¹)')
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.7)

    # Create an annotation box for hover information (initially empty)
    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20),
                      textcoords="offset points",
                      bbox=dict(boxstyle="round,pad=0.5", fc="yellow", alpha=0.8),
                      arrowprops=dict(arrowstyle="->"),
                      fontsize=9,
                      visible=False)

    # Function to handle hover events
    def hover(event):
        if event.inaxes == ax:
            # Check if we're near any of the spots
            for spot, data in spots_dict.items():
                contains, _ = spot.contains(event)
                if contains:
                    x, y = data['position']
                    hkl = data['hkl']
                    intensity = data['intensity']
                    d_spacing = data['d_spacing']
                    
                    # Update annotation with spot information
                    annot.xy = (x, y)
                    text = f"hkl: {hkl}\nIntensity: {intensity:.4f}\nd-spacing: {d_spacing:.4f} Å"
                    annot.set_text(text)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                    return
            
            # If we're not hovering over any spot, hide the annotation
            if annot.get_visible():
                annot.set_visible(False)
                fig.canvas.draw_idle()

    # Connect the hover event handler
    fig.canvas.mpl_connect("motion_notify_event", hover)

    return fig, ax, pattern

def get_interactive_diffraction_pattern_from_cif(cif_file, zone_axis, voltage=200, max_range=20):
    """
    Generates an interactive TEM diffraction pattern from a local CIF file for a specific zone axis,
    limited to a maximum range from the central beam.
    Allows mouse hover to display intensity and d-spacing information.

    Args:
        cif_file (str): Path to the CIF file
        zone_axis (list): The zone axis as a list of three integers [h, k, l]
        voltage (float, optional): Electron beam energy in keV. Default is 200 keV.
        max_range (float, optional): Maximum range from center in nm⁻¹. Default is 20 nm⁻¹.

    Returns:
        tuple: (fig, ax, pattern) matplotlib figure, axis objects, and diffraction pattern data
    """
    # Convert zone axis to tuple of integers
    zone_axis = tuple(map(int, zone_axis))

    # Load structure from CIF file
    structure = Structure.from_file(cif_file)
    print(f"Loaded structure from file: {structure.formula}")

    # Create TEM calculator
    tem_calculator = TEMCalculator(voltage=voltage, beam_direction=zone_axis)

    # Generate the diffraction pattern
    pattern = tem_calculator.get_pattern(structure)

    # Create figure for the diffraction pattern
    fig, ax = plt.subplots(figsize=(10, 10))

    # Extract pattern data
    pattern_data = pattern.to_dict('records')
    
    # Filter based on the actual distance
    filtered_pattern_data = []
    for dot in pattern_data:
        x, y = dot['Position']
        distance_from_center = np.sqrt(x**2 + y**2)
        if distance_from_center <= max_range:
            filtered_pattern_data.append(dot)
    
    print(f"Total spots in pattern: {len(pattern_data)}")
    print(f"Spots within {max_range} nm⁻¹ range: {len(filtered_pattern_data)}")
    
    pattern_data = filtered_pattern_data

    # Scale factor for spot sizes based on intensity
    intensities = [dot['Intensity (norm)'] for dot in pattern_data]
    max_intensity = max(intensities) if intensities else 1
    intensity_scale = 200 / max_intensity

    # Dictionary to store the spots with their data (for hover functionality)
    spots_dict = {}

    # Plot the diffraction spots
    for dot in pattern_data:
        # Skip the direct beam for now (it will be added separately)
        hkl = dot['(hkl)']
        if hkl == '(0, 0, 0)':
            continue

        # Extract position and data
        x, y = dot['Position']
        intensity = dot['Intensity (norm)']
        
        # Calculate d-spacing from the reciprocal lattice vector magnitude
        # d = 1/|g| where |g| is the magnitude of the reciprocal lattice vector
        g_magnitude = np.sqrt(x**2 + y**2)  # nm^-1
        d_spacing = 1.0/g_magnitude if g_magnitude > 0 else float('inf')  # in nm
        
        # Convert d-spacing from nm to Ångstroms (1 nm = 10 Å)
        d_spacing_angstrom = d_spacing * 10.0
        
        # Plot the spot
        spot = ax.scatter(x, y, s=intensity*intensity_scale, c='black', picker=5)
        
        # Store the spot data for hover functionality
        spots_dict[spot] = {
            'hkl': hkl,
            'intensity': intensity,
            'd_spacing': d_spacing_angstrom,
            'position': (x, y)
        }

        # Label the spots with their Miller indices
        hkl_str = ''.join(c for c in str(hkl) if c.isdigit() or c == '-')
        ax.text(x, y+0.05*max_range/20, hkl_str, ha='center', fontsize=8)

    # Add a central spot (direct beam)
    central_spot = ax.scatter(0, 0, s=200, c='black', marker='o', picker=5)
    ax.text(0, 0.02*max_range/20, '000', ha='center', fontsize=8)
    spots_dict[central_spot] = {
        'hkl': '(0, 0, 0)',
        'intensity': 1.0,  # Normalized intensity for central beam
        'd_spacing': float('inf'),  # Infinite d-spacing for central beam
        'position': (0, 0)
    }

    # Set plot limits to the specified max_range
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)

    # Add labels and title
    ax.set_xlabel('Distance (nm⁻¹)')
    ax.set_ylabel('Distance (nm⁻¹)')
    cif_name = cif_file.split('/')[-1].split('.')[0]
    ax.set_title(f'TEM Diffraction Pattern for {cif_name}\nZone Axis: [{zone_axis[0]} {zone_axis[1]} {zone_axis[2]}] (Limited to {max_range} nm⁻¹)')
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.7)

    # Create an annotation box for hover information (initially empty)
    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20),
                      textcoords="offset points",
                      bbox=dict(boxstyle="round,pad=0.5", fc="yellow", alpha=0.8),
                      arrowprops=dict(arrowstyle="->"),
                      fontsize=9,
                      visible=False)

    # Function to handle hover events
    def hover(event):
        if event.inaxes == ax:
            # Check if we're near any of the spots
            for spot, data in spots_dict.items():
                contains, _ = spot.contains(event)
                if contains:
                    x, y = data['position']
                    hkl = data['hkl']
                    intensity = data['intensity']
                    d_spacing = data['d_spacing']
                    
                    # Update annotation with spot information
                    annot.xy = (x, y)
                    text = f"hkl: {hkl}\nIntensity: {intensity:.4f}\nd-spacing: {d_spacing:.4f} Å"
                    annot.set_text(text)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                    return
            
            # If we're not hovering over any spot, hide the annotation
            if annot.get_visible():
                annot.set_visible(False)
                fig.canvas.draw_idle()

    # Connect the hover event handler
    fig.canvas.mpl_connect("motion_notify_event", hover)

    return fig, ax, pattern

# Example usage
if __name__ == "__main__":
    # Example: Get interactive diffraction pattern directly from Materials Project
    material_id = "mp-30"  # Silicon
    zone_axis = [1, 0, 1]   # [101] zone axis
    api_key = "daOUQsZxLXFDwpCXnB0uBMoXiicXZ8nq"  # Replace with your Materials Project API key
    
    # Generate interactive diffraction pattern
    fig, ax, _ = get_interactive_diffraction_pattern(material_id, zone_axis, api_key, max_range=10)
    plt.show()
    
    # For CIF files, use:
    # fig, ax, _ = get_interactive_diffraction_pattern_from_cif("path/to/your/file.cif", zone_axis, max_range=10)
    # plt.show()

    #daOUQsZxLXFDwpCXnB0uBMoXiicXZ8nq