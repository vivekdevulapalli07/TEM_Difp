import numpy as np
import matplotlib.pyplot as plt
from mp_api.client import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.diffraction.tem import TEMCalculator
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

def get_diffraction_pattern(material_id, zone_axis, api_key=None, voltage=200, max_range=20):
    """
    Retrieves a structure from Materials Project and generates a TEM diffraction pattern
    for a specific zone axis, limited to a maximum range from the central beam.
    The pattern is correctly limited to the specified max_range.

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

    # CRITICALLY IMPORTANT:
    # We DO NOT scale the positions. Instead, we strictly filter based on the actual distance.
    # This ensures that only spots within the true max_range are shown.
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

    # Plot the diffraction spots
    for dot in pattern_data:
        # Skip the direct beam for now (it will be added separately)
        hkl = dot['(hkl)']
        if hkl == '(0, 0, 0)':
            continue

        # Extract position
        x, y = dot['Position']
        intensity = dot['Intensity (norm)']

        # Plot the spot
        ax.scatter(x, y, s=intensity*intensity_scale, c='black')

        # Label the spots with their Miller indices
        # Convert to string first, then extract digits and minus signs safely
        hkl_str = ''.join(c for c in str(hkl) if c.isdigit() or c == '-')
        ax.text(x, y+0.05*max_range/20, hkl_str, ha='center', fontsize=8)

    # Add a central spot (direct beam)
    ax.scatter(0, 0, s=200, c='black', marker='o')
    ax.text(0, 0.02*max_range/20, '000', ha='center', fontsize=8)

    # Set plot limits to the specified max_range
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)

    # Add labels and title
    ax.set_xlabel('Distance (nm⁻¹)')
    ax.set_ylabel('Distance (nm⁻¹)')
    ax.set_title(f'TEM Diffraction Pattern for {structure.formula}\nZone Axis: [{zone_axis[0]} {zone_axis[1]} {zone_axis[2]}] (Limited to {max_range} nm⁻¹)')
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.7)

    return fig, ax, pattern

def get_diffraction_pattern_from_cif(cif_file, zone_axis, voltage=200, max_range=20):
    """
    Generates a TEM diffraction pattern from a local CIF file for a specific zone axis,
    limited to a maximum range from the central beam.
    The pattern is correctly limited to the specified max_range.

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
    
    # CRITICALLY IMPORTANT:
    # We DO NOT scale the positions. Instead, we strictly filter based on the actual distance.
    # This ensures that only spots within the true max_range are shown.
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

    # Plot the diffraction spots
    for dot in pattern_data:
        # Skip the direct beam for now (it will be added separately)
        hkl = dot['(hkl)']
        if hkl == '(0, 0, 0)':
            continue

        # Extract position
        x, y = dot['Position']
        intensity = dot['Intensity (norm)']

        # Plot the spot
        ax.scatter(x, y, s=intensity*intensity_scale, c='black')

        # Label the spots with their Miller indices
        # Convert to string first, then extract digits and minus signs safely
        hkl_str = ''.join(c for c in str(hkl) if c.isdigit() or c == '-')
        ax.text(x, y+0.05*max_range/20, hkl_str, ha='center', fontsize=8)

    # Add a central spot (direct beam)
    ax.scatter(0, 0, s=200, c='black', marker='o')
    ax.text(0, 0.02*max_range/20, '000', ha='center', fontsize=8)

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

    return fig, ax, pattern

# Add a utility function to check the number of spots at different ranges
def analyze_diffraction_spots(material_id, zone_axis, api_key=None, voltage=200, ranges=[5, 10, 20, 50]):
    """
    Analyzes how many diffraction spots are visible at different max_range values.
    This helps understand the right range to use for a given material.
    
    Args:
        material_id (str): The Materials Project ID
        zone_axis (list): Zone axis as a list of three integers
        api_key (str, optional): Materials Project API key
        voltage (float, optional): Electron beam energy in keV
        ranges (list, optional): List of ranges to analyze
        
    Returns:
        dict: Dictionary mapping ranges to spot counts
    """
    zone_axis = tuple(map(int, zone_axis))
    
    with MPRester(api_key) as mpr:
        structure = mpr.get_structure_by_material_id(material_id)
        print(f"Analyzing diffraction pattern for {structure.formula}, zone axis {zone_axis}")
    
    # Create TEM calculator
    tem_calculator = TEMCalculator(voltage=voltage, beam_direction=zone_axis)
    
    # Generate the diffraction pattern
    pattern = tem_calculator.get_pattern(structure)
    pattern_data = pattern.to_dict('records')
    
    # Count spots at each range
    results = {}
    for max_range in ranges:
        count = 0
        for dot in pattern_data:
            x, y = dot['Position']
            distance = np.sqrt(x**2 + y**2)
            if distance <= max_range:
                count += 1
        results[max_range] = count
        print(f"Range {max_range} nm⁻¹: {count} spots")
    
    return results

# Example usage
if __name__ == "__main__":
    # Example: Get diffraction pattern directly from Materials Project
    material_id = "mp-30"  # Silicon
    zone_axis = [1, 0, 1]   # [001] zone axis
    api_key = "daOUQsZxLXFDwpCXnB0uBMoXiicXZ8nq"  # Replace with your Materials Project API key
    
    # First analyze how many spots are at different ranges to help choose appropriate values
    analyze_diffraction_spots(material_id, zone_axis, api_key, ranges=[5, 10, 15, 20, 30, 50])
    
    # Generate diffraction patterns at different ranges
    fig10, ax10, _ = get_diffraction_pattern(material_id, zone_axis, api_key, max_range=10)
    #plt.savefig(f"{material_id}_zone_{zone_axis[0]}{zone_axis[1]}{zone_axis[2]}_10nm-1.png", dpi=300)

    #daOUQsZxLXFDwpCXnB0uBMoXiicXZ8nq