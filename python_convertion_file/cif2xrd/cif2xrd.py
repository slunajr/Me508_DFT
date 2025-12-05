# General python modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Materials science python modules
!pip install --upgrade pymatgen
import pymatgen
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Change directory to folder where you would like to create a folder to store the generated figures
target_directory = # Insert path to target folder
os.chdir(target_directory) # Change current directory to target folder

# Get cif file
cif_file = # Insert path to cif file for structure of interest

# Create folder to save figures
figure_folder = 'Figures' # Create name for folder to hold generated figures
os.makedirs(figure_folder, exist_ok=True) # If folder does not exist, create folder
full_path = os.path.join(os.getcwd(), figure_folder) # Get path to created folder

# Get final structure from cif file
structure = Structure.from_file(cif_file)

# Get string for compound
compound = # Insert structure compound name here

# Initialize symmetry analyzer
sym_analyzer = SpacegroupAnalyzer(structure)

# Get symmetry and lattice info for structure
space_group_symbol = sym_analyzer.get_space_group_symbol()
space_group_number = sym_analyzer.get_space_group_number()
point_group_symbol = sym_analyzer.get_point_group_symbol()
crystal_system = sym_analyzer.get_crystal_system()
conventional_structure = sym_analyzer.get_conventional_standard_structure()

# Print the properties
print("-" * 30)
print(f"Symmetry Properties for the Structure:")
print(f"  Space Group Symbol (International): {space_group_symbol}")
print(f"  Space Group Number (International): {space_group_number}")
print(f"  Point Group Symbol: {point_group_symbol}")
print(f"  Crystal System: {crystal_system}")
print("-" * 30)
print(f"\nConventional standard structure formula: {conventional_structure.formula}")
print(f"Conventional lattice parameters:")
print(f"  a, b, c: {conventional_structure.lattice.abc}")
print(f"  alpha, beta, gamma: {conventional_structure.lattice.angles}")

# Parameters
wl = 1.5418; # Wavelength for Cu K-alpha radiation, replace with desired radiation wavelength if necessary

# Initialize the XRD calculator (default wavelength is Cu K-alpha, lambda=1.5418 A)
xrd_calculator = XRDCalculator(wavelength=wl)

# Calculate the diffraction pattern
pattern = xrd_calculator.get_pattern(structure, scaled=True, two_theta_range=(0, 90))
two_thetas = pattern.x  # Get 2theta values for peaks
intensities = pattern.y  # Get intensity values for peaks
d_hkls = pattern.d_hkls # Get list of hkl information

# Get peak hkls
hkls = pattern.hkls
hkl_indices_only = []

# Print peak information
for i, two_theta in enumerate(two_thetas):
    intensity = intensities[i]
    hkl_list = hkls[i]
    d_hkl = d_hkls[i]

    print(f"\nPeak at 2θ = {two_theta:.2f}° (Intensity: {intensity:.2f}, d-spacing: {d_hkl:.3f} Å):")
    for hkl_info in hkl_list:
        hkl_tuple = hkl_info['hkl']
        hkl = "".join(map(str,hkl_tuple))
        multiplicity = hkl_info['multiplicity']
        print(f"  - HKL: {hkl}, Multiplicity: {multiplicity}")

# Plot peaks only
peaks_plot = xrd_calculator.get_plot(structure)
peaks_plot.set_title(f"Post-DFT XRD Pattern for {compound}", fontsize=36)

# Save ideal xrd plot to Figures folder
ideal_filename = f"post_dft_ideal_xrd_{compound}.png"
plt.savefig(f"{full_path}/{ideal_filename}")

# Show plot
plt.show()

# Create broadening function
def apply_broadening(two_theta_ideal, intensities_ideal, fwhm, two_theta_range=(0, 90), step=0.01):
    """
    Applies Gaussian broadening to an ideal XRD pattern.

    Args:
        two_theta_ideal: Ideal 2-theta peak positions.
        intensities_ideal: Ideal peak intensities.
        fwhm: Full Width at Half Maximum for all peaks (constant broadening for simplicity).
        two_theta_range: The range for the output pattern.
        step: The step size for the output 2-theta array.
    """
    two_theta_broad = np.arange(two_theta_range[0], two_theta_range[1], step)
    broadened_intensities = np.zeros(two_theta_broad.shape)

    # Gaussian broadening formula
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2))) # Convert FWHM to sigma

    for two_t, intensity in zip(two_theta_ideal, intensities_ideal):
        # Add contribution of each peak to the overall pattern
        broadened_intensities += intensity * np.exp(-(two_theta_broad - two_t)**2 / (2 * sigma**2))

    return two_theta_broad, broadened_intensities

# Apply broadening
fwhm = 0.2
two_theta_broad, intensities_broad = apply_broadening(two_thetas, intensities, fwhm)

# Plot peaks with broadening with hkl labels
plt.figure(figsize=(10, 6))
plt.plot(two_theta_broad, intensities_broad, label=f'Broadened Peaks (FWHM={fwhm}°)', color='b')
plt.stem(two_thetas, intensities, linefmt='r-', markerfmt=' ', basefmt=' ', label='Ideal Peaks') # Ideal peaks as vertical lines

plt.title(f"XRD Pattern with Peak Broadening ({compound}, Post-DFT)")
plt.xlabel('2$\\theta$ (degrees)')
plt.ylabel('Intensity')
plt.legend()
plt.xlim(0, 90)
plt.ylim(0, (max(intensities)+20))
plt.grid(True)

# Add peak hkl labels
for i in range(len(two_thetas)):
    hkl_list = hkls[i]
    hkl_labels = []
    for hkl_info in hkl_list:
        hkl_tuple = hkl_info['hkl']
        hkl_labels.append("".join(map(str,hkl_tuple)))
        angle = two_thetas[i]
        intensity = intensities_broad[i]
    hkl_labels = ", ".join(hkl_labels)
    label_text = f'{hkl_labels}'
    plt.text(angle, intensity+1, label_text, rotation = 90, color = 'r', ha = 'center', va = 'bottom', fontsize = 8)

# Save xrd plot with broadening to Figures folder
broadening_filename = f"post_dft_broadening_xrd_{compound}.png"
plt.savefig(f"{full_path}/{broadening_filename}")

plt.show()

XRD Comparison Code:
Description: This code generates xrd spectra from the pre- and post-DFT .cif files and plots them on the same figure for convenient comparison. It creates an xrd plot of just the ideal peaks from pymatgen as well as a plot with some peak broadening to simulate a more realistic output that you would get from experiments.

# General python modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Materials science python modules
!pip install --upgrade pymatgen
import pymatgen
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Change directory to folder where you would like to create a folder to store the generated figures
target_directory = # Insert path to target folder
os.chdir(target_directory) # Change current directory to target folder

# Get cif files
pre_cif_file = '/content/drive/MyDrive/MS 508/MS508-SemesterProject /XRD/Structure 6/Pre DFT/rank05_KVO3.cif'
post_cif_file = '/content/drive/MyDrive/MS 508/MS508-SemesterProject /XRD/Structure 6/Post DFT/KVO3.cif'

# Create folder to save figures
figure_folder = 'Figures' # Create name for folder to hold generated figures
os.makedirs(figure_folder, exist_ok=True) # If folder does not exist, create folder
full_path = os.path.join(os.getcwd(), figure_folder) # Get path to created folder

# Get final structures from cif files
pre_structure = Structure.from_file(pre_cif_file)
post_structure = Structure.from_file(post_cif_file)

# Get string for compound
compound = # Insert structure compound name here

# Parameters
wl = 1.5418; # Wavelength for Cu K-alpha radiation, replace with desired radiation wavelength if necessary

# Initialize the XRD calculator (default wavelength is Cu K-alpha, lambda=1.5418 A)
xrd_calculator = XRDCalculator(wavelength=wl)

# Calculate the diffraction pattern pre dft
pre_pattern = xrd_calculator.get_pattern(pre_structure, scaled=True, two_theta_range=(0, 90))
pre_two_thetas = pre_pattern.x  # get 2theta values for peaks
pre_intensities = pre_pattern.y  # get intensity values for peaks
pre_d_hkls = pre_pattern.d_hkls

# Get peak hkls
pre_hkls = pre_pattern.hkls
pre_hkl_indices_only = []

# Calculate the diffraction pattern post dft
post_pattern = xrd_calculator.get_pattern(post_structure, scaled=True, two_theta_range=(0, 90))
post_two_thetas = post_pattern.x  # get 2theta values for peaks
post_intensities = post_pattern.y  # get intensity values for peaks
post_d_hkls = post_pattern.d_hkls

# Get peak hkls
post_hkls = post_pattern.hkls
post_hkl_indices_only = []

# Plot comparison with peaks only with labels
comparison_figure, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize =(12,8), sharey=True)
comparison_pre_plot = xrd_calculator.get_plot(pre_structure, ax=ax1, fontsize=10)
comparison_post_plot = xrd_calculator.get_plot(post_structure, ax=ax2, fontsize=10)
comparison_figure.suptitle(f"Pre- and Post- DFT XRD Pattern Comparison for {compound}", fontsize=20)
comparison_figure.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.35)
ax1.set_title('Pre DFT', fontsize=18)
ax2.set_title('Post DFT', fontsize=18)
ax1.set_xlabel('2$\\theta$ (degrees)', fontsize=16)
ax2.set_xlabel('2$\\theta$ (degrees)', fontsize=16)
ax1.set_ylabel('Intensity', fontsize=16)
ax2.set_ylabel('Intensity', fontsize=16)

# Save ideal xrd plot to Figures folder
ideal_filename = f"comparison_ideal_xrd_{compound}.png"
plt.savefig(f"{full_path}/{ideal_filename}")

# Show plot
plt.show()

# Plot comparison with peaks only without labels
comparison_figure, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize =(12,8), sharey=True)
comparison_pre_plot = xrd_calculator.get_plot(pre_structure, ax=ax1, fontsize=10, annotate_peaks=None)
comparison_post_plot = xrd_calculator.get_plot(post_structure, ax=ax2, fontsize=10, annotate_peaks=None)
comparison_figure.suptitle(f"Pre- and Post- DFT XRD Pattern Comparison for {compound}", fontsize=20)
comparison_figure.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace=0.35)
ax1.set_title('Pre DFT', fontsize=18)
ax2.set_title('Post DFT', fontsize=18)
ax1.set_xlabel('2$\\theta$ (degrees)', fontsize=16)
ax2.set_xlabel('2$\\theta$ (degrees)', fontsize=16)
ax1.set_ylabel('Intensity', fontsize=16)
ax2.set_ylabel('Intensity', fontsize=16)

# Save ideal xrd plot to Figures folder
ideal_filename = f"comparison_ideal_xrd_no_hkls{compound}.png"
plt.savefig(f"{full_path}/{ideal_filename}")

# Show plot
plt.show()

# Create broadening function
def apply_broadening(two_theta_ideal, intensities_ideal, fwhm, two_theta_range=(0, 90), step=0.01):
    """
    Applies Gaussian broadening to an ideal XRD pattern.

    Args:
        two_theta_ideal: Ideal 2-theta peak positions.
        intensities_ideal: Ideal peak intensities.
        fwhm: Full Width at Half Maximum for all peaks (constant broadening for simplicity).
        two_theta_range: The range for the output pattern.
        step: The step size for the output 2-theta array.
    """
    two_theta_broad = np.arange(two_theta_range[0], two_theta_range[1], step)
    broadened_intensities = np.zeros(two_theta_broad.shape)

    # Gaussian broadening formula
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2))) # Convert FWHM to sigma

    for two_t, intensity in zip(two_theta_ideal, intensities_ideal):
        # Add contribution of each peak to the overall pattern
        broadened_intensities += intensity * np.exp(-(two_theta_broad - two_t)**2 / (2 * sigma**2))

    return two_theta_broad, broadened_intensities

# Apply broadening
fwhm = 0.2
pre_two_theta_broad, pre_intensities_broad = apply_broadening(pre_two_thetas, pre_intensities, fwhm)
post_two_theta_broad, post_intensities_broad = apply_broadening(post_two_thetas, post_intensities, fwhm)

# Plot peaks with broadening without hkl labels
plt.figure(figsize=(10, 6))
plt.plot(pre_two_theta_broad, pre_intensities_broad, label=f'Pre-DFT Broadened Peaks (FWHM={fwhm}°)', color='darkorange')
plt.plot(post_two_theta_broad, post_intensities_broad, label=f'Post-DFT Broadened Peaks (FWHM={fwhm}°)', color='darkgreen')

plt.title(f"Comparison of Pre- and Post- DFT XRD Patterns with Peak Broadening for {compound}")
plt.xlabel('2$\\theta$ (degrees)')
plt.ylabel('Intensity')
plt.legend()
plt.xlim(0, 90)
plt.ylim(0, (max(max(pre_intensities), max(post_intensities))+20))
plt.grid(True)

# Save xrd plot with broadening to Figures folder
broadening_filename = f"comparison_broadening_xrd_no_hkls_{compound}.png"
plt.savefig(f"{full_path}/{broadening_filename}")

plt.show()

# Plot peaks with broadening with hkl labels
plt.figure(figsize=(10, 6))
plt.plot(pre_two_theta_broad, pre_intensities_broad, label=f'Pre-DFT Broadened Peaks (FWHM={fwhm}°)', color='darkorange')

# Add peak hkl labels for pre dft
for i in range(len(pre_two_thetas)):
    pre_hkl_list = pre_hkls[i]
    pre_hkl_labels = []
    for pre_hkl_info in pre_hkl_list:
        pre_hkl_tuple = pre_hkl_info['hkl']
        pre_hkl_labels.append("".join(map(str,pre_hkl_tuple)))
        pre_angle = pre_two_thetas[i]
        pre_intensity = pre_intensities_broad[i]
    pre_hkl_labels = ", ".join(pre_hkl_labels)
    pre_label_text = f'{pre_hkl_labels}'
    plt.text(pre_angle, pre_intensity+1, pre_label_text, rotation = 90, color = 'darkorange', ha = 'center', va = 'bottom', fontsize = 8)

plt.plot(post_two_theta_broad, post_intensities_broad, label=f'Post-DFT Broadened Peaks (FWHM={fwhm}°)', color='darkgreen')

# Add peak hkl labels for post dft
for i in range(len(post_two_thetas)):
    post_hkl_list = post_hkls[i]
    post_hkl_labels = []
    for post_hkl_info in post_hkl_list:
        post_hkl_tuple = post_hkl_info['hkl']
        post_hkl_labels.append("".join(map(str,post_hkl_tuple)))
        post_angle = post_two_thetas[i]
        post_intensity = post_intensities_broad[i]
    post_hkl_labels = ", ".join(post_hkl_labels)
    post_label_text = f'{post_hkl_labels}'
    plt.text(post_angle, post_intensity+1, post_label_text, rotation = 90, color = 'darkgreen', ha = 'center', va = 'bottom', fontsize = 8)

plt.title(f"Comparison of Pre- and Post- DFT XRD Patterns with Peak Broadening for {compound}")
plt.xlabel('2$\\theta$ (degrees)')
plt.ylabel('Intensity')
plt.legend()
plt.xlim(0, 90)
plt.ylim(0, (max(max(pre_intensities), max(post_intensities))+20))
plt.grid(True)

# Save xrd plot with broadening to Figures folder
broadening_filename = f"comparison_broadening_xrd_{compound}.png"
plt.savefig(f"{full_path}/{broadening_filename}")

plt.show()

Thanks, 
Abigail

--
Abigail Bemis (she/her)
PhD Student | Division of Materials Science and Engineering
Boston University
abemis@bu.edu
