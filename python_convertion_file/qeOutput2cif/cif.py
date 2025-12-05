# ==============================
# QE to CIF Converter (Final Relaxed Structure)
# ==============================

# Import ASE functions for reading and writing structures
from ase.io import read, write

# -----------------------------
# Step 1: Read QE output file
# -----------------------------
# Replace "pw.out" with your QE output filename.
# ASE automatically reads the last frame by default, which corresponds
# to the final relaxed structure.
structure = read("geom.out")  

# If you want to be extra sure, you can read all frames and pick the last one:
# all_steps = read("pw.out", index=":")  # read all frames
# structure = all_steps[-1]               # select the last step (final geometry)

# -----------------------------
# Step 2: Write to CIF
# -----------------------------
# Writes the structure to a CIF file.
write("final_output.cif", structure)

# -----------------------------
# Step 3: Confirmation message
# -----------------------------
print("CIF file successfully written as 'final_output.cif'.")
print("This represents the final relaxed geometry from your QE calculation.")
