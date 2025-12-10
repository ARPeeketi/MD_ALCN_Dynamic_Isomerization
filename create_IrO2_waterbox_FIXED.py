"""
IrO2(110) Surface with Water Box and Ions - FIXED VERSION
=============================================================
Uses proper molify packing with cropping approach:
1. Pack in larger box at correct density
2. Crop to exact dimensions
3. Include Na+ directly in packing

Creates:
- IrO2(110) slab (4x4 supercell, 5 layers)
- 8 adsorbates (*O, *OH, *H2O) at Ir atop sites
- ~15 Å water box with correct density
- Na+ ion (packed with water)

Requirements:
- ase (pip install ase)
- molify (pip install molify)
- numpy
"""

from ase.build import surface, add_adsorbate, molecule
from ase.io import write
from ase import Atom, Atoms
import numpy as np

# 1. Create IrO2 (110) slab
print("Creating IrO2(110) slab...")
alat, clat = 4.498, 3.154
irO2_bulk = Atoms('Ir2O4',
                  cell=[(alat, 0, 0), (0, alat, 0), (0, 0, clat)],
                  scaled_positions=[
                      (0.0, 0.0, 0.0),      # Ir
                      (0.5, 0.5, 0.5),      # Ir
                      (0.3, 0.3, 0.0),      # O
                      (0.7, 0.7, 0.0),      # O
                      (0.2, 0.8, 0.5),      # O
                      (0.8, 0.2, 0.5),      # O
                  ],
                  pbc=True)
slab_unit = surface(irO2_bulk, (1, 1, 0), layers=5, vacuum=15.0)
slab = slab_unit.repeat((4, 4, 1))
slab.center(axis=2, vacuum=15.0)

# 2. Add 8 adsorbates at Ir atop sites
print("Adding adsorbates...")
adsorbates = ['O', 'OH', 'H2O', 'O', 'OH', 'OH', 'OH', 'O']
adz_height = 2.0

positions = slab.get_positions()
elements = slab.get_chemical_symbols()
top_z = np.percentile(positions[:, 2], 95)
ir_indices = [i for i, el in enumerate(elements) if el == 'Ir' and positions[i, 2] >= top_z]
top_sites = [positions[i] for i in ir_indices]
top_sites = sorted(top_sites, key=lambda p: (round(p[1], 1), round(p[0], 1)))[:8]

for i, site in enumerate(top_sites):
    kind = adsorbates[i % len(adsorbates)]
    if kind == 'O':
        add_adsorbate(slab, 'O', adz_height, position=site[:2])
    elif kind == 'OH':
        oh = molecule('OH')
        oh.rotate(180, 'x')
        oh.translate(site + [0, 0, adz_height])
        slab += oh
    elif kind == 'H2O':
        h2o = molecule('H2O')
        h2o.translate(site + [0, 0, adz_height])
        slab += h2o

# 3. Pack water + Na+ in larger box, then crop
from molify import pack, smiles2conformers

# Target dimensions
cell = slab.get_cell()
target_x = cell[0, 0]  # ~25.44 Å
target_y = cell[1, 1]  # ~12.62 Å  
target_z = 15.0        # Water height

z_top_slab = slab.get_positions()[:, 2].max()
z_bot = z_top_slab + 1.0  # Start 1 Å above slab

# Pack in larger box (1.3x in all dimensions) - smaller for faster packing
pack_factor = 1.3
pack_x = target_x * pack_factor
pack_y = target_y * pack_factor
pack_z = target_z * pack_factor

# Calculate molecules for LARGER box at density = 1000 kg/m³
volume_m3 = (pack_x * pack_y * pack_z) * 1e-30  # Å³ to m³
water_mw = 18.015e-3  # kg/mol
na_mw = 22.990e-3
avogadro = 6.022e23

# Mostly water, a few Na+ ions
n_na = 5  # Pack more Na+ to ensure 1-2 survive crop
na_mass = n_na * na_mw / avogadro
water_mass = 1000 * volume_m3 - na_mass  # Total mass - ion mass
n_waters = int(water_mass / (water_mw / avogadro))

print(f"\nPacking {n_waters} H2O + {n_na} Na+ in box: {pack_x:.1f} x {pack_y:.1f} x {pack_z:.1f} Å")

# Generate molecules
waters = [smiles2conformers("O", 1)[0] for _ in range(n_waters)]
na_ions = [Atoms([Atom('Na', position=[0, 0, 0])]) for _ in range(n_na)]

# Pack together at correct density
print("Running PACKMOL...")
packed = pack([waters, na_ions], [n_waters, n_na], 
              density=1000,  # kg/m³
              ratio=(1.0, pack_y/pack_x, pack_z/pack_x),
              tolerance=2.0,
              verbose=False)

print(f"Packed: {len(packed)} atoms")

# The packed box will have approximately the right aspect ratio
# We don't scale - we just crop to exact dimensions

# Shift box to start at origin
pos = packed.get_positions()
pos -= pos.min(axis=0)
packed.set_positions(pos)

# Use CENTERED crop to maximize chance of keeping ions
pos = packed.get_positions()
box_size = pos.max(axis=0) - pos.min(axis=0)

# Center of packed box
center = pos.min(axis=0) + box_size / 2

# Define target box centered around packed box center
x_min = center[0] - target_x / 2
x_max = center[0] + target_x / 2
y_min = center[1] - target_y / 2
y_max = center[1] + target_y / 2
z_min = 0  # Keep z from bottom
z_max = target_z

mask = ((pos[:, 0] >= x_min) & (pos[:, 0] <= x_max) &
        (pos[:, 1] >= y_min) & (pos[:, 1] <= y_max) &
        (pos[:, 2] >= z_min) & (pos[:, 2] <= z_max))

water_box = packed[mask]

# Shift cropped box to origin
pos = water_box.get_positions()
pos[:, 0] -= x_min
pos[:, 1] -= y_min
water_box.set_positions(pos)

print(f"After cropping to {target_x:.1f} x {target_y:.1f} x {target_z:.1f} Å: {len(water_box)} atoms")

# Check Na+ survived
symbols = water_box.get_chemical_symbols()
n_na_final = symbols.count('Na')
print(f"Na+ ions after crop: {n_na_final}")

# Set correct cell and PBC
water_box.set_cell([target_x, target_y, target_z])
water_box.set_pbc([True, True, False])

# Position above slab
water_box.translate([0, 0, z_bot])

# 4. Merge slab + water
final_cell = slab.get_cell().copy()
z_total = water_box.get_positions()[:, 2].max() - slab.get_positions()[:, 2].min() + 2.0
final_cell[2, 2] = z_total

final_atoms = slab + water_box
final_atoms.set_cell(final_cell)
final_atoms.set_pbc([True, True, False])
final_atoms.center(axis=2, vacuum=10.0)
final_atoms.wrap()

# 5. Export
write("IrO2_110_waterbox.data", final_atoms, format="lammps-data")
write("IrO2_110_waterbox.xyz", final_atoms, format="extxyz")
write("IrO2_110_waterbox.cif", final_atoms, format="cif")

print(f"\n✓ System created:")
print(f"  Total atoms: {len(final_atoms)}")
print(f"  Slab atoms: {len(slab)}")
print(f"  Water box atoms: {len(water_box)}")
print(f"  Cell: {final_atoms.get_cell().cellpar()[:3]}")

# 6. Save ClO4- template
try:
    clo4 = smiles2conformers("[O-]Cl(=O)(=O)=O", 1)[0]
    write("clo4.xyz", clo4)
    write("clo4.pdb", clo4)
    print("✓ ClO4- template saved")
except Exception as e:
    print(f"Note: ClO4- generation failed: {e}")
