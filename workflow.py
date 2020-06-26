import mbuild as mb
import foyer
import antefoyer

import gmso
import gmso.external
import gmso.formats

from run_simulation import run_energy_minimization, run_nvt

import panedr

# System building with mbuild
path_to_mol2 = input("Enter the path to your molecule file (str): ")
box_size = float(input("Enter your box size (int): "))
n_compounds = input("Enter the numbe of compounds in the box (int): ")


compound = mb.load(path_to_mol2)
box_of_compounds = mb.fill_box(compound=compound,
                    n_compounds=int(n_compounds),
                    box=[box_size]*3)

# Atomtyping with foyer (GAFF forcefield)

def apply_charges(box_structure, single_compound, n_atoms, ff):
    single_mol_struct_charge = antefoyer.ante_charges(
            single_compound, 'bcc', net_charge=0.00, multiplicity=1)

    for index, atom in enumerate(box_structure.atoms):
        atom.charge = single_mol_struct_charge.atoms[index%n_atoms].charge
    return box_structure

gaff_ff = foyer.forcefields.load_GAFF()
typed_compound = gaff_ff.apply(box_of_compounds,
                          assert_dihedral_params=False)

charge_structure = apply_charges(box_structure=typed_compound,
                                 single_compound=compound,
                                 n_atoms=compound.n_particles,
                                 ff=gaff_ff)

# Handed back to our backend GMSO
topology = gmso.external.from_parmed(charge_structure)
topology.name = compound.name

gmso.formats.write_top(topology, "simulation/topol.top",
                       top_vars={"fudgeLJ": 0.5, "fudgeQQ": 0.8, "comb-rule": "geometric"})

gmso.formats.write_gro(topology, "simulation/conf.gro")

# Run the simulation with gromacs
run_energy_minimization()
run_nvt()

# Data analysis
sim_data = panedr.edr_to_df("simulation/ener.edr")

