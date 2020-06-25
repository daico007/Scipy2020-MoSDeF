import mbuild as mb
import foyer
import antefoyer
import gmso
import gmso.external
import gmso.formats
from run_simulation import run_energy_minimization, run_nvt


# LOADING IN COMPOUND TO MBUILD
# Note: we will replace this with the jupyter widget later.
smiles_string = input("Enter SMILES String: ")

compound = mb.load(smiles_string, smiles=True)
system = mb.fill_box(compound, n_compounds=10, density=10)

ff = foyer.forcefields.load_GAFF()

structure = ff.apply(system)

def apply_charges(structure, compound, n_compounds):
    single_mol_struct = ff.apply(compound)
    single_mol_struct_charge = antefoyer.ante_charges(
            single_mol_struct, 'bcc', net_charge=0.00, multiplicity=1)
    for index, atom in enumerate(structure.atoms):
        atom.charge = single_mol_struct_charge.atoms[index%n_compounds].charge
    return structure

structure_charge = apply_charges(structure, compound, 10)

topology = gmso.external.from_parmed(structure_charge)
topology.name = compound.name

gmso.formats.write_top(topology, "simulation/topol.top")
gmso.formats.write_gro(topology, "simulation/conf.gro")


