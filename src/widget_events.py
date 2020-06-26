import glob
import re
import ipywidgets as widgets
from IPython.display import clear_output

import mbuild as mb
import warnings
warnings.filterwarnings('ignore')

# Read in and create a dictionary of all the compounds
files = glob.glob('./molecules-mol2/*.mol2')
compounds = {'Select your molecule': None}
reverse_path = dict()

for file in files:
    molname = re.sub('.mol2','',re.sub('^.*/','',file))
    compound = mb.load(file, use_parmed=True)
    compounds[molname] = compound
    reverse_path[compound] = file

# Define some globale variable
# These variables will be update when the widget event handlers are called
COMPOUND = None
BOX_OF_COMPOUNDS = None

# Define some widget styles and layout
style = {'description_width': 'initial'}

# Create all the widgets
compound_dropdown = widgets.Dropdown(options=compounds,
                                     description='Select your molecule',
                                     style=style)
box_slider = widgets.FloatSlider(min=1, max=10,
                description="Dimension of box (nm)",
                value=5, style=style)
n_slider = widgets.IntSlider(min=1, max=100,
                description="Number of particles",
                value=50, style=style)

# Define output

out_mol = widgets.Output()
out_box = widgets.Output()

def visualize(compound):
    # A method to filter out numeric string in particle name
    # before being visualize by py3dmol (so the color scheme
    # could be consistent)
    vis_compound = mb.clone(compound)
    for particle in vis_compound.particles():
        particle.name = ''.join(filter(str.isalpha,
                                    particle.name))

    display(vis_compound.visualize())

def compound_handler(obj):
    global COMPOUND
    COMPOUND = compound_dropdown.value
    with out_mol:
        clear_output()
        path = reverse_path[COMPOUND]
        print('compound = mb.load("{}")'.format(path))
        print("compound.visualize()")
        if compound_dropdown.value:
            visualize(compound_dropdown.value)


def box_handler(obj):
    box = [box_slider.value] * 3
    n_compounds = n_slider.value
    compound = compound_dropdown.value

    box_of_compound = mb.fill_box(compound=compound,
                                  n_compounds=n_compounds,
                                  box=box)
    global BOX_OF_COMPOUNDS
    BOX_OF_COMPOUNDS = box_of_compound
    with out_box:
        clear_output()
        print("box = {}".format(box))
        print("box_of_compound = mb.fill_box(compound=compoound, "
              "n_compounds={}, "
              "box={})".format(n_compounds, box))
        print("box_of_compound.visualize()")
        visualize(box_of_compound)

# Define widget event trigger
compound_dropdown.observe(compound_handler, names='value')
compound_dropdown.observe(box_handler, names='value')
box_slider.observe(box_handler, names='value')
n_slider.observe(box_handler, names='value')


# Experiment with SMILES string

out_smiles = widgets.Output()

smiles_box = widgets.Text(
    placeholder='Enter a SMILES string',
    description='SMILES:',
    disabled=False,
    style=style)

def smiles_handler(obj):
    try:
        compound = mb.load(smiles_box.value, smiles=True)
    except:
        compound = None
    with out_smiles:
        if compound:
            clear_output()
            print('smiles_compound = mb.load("{}", smiles=True)'.format(
                                                        smiles_box.value))
            print('smiles_compound.visualize()')
            visualize(compound)
        else:
            clear_output()
            print('Invalid SMILES string')

smiles_box.observe(smiles_handler, names='value')
