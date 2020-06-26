from distutils.spawn import find_executable
from subprocess import PIPE, Popen
from os import path

GMX = find_executable("gmx")

def run_energy_minimization(wd="./simulation/"):
    print("Preprocessing energy minimization input...")
    mdpfile = path.join(wd, "em.mdp")
    topfile = path.join(wd, "topol.top")
    grofile = path.join(wd, "conf.gro")
    try:
        assert all([path.exists(f) for f in [mdpfile, topfile, grofile]])
    except:
        raise FileNotFoundError(f"em.mdp, topol.top, and conf.gro were "
                                 f"not found in {wd}")
    grompp_command = (f"gmx grompp -f em.mdp -c conf.gro -p topol.top -maxwarn 2")
    proc = Popen(grompp_command, stdout=PIPE, stderr=PIPE,
                 universal_newlines=True, shell=True, cwd=wd)

    out, err = proc.communicate()

    if "Fatal error" in err:
        print(err)
        raise RuntimeError("Gromacs GROMMP command failed")

    print("Running energy minimization...")
    mdrun_command = "gmx mdrun"

    proc = Popen(mdrun_command, stdout=PIPE, stderr=PIPE,
                 universal_newlines=True, shell=True, cwd=wd)

    out, err = proc.communicate()


def run_nvt(wd="./simulation/"):
    print("Preprocessing simulation input...")
    mdpfile = path.join(wd, "nvt.mdp")
    topfile = path.join(wd, "topol.top")
    grofile = path.join(wd, "confout.gro")
    try:
        assert all([path.exists(f) for f in [mdpfile, topfile, grofile]])
    except:
        raise FileNotFoundError(f"nvt.mdp, topol.top, and confout.gro were "
                                 f"not found in {wd}")
    grompp_command = (f"gmx grompp -f nvt.mdp -c confout.gro -p topol.top -maxwarn 2")
    proc = Popen(grompp_command, stdout=PIPE, stderr=PIPE,
                 universal_newlines=True, shell=True, cwd=wd)

    out, err = proc.communicate()

    if "Fatal error" in err:
        print(err)
        raise RuntimeError("Gromacs GROMMP command failed")

    print("Running simulation...")
    mdrun_command = "gmx mdrun -v"

    proc = Popen(mdrun_command, stdout=PIPE, stderr=PIPE,
                 universal_newlines=True, shell=True, cwd=wd)
    while True:
        output = proc.stderr.readline()
        if output == '' and proc.poll() is not None:
            break
        if output:
            if ("will finish" in output or
                "remaining wall" in output):
                print(output.strip(), end="\r")
            else:
                pass
                #print(output.strip())

    print()
    rc = proc.poll()


def visualize_trajectory(wd="./simulation/"):
    xtcfile = path.join(wd, "traj_comp.xtc")
    wrappedxtcfile = path.join(wd, "wrapped.xtc")
    grofile = path.join(wd, "confout.gro")
    try:
        assert all([path.exists(f) for f in [xtcfile, grofile]])
    except:
        raise FileNotFoundError(f"traj_comp.xtc, and confout.gro were "
                                 f"not found in {wd}")
    grompp_command = (f"echo 0 | gmx trjconv -f traj_comp.xtc -o wrapped.xtc -pbc whole")
    proc = Popen(grompp_command, stdout=PIPE, stderr=PIPE,
                 universal_newlines=True, shell=True, cwd=wd)

    out, err = proc.communicate()

    if "Fatal error" in err:
        print(err)
        raise RuntimeError("Gromacs TRJCONV command failed")

    import nglview as nv
    import mdtraj as md
    traj = md.load(wrappedxtcfile, top=grofile)
    ngl_widget = nv.show_mdtraj(traj)
    return ngl_widget

