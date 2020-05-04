from flow import FlowProject
import signac
import flow
import matplotlib.pyplot as plt
import mbuild as mb
import mdtraj as md
import numpy as np
from foyer import Forcefield
from get_sol_il_xml import GetSolv, GetIL, Get_ff_path
from mtools.gromacs.gromacs import make_comtrj
import os
import environment
import scipy.constants as constants
from ilforcefields.utils.utils import get_il, get_ff_path


def workspace_command(cmd):
    """Simple command to always go to the workspace directory"""
    return " && ".join(
        [
            "cd {job.ws}",
            cmd if not isinstance(cmd, list) else " && ".join(cmd),
            "cd ..",
        ]
    )


init_file = "system.gro"
em_file = "em.gro"
nvt_file = "nvt.gro"
npt_file = "npt.gro"
equil_file = "equil.gro"
sample_file = "sample.gro"
unwrapped_file = "sample_unwrapped.xtc"
rdf_file = "rdf-solvent-solvent.txt"


class Project(FlowProject):
    pass


@Project.label
def initialized(job):
    return job.isfile(init_file)


@Project.label
def minimized(job):
    return job.isfile(em_file)


@Project.label
def nvt_done(job):
    return job.isfile(nvt_file)


@Project.label
def npt_done(job):
    return job.isfile(npt_file)

@Project.label
def equil_done(job):
    return job.isfile(equil_file)


@Project.label
def sampled(job):
    return job.isfile(sample_file)


@Project.label
def prepared(job):
    return job.isfile(unwrapped_file)


@Project.label
def rdf_done(job):
    return job.isfile(rdf_file)


@Project.operation
@Project.post.isfile(init_file)
def initialize(job):
    with job:
        print(job.id())
        print("Setting up packing ...")
        il_conc = job.statepoint()[
            "concIL"
        ]  # mol/L , concIL is the concentration of either cation or anion
        acn_conc = job.statepoint()["concacn"]
        pack_dim = 6
        packing_box = mb.Box([pack_dim, pack_dim, pack_dim])
        unit_n = constants.Avogadro / (
            10 ** (24)
        )  # mol/L to  number of molecules / nm^3
        n_acn = round(acn_conc * unit_n * pack_dim ** (3))
        n_ions = round(il_conc * unit_n * pack_dim ** (3))

        # Load in mol2 files as mb.Compound
        anion = GetIL("tfsi")
        anion.name = "tfsi"
        cation = GetIL(job.statepoint()["cation"])
        cation.name = job.statepoint()["cation"]
        acn = GetSolv("acn")
        acn.name = "acn"

        il = mb.Compound()
        solvent = mb.Compound()
        for child in system.children:
            if child.name in [job.statepoint()["cation"], "tfsi"]:
                il.add(mb.clone(child))
            elif child.name == "acn":
                solvent.add(mb.clone(child))

        if n_ions == 0:
            system = mb.fill_box(compound=[acn], n_compounds=[n_acn], box=packing_box)
            opls = Forcefield(name="oplsaa")
            solventPM = opls.apply(solvent)
            print("Saving .gro, .pdb and .top ... ")
            systemPM.save("system.top", combine="all", overwrite=True)
            system.save("system.gro", combine="all", overwrite=True, residues=["acn"])
        elif n_acn == 0:
            system = mb.fill_box(
                compound=[cation, anion], n_compounds=[n_ions, n_ions], box=packing_box
            )
            ilPM = Forcefield(get_ff_path(kpl))
            il_system = ilPM.apply(il)
            il_system.save(
                "system.top",
                combine="all",
                overwrite=True,
            )
            system.save(
                "system.gro",
                combine="all",
                overwrite=True,
                residues=[job.statepoint()["cation"], "tfsi"],
            )
        else:
            system = mb.fill_box(
                compound=[acn, cation, anion],
                n_compounds=[n_acn, n_ions, n_ions],
                box=packing_box,
            )
            ilPM = "/raid6/homes/linx6/signac/oak_ridge/src/util/lib/kpl.xml"  # ions' forcefield xml file
            kpl = Forcefield(ilPM)
            il_system = kpl.apply(il)
            opls = Forcefield(name="oplsaa")
            solventPM = opls.apply(solvent)
            system_total = il_system + solventPM
            print("Saving .gro, .pdb and .top ... ")
            system_total.save(
                "system.top", combine="all", overwrite=True
            )  # for gromacs
            system.save(
                "system.gro",
                combine="all",
                overwrite=True,
                residues=[job.statepoint()["cation"], "tfsi", "acn"],
            )


@Project.operation
@Project.pre.isfile(init_file)
@Project.post.isfile(em_file)
@flow.cmd
def em(job):
    return _gromacs_str("em", "init", "init", job)


@Project.operation
@Project.pre.isfile(em_file)
@Project.post.isfile(nvt_file)
@flow.cmd
def nvt(job):
    return _gromacs_str("nvt", "em", "init", job)


@Project.operation
@Project.pre.isfile(nvt_file)
@Project.post.isfile(npt_file)
@flow.cmd
def npt(job):
    return _gromacs_str("npt", "nvt", "init", job)

@Project.operation
@Project.pre.isfile(npt_file)
@Project.post.isfile(equil_file)
@flow.cmd
def npt(job):
    return _gromacs_str("equil", "npt", "init", job)


@Project.operation
@Project.pre.isfile(equil_file)
@Project.post.isfile(sample_file)
@flow.cmd
def sample(job):
    return _gromacs_str("sample", "equil", "init", job)


@Project.operation
@Project.pre.isfile(sample_file)
@Project.post.isfile(unwrapped_file)
def prepare(job):
    trr_file = os.path.join(job.workspace(), "sample.trr")
    xtc_file = os.path.join(job.workspace(), "sample.xtc")
    gro_file = os.path.join(job.workspace(), "sample.gro")
    tpr_file = os.path.join(job.workspace(), "sample.tpr")
    if os.path.isfile(trr_file) and os.path.isfile(gro_file):
        unwrap_trj("sample.trr")
        unwrapped_trj = os.path.join(job.workspace(), "sample_unwrapped.trr")
        trj = md.load(unwrapped_trj, top=gro_file)

        comtrj = make_comtrj(trj)
        comtrj.save_xtc(com_trj)

        comtrj[-1].save_gro(os.path.join(job.workspace(), "com.gro"))
        print("made comtrj ...")
        # if os.path.isfile(com_trj) and not os.path.isfile(unwrapped_com_trj)    :
        os.system(
            "gmx trjconv -f {0} -o {1} -pbc nojump".format(com_trj, unwrapped_com_trj)
        )


@Project.operation
# @Project.pre.isfile(msd_file)
def run_rdf(job):
    print("Loading trj {}".format(job))
    if os.path.exists(os.path.join(job.workspace(), "com.gro")):
        top_file = os.path.join(job.workspace(), "com.gro")
        trj_file = os.path.join(job.workspace(), "sample_com.xtc")
        trj = md.load(trj_file, top=top_file, stride=10)

        selections = dict()
        selections["anion"] = trj.topology.select("name tfsi")
        selections["cation"] = trj.topology.select(
            "resname {}".format(job.statepoint()["cation"])
        )
        selections["acn"] = trj.topology.select("resname acn")

        combos = [
            ("cation", "anion"),
            ("cation", "cation"),
            ("anion", "anion"),
            ("acn", "anion"),
            ("acn", "cation"),
        ]
        for combo in combos:
            fig, ax = plt.subplots()
            print(
                "running rdf between {0} ({1}) and\t{2} ({3})\t...".format(
                    combo[0],
                    len(selections[combo[0]]),
                    combo[1],
                    len(selections[combo[1]]),
                )
            )
            r, g_r = md.compute_rdf(
                trj,
                pairs=trj.topology.select_pairs(
                    selections[combo[0]], selections[combo[1]]
                ),
                r_range=((0.0, 2.0)),
            )

            data = np.vstack([r, g_r])
            np.savetxt(
                os.path.join(
                    job.workspace(), "rdf-{}-{}.txt".format(combo[0], combo[1])
                ),
                np.transpose(np.vstack([r, g_r])),
                header="# r (nm)\tg(r)",
            )
            ax.plot(r, g_r)
            plt.xlabel("r (nm)")
            plt.ylabel("G(r)")
            plt.savefig(os.path.join(job.workspace(), f"rdf-{combo[0]}-{combo[1]}.pdf"))
            print(" ... done\n")


def _gromacs_str(op_name, gro_name, sys_name, job):
    """Helper function, returns grompp command string for operation """
    if op_name == "em":
        mdp = signac.get_project().fn("src/util/mdp_files/{}.mdp".format(op_name))
        cmd = "gmx grompp -f {mdp} -c system.gro -p system.top -o {op}.tpr --maxwarn 1 && gmx mdrun -deffnm {op} -ntmpi 1"
    else:
        mdp = signac.get_project().fn(
            "src/util/mdp_files/{}-{}.mdp".format(op_name, job.sp.T)
        )
        cmd = "gmx grompp -f {mdp} -c {gro}.gro -p system.top -o {op}.tpr --maxwarn 1 && gmx mdrun -deffnm {op} -ntmpi 1"
    return workspace_command(
        cmd.format(mdp=mdp, op=op_name, gro=gro_name, sys=sys_name)
    )


if __name__ == "__main__":
    Project().main()
