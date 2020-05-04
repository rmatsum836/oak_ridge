#!/usr/bin/env python
"""Initialize the project's data space.

Iterates over all defined state points and initializes
the associated job workspace directories."""
import logging
import argparse
import signac
from sympy.solvers import solve
from sympy import Symbol

MW = {'emim_tfsi': 391.31, 'bmim_tfsi': 419.36, 'hmim_tfsi': 447.42, 'acn': 41.05} # g/mol
density = {'emim_tfsi': 1.5235, 'bmim_tfsi': 1.4431, 'hmim_tfsi': 1.37, 'acn': 0.786} # g/cm^(-3)
# How many Liters for 1 mol emim_tfsi, and unit has been changed to L/mol
E_1mol_V = MW['emim_tfsi']/density['emim_tfsi']*10**(-3) #L/mol
B_1mol_V = MW['bmim_tfsi']/density['bmim_tfsi']*10**(-3)
H_1mol_V = MW['hmim_tfsi']/density['hmim_tfsi']*10**(-3)
acn_1mol_V = MW['acn']/density['acn']*10**(-3)
concentration = [0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 10] # x mol IL in 24.36 mol acn

x = Symbol('x')
ES_conc =[] #[[IL mol/L, acn mol/L],,,,] for emim_tfsi in acn at different concentration
BS_conc =[]
HS_conc =[]  
for IL_1mol_V in [E_1mol_V,B_1mol_V,H_1mol_V]:
    con_n = 0
    IL_conc = []
    for con_n in range(len(concentration)):
        if concentration[con_n] == 0:
            IL_conc.append([0,1/acn_1mol_V])
        elif concentration[con_n] == 10:
            IL_conc.append([1/IL_1mol_V,0])
        else:
            # Calculate how many mols of IL should be in 1 L volume based on the concentration ratio, x mol* E_1mol_V + x * 24.36/concentration[con_n] * acn_1mol_V = 1L
            s = solve(x * IL_1mol_V + x * 24.36/concentration[con_n] * acn_1mol_V - 1, x) 
            IL_conc.append([s[0],s[0] * 24.36/concentration[con_n]]) # IL_conc[[IL conc, acn conc], ....]
        if IL_1mol_V==E_1mol_V:
            ES_conc=IL_conc
        elif IL_1mol_V==B_1mol_V:
            BS_conc=IL_conc
        else:
            HS_conc=IL_conc   


def main(args):
    project = signac.init_project('oak_ridge_ray')
    statepoints_init = []
    temperatures = [283, 288, 293, 298, 303, 308, 313, 318, 323, 328, 333]
    cations= ['emim', 'bmim', 'hmim']
    for cation in cations:
        for temp in temperatures:
            if cation == 'emim':
                for conc in ES_conc:
                    statepoint = dict(
                            cation=cation,
                            T=temp,
                            concIL=float(conc[0]),
                            concacn=float(conc[1])
                            )
                    project.open_job(statepoint).init()
                    statepoints_init.append(statepoint)
            if cation == 'bmim':
                for conc in BS_conc:
                    statepoint = dict(
                            cation=cation,
                            T=temp,
                            concIL=float(conc[0]),
                            concacn=float(conc[1])
                            )
                    project.open_job(statepoint).init()
                    statepoints_init.append(statepoint)
            if cation == 'hmim':
                for conc in HS_conc:
                    statepoint = dict(
                            cation=cation,
                            T=temp,
                            concIL=float(conc[0]),
                            concacn=float(conc[1])
                            )
                    project.open_job(statepoint).init()
                    statepoints_init.append(statepoint)
    # Writing statepoints to hash table as a backup
    project.write_statepoints(statepoints_init)


if __name__ == '__main__':
     parser = argparse.ArgumentParser(
         description="Initialize the data space.")
     parser.add_argument(
         '-n', '--num-replicas',
         type=int,
         default=1,
         help="Initialize multiple replications.")
     args = parser.parse_args()
  
     logging.basicConfig(level=logging.INFO)
     main(args)

        
