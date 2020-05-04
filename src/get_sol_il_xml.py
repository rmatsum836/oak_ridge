import inspect
import os
import mbuild as mb
from pkg_resources import resource_filename

def GetSolv(solv_name):
    cache_dir = '/raid6/homes/linx6/signac/oak_ridge/src/util/lib/mol2/'
    filename = '{}.mol2'.format(solv_name)
    if any(file == filename for file in os.listdir(cache_dir)):
        solv = mb.load(os.path.join(cache_dir, filename))
        solv.name = solv_name

    return solv

def GetIL(il_name):
    cache_dir = '/raid6/homes/linx6/signac/oak_ridge/src/util/lib/mol2/'
    filename = '{}.mol2'.format(il_name)
    if any(file == filename for file in os.listdir(cache_dir)):
        il = mb.load(os.path.join(cache_dir, filename))
        il.name = il_name

    return il

def Get_ff_path(ff_name):
    """Get the path to a force field xml file """
    """in a directory of the same name."""
    cache_dir = '/raid6/homes/linx6/signac/oak_ridge/src/util/lib/'
    ff_path = os.path.join(cache_dir, ff_name + '.xml')
    return ff_path