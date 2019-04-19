from base64 import encodebytes
from CGRdb import load_schema, Molecule
from CGRdb.database import MoleculeReaction
from json import dumps
from logging import info, basicConfig, INFO
from networkx import DiGraph
from pickle import load
from pony.orm import db_session
from CGRtools.containers import MoleculeContainer

def substructure(target):
    molecules = list()
    rctns = None
    substructures = Molecule._structure_query(target, 'substructure')
    _target = Molecule.find_structure(target)
    if _target:
        if MoleculeReaction.exists(molecule=_target, is_product=True):
            molecules.append((_target, 1, _target.reactions_entities(pagesize=100, product=True)))
    for pair in substructures:
        m = None
        sub, tan = pair
        if sub != _target:
            for mapping in target.get_substructure_mapping(sub.structure, limit=0):
                if MoleculeReaction.exists(molecule=sub, is_product=True):
                    if rctns and len(molecules) < 10 and m:
                        molecules.append((sub, round(tan, 2), rctns))
                    reactions = sub.reactions_entities(pagesize=100, product=True)
                    rctns = []
                    for r in reactions:
                        for mppng in select(x for x in MoleculeReaction.entity if x.molecule == sub and x.reaction == r):
                            mppng = mppng.mapping
                            if not set(r.cgr.center_atoms).isdisjoint(mppng.get(x, x) for x in mapping.values()):
                                m = sub
                                rctns.append(r.id)

    if molecules:
        return molecules
    else:
        info('Substructures not found as product')
        raise Exception
