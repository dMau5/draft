from base64 import encodebytes
from CGRdb import load_schema, Molecule, Reaction
from CGRdb.database import MoleculeReaction
from CGRdbUser import User
from json import dumps
from logging import info
from networkx import DiGraph
from pickle import load
from pony.orm import db_session

load_schema('patents', )
with open('zinc.pickle', 'rb') as z:
    zinc = load(z)


def b64(m):
    pass


def visualization(G, target):
    pass


def search_synth_paths(target, stages=3, first_number_of_paths=3):
    target.standardize()
    target.aromatize()
    _target = Molecule.find_structure(target)
    if _target is None:
        with db_session:
            info('find_similar')
            similar_molecules = Molecule.find_similar(target, page=1, pagesize=10)
            info('done')
            for pair in similar_molecules:
                molecule = pair[0]
                if MoleculeReaction.exists(molecule=molecule, product=True):
                    _target = molecule
                    info(f'Target = {molecule.structure}, Index Tanimoto = {pair[-1]}')
                    break
    if _target is None:
        info('Similar molecules not found as product')
        raise Exception

    seen = {_target}
    stack = [(_target, stages)]
    g = DiGraph()
    n = 0
    added_nodes = {}
    while stack:
        mt, st = stack.pop(0)
        st -= 1
        with db_session:
            reactions = mt.reactions_entities(pagesize=first_number_of_paths, product=True)
            for r in reactions:
                if r.id not in added_nodes:
                    g.add_node(n, data=None)
                    added_nodes[r.id] = n
                    n += 1
                for m in r.molecules:
                    obj_mol = m.molecule
                    if obj_mol not in seen:
                        seen.add(obj_mol)
                        structure = obj_mol.structure
                        if not m.is_product:
                            if st:
                                stack.append((obj_mol, st))
                            g.add_edge(structure, n)
                        else:
                            g.add_edge(n, structure)
                        if bytes(structure) in zinc:
                            g.nodes[structure]['zinc'] = 1
    return g
