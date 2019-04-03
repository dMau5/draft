from base64 import encodebytes
from CGRdb import load_schema, Molecule, Reaction
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
        info('find_similar')
        with db_session:
            similar_molecules = Molecule.find_similar(target, page=1, pagesize=1)
            info('done')
            for pair in similar_molecules:
                molecule = pair[0]
                if Molecule.reactions_entities(molecule, pagesize=10, product=True):
                    _target = molecule
                    info(f'Target = {molecule.structure}, Index Tanimoto = {pair[-1]}')
                    break
    if _target is None:
        info('Similar molecules not found as product')
        raise Exception

    sign = {_target}
    stack = [(_target, stages)]
    g = DiGraph()
    n = 0
    while stack:
        item = stack.pop(0)
        st = item[-1] - 1
        if st:
            with db_session:
                reactions = Molecule.reactions_entities(item[0], pagesize=first_number_of_paths, product=True)
                for r in reactions:
                    g.add_node(n, data=None)
                    for m in r.molecules:
                        if m not in sign:
                            number_of_mol = m.molecule
                            sign.add(number_of_mol)
                            structure = number_of_mol.structure
                            if not m.is_product:
                                stack.append((number_of_mol, st))
                                g.add_edge(structure, n)
                            g.add_edge(n, structure)
                            if bytes(structure) in zinc:
                                g.nodes[structure]['zinc'] = 1
        else:
            break
    return g
