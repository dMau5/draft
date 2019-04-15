from base64 import encodebytes
from CGRdb import load_schema, Molecule
from CGRdb.database import MoleculeReaction
from json import dumps
from logging import info, basicConfig, INFO
from networkx import DiGraph
from pickle import load
from pony.orm import db_session
from CGRtools.containers import MoleculeContainer

basicConfig(level=INFO)
load_schema('patents', user='postgres', password='jyvt0n3', host='localhost', database='postgres')
with open('zinc.pickle', 'rb') as z:
    zinc = load(z)


def search_synthetic_paths(target, stages=3, number_of_paths=3, substructure=None):
    target.standardize()
    target.aromatize()
    error = False
    with db_session:
        _target = Molecule.find_structure(target)
        if _target is None:
            if not substructure:
                info('find_similar')
                similar_molecules = Molecule.find_similar(target, page=1, pagesize=10)
                info('done')
                for pair in similar_molecules:
                    molecule, tan = pair
                    if MoleculeReaction.exists(molecule=molecule, is_product=True):
                        _target = molecule
                        break
                else:
                    error = True
            else:
                substructures = Molecule._structure_query(target, substructure)
                for sub in substructures:
                    if target.get_substructure_mapping(sub.structure):
                        if MoleculeReaction.exists(molecule=sub, is_product=True):
                            _target = sub
                            break
        else:
            tan = None
            molecule = None
            if not MoleculeReaction.exists(molecule=_target, is_product=True):
                error = True

        if error:
            info('Molecule-target and similar molecules not found as product')
            raise Exception

        seen = {_target}
        stack = [(_target, stages)]
        g = DiGraph()
        n = 0
        added_reactions = set()
        while stack:
            mt, st = stack.pop(0)
            st -= 1
            reactions = mt.reactions_entities(pagesize=number_of_paths, product=True)
            for r in reactions:
                if r.id not in added_reactions:
                    data = list(r.metadata)[0].data
                    g.add_node(n, data=f"{data['source_id']}, {data['text']}".replace('+\n', ''))
                    g.add_edge(n, mt.structure)
                    added_reactions.add(n)
                else:
                    continue
                for m in r.molecules:
                    obj_mol = m.molecule
                    structure = obj_mol.structure
                    if not m.is_product:
                        if st and obj_mol not in seen:
                            seen.add(obj_mol)
                            stack.append((obj_mol, st))
                        g.add_edge(structure, n)
                    else:
                        g.add_edge(n, structure)
                    if bytes(structure) in zinc:
                        g.nodes[structure]['zinc'] = 1
                n += 1

    html = visualization(g, _target.structure)
    info('Finish Him!')
    return html, molecule, tan
