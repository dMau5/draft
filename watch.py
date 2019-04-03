from base64 import encodebytes
from CGRdb import load_schema, Molecule, Reaction
from CGRdb.database import MoleculeReaction
from CGRtools.files import SDFread
from CGRdbUser import User
from json import dumps
from logging import info, INFO, basicConfig
from networkx import DiGraph
from pickle import load
from pony.orm import db_session

basicConfig(level=INFO)
load_schema('patents', )
with open('zinc.pickle', 'rb') as z:
    zinc = load(z)


def b64(m):
    return encodebytes(m.depict().encode()).decode().replace('\n', '')


def search_synth_paths(target, stages=3, first_number_of_paths=3):
    target.standardize()
    target.aromatize()
    error = False
    print(target)
    with db_session:
        _target = Molecule.find_structure(target)
        if _target is None:
            info('find_similar')
            similar_molecules = Molecule.find_similar(target, page=1, pagesize=10)
            info('done')
            for pair in similar_molecules:
                molecule, tan = pair[0]
                print(molecule.structure, tan)
                if MoleculeReaction.exists(molecule=molecule, is_product=True):
                    info(f'Target = {molecule.structure}, Index Tanimoto = {tan}')
                    _target = molecule
                    break
            else:
                error = True
        else:
            if not MoleculeReaction.exists(molecule=_target, is_product=True):
                error = True

        if error:
            info('Molecule-target and similar molecules not found as product')
            raise Exception

        seen = {_target}
        stack = [(_target, stages)]
        g = DiGraph()
        n = 0
        added_nodes = {}
        while stack:
            mt, st = stack.pop(0)
            print(mt, st)
            st -= 1
            reactions = mt.reactions_entities(pagesize=first_number_of_paths, product=True)
            print('reactions', reactions)
            for r in reactions:
                if r.id not in added_nodes:
                    data = [r.metadata][0].data[0]
                    g.add_node(n, data=f"{data['source_id']}, {data['text']}")
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
    with open('index.html', 'w') as w:
        w.write(visualization(g, target))


with SDFread('s.sdf') as f:
    mol = next(f)
    search_synth_paths(mol)
