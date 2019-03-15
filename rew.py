import logging
import pickle
from pony.orm import db_session
from CGRtools.files import SDFread, RDFread, SDFwrite
from digraph import CGRdbDigraph
import networkx as nx
from multiprocessing import Process, Queue
from CGRdbUser import User
from CGRdb import load_schema, Molecule
from itertools import islice
logging.basicConfig(level=logging.ERROR)


# txt = 'US20040096399A1'
# n = 0
# with open('ex2.rdf', 'r', encoding='utf-8') as f:
#     reactions = RDFread(f)
#     for r in reactions:
#         for reactant in r.reactants:
#             _reactants = reactant.split()
#             if any(sum(atom.charge for t, atom in molecule.atoms()) for molecule in _reactants):
#                 pass
#             else:
#                 print('reactant', reactant)
#                 print('_reactants', [str(x) for x in _reactants])
#                 print('reagents', [str(x) for x in r.reagents])
#         for product in r.products:
#             _products = product.split()
#             if any(sum(atom.charge for t, atom in molecule.atoms()) for molecule in _products):
#                 pass
#             else:
#                 print('product', product)
#                 print('_products', [str(x) for x in _products])
#                 print('reagents', [str(x) for x in r.reagents])
#
# quit()

# with open('up_graph.pickle', 'rb') as gr:
#     print('load graph')
#     Gr = pickle.load(gr)
#     print('graph loaded')
# zinc = set()
# for x in range(1, 1872):
#     try:
#         with SDFread(f'zinc/{x}.sdf') as g:
#             print(x)
#             for mol in g:
#                 mol.standardize()
#                 zinc.add(bytes(mol))
#     except: continue
# with open('zinc.pickle', 'wb') as d:
#     pickle.dump(zinc, d)


def worker(input_queue):
    for r in iter(input_queue.get, 'STOP'):
        try:
            r.standardize()
            with db_session:
                for reactant in r.reactants:
                    _reactants = reactant.split()
                    if all(sum(atom.charge for t, atom in molecule.atoms()) for molecule in _reactants) \
                            or len(_reactants) == 1:
                        if not Molecule.structure_exists(reactant):
                            Molecule(reactant, User[1])
                    else:
                        other = []
                        for molecule in _reactants:
                            if not sum(atom.charge for t, atom in molecule.atoms()):
                                if not Molecule.structure_exists(molecule):
                                    Molecule(molecule, User[1])
                            else:
                                other.append(molecule)
                        for mol in other:
                            mol |= mol
                        if not Molecule.structure_exists(mol):
                            Molecule(mol, User[1])
                
                for product in r.products:
                    _products = product.split()
                    if any(sum(atom.charge for t, atom in molecule.atoms()) for molecule in _products):
                        if not Molecule.structure_exists(product):
                            Molecule(product, User[1])
                    else:
                        for molecule in _products:
                            if not Molecule.structure_exists(molecule):
                                Molecule(molecule, User[1])
        except:
            continue


if __name__ == '__main__':
    i = 0
    num = []
    db = load_schema('sandbox',)
    with open('final2.rdf', 'r', encoding='utf-8') as f:
        reactions = RDFread(f)
        print('work')
        # g = CGRdbDigraph()
        inp = Queue()
        for _ in range(12):
            Process(target=worker, args=(inp,)).start()
        for x in islice(reactions, 200):
            print(f'--------{i} done--------')
            inp.put(x)
            i += 1
        for _ in range(12):
            inp.put('STOP')
