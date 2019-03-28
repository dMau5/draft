
import pickle
from pony.orm import db_session
from CGRtools.containers import ReactionContainer
from CGRtools.files import SDFread, RDFread, SDFwrite, RDFwrite
from digraph import CGRdbDigraph
import networkx as nx
from multiprocessing import Process, Queue
from CGRdbUser import User
from CGRdb import load_schema, Molecule
from logging import warning, basicConfig, ERROR
from itertools import islice
from CGRtools.algorithms import union
# basicConfig(level=ERROR)

db = load_schema('sandbox',)
w = SDFwrite('exists.sdf')
we = RDFwrite('similar.rdf')
zinc = set()
for x in range(1, 1872):
    try:
        with SDFread(f'zinc/{1}.sdf', indexable=True) as f:
            print(f'do zinc/{x}.sdf')
            l = len(f)
            _3 = f[3]
            for z_mol in f:
                z_mol.standardize()
                z_mol.aromatize()
                z_mol.implicify_hydrogens()
                zinc.add(bytes(z_mol))
                with db_session:
                    if Molecule.structure_exists(z_mol):
                        print('exists')
                        w.write(z_mol)
                        continue
                    mols = Molecule.find_similar(z_mol, page=1, pagesize=1)
                    if mols and mols[0][-1] >= .95:
                        print('similar')
                        mol = [mols[0][0].structure]
                        r = ReactionContainer(reactants=[z_mol], products=mol)
                        we.write(r)
    except:
        warning(f'open  ne mojet')
with open('zinc.pickle', 'wb') as z:
    pickle.dump(zinc, z)

# with open('2_5_graph.pickle', 'rb') as gr:
#     print('load graph')
#     Gr = pickle.load(gr)
#     print('graph loaded')
# n = 0
# with db_session:
#     for node in Gr.nodes():
#         if not isinstance(node, int):
#             if 'zinc' in node:
#                 n += 1
# print('n', n)
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


# def worker(input_queue):
#     for r in iter(input_queue.get, 'STOP'):
#         try:
#             r.standardize()
#             with db_session:
#                 for reactant in r.reactants:
                    # _reactants = reactant.split()
                    # if all(sum(atom.charge for t, atom in molecule.atoms()) for molecule in _reactants) \
                    #         or len(_reactants) == 1:
                    # if not Molecule.structure_exists(reactant):
                    #     Molecule(reactant, User[1])
                    # else:
                    #     other = []
                    #     for molecule in _reactants:
                    #         if not sum(atom.charge for t, atom in molecule.atoms()):
                    #             if not Molecule.structure_exists(molecule):
                    #                 Molecule(molecule, User[1])
                    #         else:
                    #             other.append(molecule)
                    #     er = other.pop()
                    #     if len(other) > 1:
                    #         for mol in other:
                    #             er |= mol
                    #     if not Molecule.structure_exists(er):
                    #         Molecule(er, User[1])

                # for product in r.products:
                    # _products = product.split()
                    # if all(sum(atom.charge for t, atom in molecule.atoms()) for molecule in _products) \
                    #         or len(_products) == 1:
                    # if not Molecule.structure_exists(product):
                    #     Molecule(product, User[1])
                    # else:
                    #     other = []
                    #     for molecule in _products:
                    #         if not sum(atom.charge for t, atom in molecule.atoms()):
                    #             if not Molecule.structure_exists(molecule):
                    #                 Molecule(molecule, User[1])
                    #         else:
                    #             other.append(molecule)
                    #     er = other.pop()
                    #     if len(other) > 1:
                    #         for mol in other:
                    #             er |= mol
                    #     if not Molecule.structure_exists(er):
                    #         Molecule(er, User[1])
        # except:
        #     continue


# if __name__ == '__main__':
#     i = 1000000
#     num = []
#     db = load_schema('old',)
#     with open('reactions_db.rdf', 'r', encoding='utf-8') as f:
#         print('init')
#         reactions = RDFread('reactions_db.rdf', indexable=True)
#         print('start seek')
#         reactions.seek(1000000)
#         print('len', len(reactions))
#         print('work')
#         # g = CGRdbDigraph()
#         inp = Queue()
#         for _ in range(12):
#             Process(target=worker, args=(inp,)).start()
#         for x in reactions:
#             if not i % 1000:
#                 print(f'--------{i} done--------')
#             inp.put(x)
#             i += 1
#         for _ in range(12):
#             inp.put('STOP')
