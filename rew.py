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


# db = load_schema()
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
# with SDFwrite('zinc_in_graph.sdf') as z:
#     with db_session:
#         for x in Gr.nodes():
#             if 'zinc' in Gr.nodes[x]:
#                 z.write(x)
#
# print('num')
# for x in range(1, 13):
#     with SDFread(f'zinc/{x}.sdf') as f1:
#         print(x)
#         for y in f1:
#             num += 1
#             if mol == y:
#                 print('have')
#         print('zinc', num)
# num = 0
# with db_session:
#     for x in Gr.nodes():
#         if 'zinc' in Gr.nodes[x]:
#             num += 1
# print('num', num)


def worker(input_queue):
    for r in iter(input_queue.get, 'STOP'):
        try:
            r.standardize()
            with db_session:
                for reactant in r.reactants:
                    # if isinstance(reactant, MoleculeContainer):
                    _reactants = reactant.split()
                    if all(sum(atom.charge for molecule in _reactants for n, atom in molecule.atoms())):
                        for molecule in _reactants:
                            if not Molecule.structure_exists(molecule):
                                Molecule(molecule, User[1])
                    else:
                        if not Molecule.structure_exists(reactant):
                            Molecule(reactant, User[1])
                    # g.add_edge(reactant, n)
                for product in r.products:
                    if not Molecule.structure_exists(product):
                        Molecule(product, User[1])
        except:
            continue


if __name__ == '__main__':
    i = 0
    db = load_schema('sandbox')
    with open('reactions_db.rdf', 'r', encoding='utf-8') as f:
        reactions = RDFread(f)
        print('work')
        # g = CGRdbDigraph()
        inp = Queue()
        for _ in range(12):
            Process(target=worker, args=(inp,)).start()
        for x in reactions:
            print(f'--------{i} done--------')
            inp.put(x)
            i += 1
        for _ in range(12):
            inp.put('STOP')

# with open('up_graph.pickle', 'wb') as f:
#     pickle.dump(Gr, f)
# print('finish')
# with db_session:
#     with open('zinc/log2.txt', 'w') as d:
#         for node in Gr.nodes():
#             if isinstance(node, int):
#                 number = Gr.node[node]['reaxys_data'][0]
#                 if number == bad:
#                     d.write(f'node {node}\n')
#                     d.write(f'bad {bad}\n')
#                     d.write(f'predec {[str(x) for x in Gr._pred[node]]}\n')
#                     d.write(f'synok {[str(x) for x in Gr._succ[node]]}\n')
#                 elif number == bad1:
#                     d.write(f'node {node}\n')
#                     d.write(f'bad {bad1}\n')
#                     d.write(f'predec {[str(x) for x in Gr._pred[node]]}\n')
#                     d.write(f'synok {[str(x) for x in Gr._succ[node]]}\n')
#                 else:
#                     continue


#             if txt not in text:
#                 text.update({txt: n})
#             else:
#                 Gr.remove_node(node)
#         else:
#             continue


