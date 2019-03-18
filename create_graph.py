from CGRdb import load_schema
from CGRtools.files import SDFwrite, RDFread
from digraph import CGRdbDigraph
import pickle
from multiprocessing import Process, Queue
from time import sleep
from itertools import islice
import logging
logging.basicConfig(level=logging.ERROR)
with open('zinc.pickle', 'rb') as d:
    zinc = pickle.load(d)


# def worker(input_queue, output_queue):
#     for r in iter(input_queue.get, 'STOP'):
#         try:
#             r.standardize()
#             if r.reagents and r.products:
#                 output_queue.put(r)
#         except:
#             continue
db = load_schema('sandbox',)
with open('final2.rdf', 'r', encoding='utf-8') as f:
    reactions = RDFread(f)
    print('work')
    g = CGRdbDigraph()
    n = 0
    for r in reactions:
        try:
            r.standardize()
            if r.reactants and r.products:
                g.add_node(n, data=[r.meta['source_id'], r.meta['text']])
                for reactant in r.reactants:
                    _reactants = reactant.split()
                    if all(sum(atom.charge for t, atom in molecule.atoms()) for molecule in _reactants) \
                            or len(_reactants) == 1:
                        g.add_edge(reactant, n)
                        if bytes(reactant) in zinc:
                            g.nodes[reactant]['zinc'] = True
                    else:
                        other = []
                        for molecule in _reactants:
                            if not sum(atom.charge for t, atom in molecule.atoms()):
                                g.add_edge(molecule, n)
                                if bytes(molecule) in zinc:
                                    g.nodes[molecule]['zinc'] = True
                            else:
                                other.append(molecule)
                        er = other.pop()
                        for mol in other:
                            er |= mol
                        g.add_edge(er, n)
                        if bytes(er) in zinc:
                            g.nodes[er]['zinc'] = True
                for product in r.products:
                    _products = product.split()
                    if all(sum(atom.charge for t, atom in molecule.atoms()) for molecule in _products) \
                            or len(_products) == 1:
                        g.add_edge(n, product)
                        if bytes(product) in zinc:
                            g.nodes[product]['zinc'] = True
                    else:
                        other = []
                        for molecule in _products:
                            if not sum(atom.charge for t, atom in molecule.atoms()):
                                g.add_edge(n, molecule)
                                if bytes(molecule) in zinc:
                                    g.nodes[molecule]['zinc'] = True
                            else:
                                other.append(molecule)
                        er = other.pop()
                        for mol in other:
                            er |= mol
                        g.add_edge(n, er)
                        if bytes(er) in zinc:
                            g.nodes[er]['zinc'] = True
                n += 1
                print(f'{n} done')
        except:
            continue
with open('experimental_graph.pickle', 'wb') as file:
    pickle.dump(g, file)
    print('URAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
# if __name__ == '__main__':
#     db = load_schema('sandbox', )
#     with open('final2.rdf', 'r', encoding='utf-8') as f:
#         reactions = RDFread(f)
#         print('work')
#         g = CGRdbDigraph()
#         inp = Queue()
#         n = 0
#         for x in islice(reactions, 100):
#             inp.put(x)
#         out = Queue()
#         for _ in range(12):
#             Process(target=worker, args=(inp, out)).start()
#         while True:
#             sleep(1)
#             q1 = out.qsize()
#             q2 = inp.qsize()
#             if q1:
#                 # print('do output', q1, q2)
#                 for _ in range(q1):
#                     r = out.get()
#                     g.add_node(n, reaxys_data=[r.meta['source_id'], r.meta['text']])
#                     for reagent in r.reactants:
#                         g.add_edge(reagent, n)
#                         if bytes(reagent) in zinc:
#                             g.nodes[reagent]['zinc'] = True
#                     for product in r.products:
#                         g.add_edge(n, product)
#                         if bytes(product) in zinc:
#                             g.nodes[product]['zinc'] = True
#                     n += 1
#                     if not n % 1000:
#                         print(f'{n} done')
#                     try:
#                         inp.put(next(reactions))
#                     except StopIteration:
#                         break
#                 else:
#                     continue
#                 break
#             elif not q2:
#                 # print('do input')
#                 for e in islice(reactions, 100):
#                     inp.put(e)
#
#         for _ in range(12):
#             inp.put('STOP')

