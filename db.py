from CGRtools import SDFRead, SDFWrite, RDFWrite
from multiprocessing import Process, Queue
from CGRdb import load_schema, Molecule, Reaction
from CGRdbData import ReactionConditions
from logging import warning, basicConfig, ERROR
from itertools import islice, chain, cycle, product, filterfalse
from collections import defaultdict
from math import ceil
from random import shuffle, sample, choice
from os import environ, listdir
from pony.orm import db_session
from pickle import load, dump
from CGRdb import Molecule, load_schema, Reaction
from pony.orm import db_session
from operator import or_
from functools import reduce
from CGRtools.containers import ReactionContainer
from CGRtools import RDFRead
db = load_schema('name', user='user', password='pass', host='host', database='database')

with open('starting_blocks_patents.pickle', 'rb') as f2:
    zinc = load(f2)
with open('set_10_binar', 'rb') as f3:
    rules = load(f3)


def worker(u, o):
    for chn in iter(u.get, 'STOP'):
        chns = []
        for y in chn:
            rct = reduce(or_, y.reactants)
            pdt = reduce(or_, y.products)
            for rule in rules:
                r_rct = rule.reactants[0]
                r_pdt = rule.products[0]
                if r_rct < rct and r_pdt < pdt:
                    chns.append(1)
                    break
        if len(chns) == len(chn):
            o.put(True)
        else:
            o.put(False)


if __name__ == '__main__':
    inp = Queue()
    out = Queue()
    for _ in range(15):
        Process(target=worker, args=(inp, out)).start()

    pths = 0
    usego = 0
    dr = set()
    with open('drugs_ids.pickle', 'rb') as f:
        drugs = load(f)
        for i, m_id in enumerate(drugs):
            n = 0
            if m_id in dr:
                continue
            try:
                with open(f'drugs_paths/chains_{m_id}.pickle', 'rb') as w:
                    chains = load(w)
            except:
                continue
            dr.add(m_id)
            print(i, 'powel')
            print('len', len(chains))
            print('drugs -->', len(dr))
            print('chains -->', usego)
            print('chains in rules -->', pths)
            seen = set()
            for path in chains:
                check = [int(x.split('_')[0]) for x in path if isinstance(x, str)]
                flag = True
                reactions = []
                with db_session:
                    for num, r in enumerate(Reaction.select(lambda x: x.id in check)):
                        if not num and r.id not in seen:
                            usego += 1
                            seen.add(r.id)
                            rctnts = []
                            for m in r._molecules:
                                if not m.is_product:
                                    rctnts.append(m.id)
                            if not all(x in zinc for x in rctnts):
                                flag = False
                        if flag:
                            reactions.append(r.structure)
                        else:
                            break
                if reactions:
                    n += 1
                    inp.put(reactions)
            for _ in range(n):
                if out.get():
                    pths += 1
    print('drugs -->', len(dr))
    print('chains -->', usego)
    print('chains in rules -->', pths)
    for _ in range(15):
        inp.put('STOP')

# with open('green/agrees.pickle', 'rb') as f1:
#     agrees = load(f1)
# dd = defaultdict(dict)
# with db_session:
#     for r in Reaction.select(lambda x: x.id in agrees):
#         pr = []
#         re = []
#         for m in r._molecules:
#             obj_mol = m.molecule
#             if not m.is_product:
#                 re.append(obj_mol.id)
#             else:
#                 pr.append(obj_mol.id)
#         for i in re:
#             for p in pr:
#                 dd[i].update({p: r.id})
# with open(f'green/graph_{2}.pickle', 'rb') as w:
#     f = load(w)
# with open('from_reactant_to_product.pickle', 'rb') as w:
#     reac_prod = load(w)

# ms = (3, 4, 5, 6, 7, 12, 13, 22)
# a = listdir('green')
# with RDFWrite('reactions_for_green_reactants.rdf') as w:
#     for name in a:
#         if 'graph' in name:
#             print(name)
#             i = name.split('_')[1].split('.')[0]
#             with open(f'green/{name}', 'rb') as f:
#                 dd = {x: 0 for x in range(1, 11)}
#                 g, source = load(f)
#                 seen = set()
#                 reactions = []
#                 for node in g.nodes():
#                     if isinstance(node, str):
#                         print(node)
#                         res = node.split('_')
#                         with db_session:
#                             r = Reaction[int(res[1])].structure
#                         r.meta['stage_in_synthesis_tree'] = res[0]
#                         w.write(r)
#             #         for suc in g._succ[node]:
#             #             if suc not in seen:
#             #                 dd[st] += 1
#             #             seen.add(suc)
#             # with open(f'green/{i}_st_stats', 'wb') as w:
#             #     dump((dd, source), w)
# exit()
# n = 0
# mols = set()
# with open('green/greens.pickle', 'rb') as f:
#     greens = load(f)
# for i, mol in enumerate(greens):
#     print(i, 'powel')
#     prods = set()
#     with db_session:
#         mol_db = Molecule.find_structure(mol)
#         stages = 10
#         if mol_db:
#             mol_id = mol_db.id
#             g = DiGraph()
#             seen = {mol_id}
#             # added_reactions = set()
#             stack = [(mol_db, stages)]
#             if mol_id in reac_prod:
#                 stages -= 1
#                 for prod_id, r_id in reac_prod[mol_id].items():
#                     r_node = f'{abs(stages-10)}_{r_id}'
#                     g.add_edge(mol_id, r_node)
#                     g.add_edge(r_node, prod_id)
#                     if prod_id not in seen:
#                         stack.append((prod_id, stages))
#                     seen.add(prod_id)
#
#                 while stack:
#                     mt, st = stack.pop(0)
#                     st -= 1
#                     if mt in reac_prod:
#                         for p, r in reac_prod[mt].items():
#                             r_node = f'{abs(st-10)}_{r}'
#                             g.add_edge(mt, r_node)
#                             g.add_edge(r_node, p)
#                             if p not in seen and st:
#                                 stack.append((p, st))
#                             seen.add(p)
#
#             if g:
#                 n += 1
#                 with open(f'green/graph_{i}.pickle', 'wb') as w:
#                     dump((g, mol_id), w)
# print('usego', n)

# # for j in ms:
# #     with open(f'green/{j}/graph_{j}_first_3st.pickle', 'rb') as f:
# #         products = load(f)
#     for i, mol in enumerate(greens):
#         # if mol in mols:
#         #     continue
#         print(i, 'powel')
#         prods = set()
#         with db_session:
#             mol_db = Molecule.find_structure(mol)
#             stages = 10
#             if mol_db:
#                 mol_id = mol_db.id
#                 g = DiGraph()
#                 seen = {mol_id}
#                 added_reactions = set()
#                 stack = [(mol_db, stages)]
#                 while stack:
#                     mt, st = stack.pop(0)
#                     st -= 1
#                     reactions = mt.reactions_entities(pagesize=10000, product=False)
#                     for r in reactions:
#                         r_id = r.id
#                         if r_id in added_reactions or r_id not in agrees:
#                             continue
#                         added_reactions.add(r.id)
#                         r_node = f'{r_id}_{abs(st - 10)}'
#                         for m in r._molecules:
#                             obj_mol = m.molecule
#                             structure = obj_mol.structure
#                             if not m.is_product:
#                                 if st and obj_mol not in seen:
#                                     seen.add(obj_mol)
#                                     stack.append((obj_mol, st))
#                                 g.add_edge(obj_mol.id, r_node)
#                             else:
#                                 if not st:
#                                     prods.add(obj_mol.id)
#                                 g.add_edge(r_node, obj_mol.id)
#
#                 if g:
#                     n += 1
#                     with open(f'green/graph_{i}_first_3st.pickle', 'wb') as w:
#                     # with open(f'green/{j}/succ/graph_{mol_id}_first_4_6st.pickle', 'wb') as w:
#                         dump((g, mol_db.structure), w)
#         if prods:
#             with open(f'green/for_{i}_th_mol_products_after_3st.pickle', 'wb') as w:
#                 dump(prods, w)
#     print('usego', n)
