import pickle
from pony.orm import db_session, select
from CGRtools.containers import ReactionContainer, MoleculeContainer
from CGRtools.files import SDFread, RDFread, SDFwrite, RDFwrite
import networkx as nx
from multiprocessing import Process, Queue
from CGRdbUser import User
from CGRdb import load_schema, Molecule, Reaction
from CGRdb.database import MoleculeStructure, MoleculeReaction
from CGRdbData import ReactionConditions, MoleculeProperties
from logging import warning, basicConfig, ERROR
from itertools import islice, chain, cycle, product, filterfalse
from CGRdb.search.fingerprints import FingerprintMolecule
from time import sleep, time
from collections import defaultdict
from math import ceil
from random import shuffle, sample
from os import environ
import numpy as np

# environ["PATH"] += ':/home/dinar/balls'
load_schema('all_patents',)


def chunks(iterable, size=10):
    iterator = iter(iterable)
    for first in iterator:
        yield chain([first], islice(iterator, size - 1))


with open('pairs_from_reactant_to_product.pickle', 'rb') as f:
    pairs = pickle.load(f)

with open('fingerprints.pickle', 'rb') as f1:
    fps = pickle.load(f1)

with open('list_of_A.pickle', 'rb') as f0:
    sigma = pickle.load(f0)

with open('similarities_array.pickle', 'rb') as f2:
    big = pickle.load(f2)

with open('old_to_new.pickle', 'rb') as f3:
    old_to_new = pickle.load(f3)

with open('new_to_old.pickle', 'rb') as f4:
    new_to_old = pickle.load(f4)


def evaluation(query, res):
    """
    оценка нод.
    возвращает танимото для пары запрос-результат.
    """
    # query_fp, res_fp = (FingerprintMolecule.get_fingerprint(x) for x in [query, res])
    qc, rc, common = len(query), len(res), len(query.intersection(res))
    return common / (qc + rc - common)


def roundrobin1(*iterables):
    num_active = len(iterables)
    iterables = cycle(iterables)
    while num_active:
        try:
            for nxt in iterables:
                yield next(nxt)
        except StopIteration:
            num_active -= 1
            iterables = cycle(islice(iterables, num_active))


def roundrobin2(iterables):
    dt = []
    for tu_sh, tu, A_sh in iterables:
        try:
            yield next(tu_sh), tu, A_sh
            dt.append((tu_sh, tu, A_sh))
        except StopIteration:
            pass

    num_active = len(dt)
    iterables = cycle(dt)
    while num_active:
        try:
            for tu_sh, tu, A_sh in iterables:
                yield next(tu_sh), tu, A_sh
        except StopIteration:
            num_active -= 1
            iterables = cycle(islice(iterables, num_active))


def unique_everseen(iterable, key=None):
    "List unique elements, preserving order. Remember all elements ever seen."
    # unique_everseen('AAAABBBCCDAABBB') --> A B C D
    # unique_everseen('ABBCcAD', str.lower) --> A B C D
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in filterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element


def worker(p, o):
    for m_A, m_A_sh, t_u_sh in iter(p.get, 'STOP'):
        fp_1 = fps[t_u_sh]
        t_unics = set(pairs[m_A]) - set(pairs[m_A_sh])
        index = 0
        for t_u in t_unics:
            fp_2 = fps[t_u]
            new_ind = evaluation(fp_1, fp_2)
            if new_ind > index:
                index = new_ind
        with db_session:
            m_A = Molecule[m_A].structure
            t_u_sh = Molecule[t_u_sh].structure
        o.put((m_A, t_u_sh, False, index))


if __name__ == '__main__':
    with open('False_pairs.txt', 'w') as fw:
        inp = Queue()
        out = Queue()
        for _ in range(20):
            Process(target=worker, args=(inp, out, )).start()

        file = 0
        rank = 1000
        data = []
        # holost = 15700
        for i, id_1 in enumerate(sigma):
            print(f'-- {i} powel --')
            new_id = old_to_new[id_1]

            tshki_1 = set(pairs[id_1])
            total = len(tshki_1)
            print('total', total)
            if total > 10000:
                tshki_1 = set(sample(tshki_1, 10000))
            # cut T to 10000
            tsh01 = []
            tsh45 = []
            tsh90 = []
            for n, x in enumerate(big[new_id]):
                if new_id != n:
                    if x <= 25.5:
                        tsh = set(pairs[new_to_old[n]])
                        if len(tsh) > 10000:
                            tsh = set(sample(tsh, 10000))
                        tsh01.append((tsh, new_to_old[n]))
                    elif 102 <= x <= 127.5:
                        tsh = set(pairs[new_to_old[n]])
                        if len(tsh) > 10000:
                            tsh = set(sample(tsh, 10000))
                        tsh45.append((tsh, new_to_old[n]))
                    elif x >= 229.5:
                        tsh = set(pairs[new_to_old[n]])
                        if len(tsh) > 10000:
                            tsh = set(sample(tsh, 10000))
                        tsh90.append((tsh, new_to_old[n]))

            t_sh01 = ((iter(tsh - tshki_1), tshki_1 - tsh, idd) for tsh, idd in tsh01)
            t_sh45 = ((iter(tsh - tshki_1), tshki_1 - tsh, idd) for tsh, idd in tsh45)
            t_sh90 = ((iter(tsh - tshki_1), tshki_1 - tsh, idd) for tsh, idd in tsh90)

            tot = ceil(total / 10000)
            for chunk in chunks(cycle((t_sh01, t_sh45, t_sh90)), 10000):
                tot -= 1
                ll = len(list(chunk))
                for triple in unique_everseen(roundrobin2(roundrobin1(*chunk)), lambda z: z[0]):
                    t_unic_sh, t_unic, m_A_sh = triple
                    inp.put((id_1, m_A_sh, t_unic_sh))
                for _ in range(ll):
                    print(f'for {i}  ||', 'inp-->', inp.qsize(), 'out-->', out.qsize(), 'rank-->', rank, 'file-->', file)
                    dt = out.get()
                    data.append(dt)
                    fw.write(f'{str(dt[0])}, {str(dt[1])}, False, {str(dt[3])}\n')
                    # holost -= 1
                    # if holost <= 0:
                    rank -= 1
                    if not rank:
                        with open(f'False_pairs/{file}.pickle', 'wb') as wq:
                            pickle.dump(data, wq)
                        data = []
                        file += 1
                        rank = 1000
                if not tot:
                    break
            # print('prowlo', time() - ini, 'sec')

        for _ in range(20):
            inp.put('STOP')

            # t_sh01 = ((iter(set(pairs[new_to_old[n]]) - tshki_1), tshki_1 - set(pairs[new_to_old[n]]), new_to_old[n])
            #           for n, x in enumerate(big[new_id]) if x <= 25.5 and new_id != n)
            # t_sh45 = ((iter(set(pairs[new_to_old[n]]) - tshki_1), tshki_1 - set(pairs[new_to_old[n]]), new_to_old[n])
            #           for n, x in enumerate(big[new_id]) if 102 <= x <= 127.5 and new_id != n)
            # t_sh90 = ((iter(set(pairs[new_to_old[n]]) - tshki_1), tshki_1 - set(pairs[new_to_old[n]]), new_to_old[n])
            #           for n, x in enumerate(big[new_id]) if x >= 229.5 and new_id != n)

# all_mols = set(pairs)
# for k, values in pairs.items():
#     for y in values:
#         all_mols.add(y)
#
# print('usego', len(all_mols))
# fps = {}
# n = 0
# with db_session:
#     for chunk in chunks(list(Molecule.select(lambda x: x.id in all_mols)), size=1000):
#         n += 1
#         for mol in chunk:
#             fp = FingerprintMolecule.get_fingerprint(mol.structure)
#             fps[mol.id] = fp
#         print('done -->', n)

# dict odnostadiynyh reactions: key == reactant from building blocks, value == products from patents reactions
# with open('from_prod_to_react.pickle', 'rb') as er:
#     df = pickle.load(er)

# sborka defaultdict: key == ids reactant from patents, value == dicts of key = product from patents, value = stages
# n = 0
# pairs = defaultdict(dict)
# trues = []
# for x in range(1, 589):
#     with open(f'bb/{x}.pickle', 'rb') as f:
#         bbs = pickle.load(f)
#         print('---', x * 1000, 'go ---')
#     for mol in bbs:
#         with db_session:
#             mol_db = Molecule.find_structure(mol)
#             # number of stages
#             stages = 10
#             if mol_db:
#                 mol_id = mol_db.id
#                 if mol_id in df:
#                     n += 1
#                     values = df[mol_id]
#                     stack = []
#                     seen = set()
#                     stages -= 1
#                     for prod in values:
#                         # update defauldict
#                         pairs[mol_id].update({prod: abs(stages - 10)})
#                         if prod not in seen:
#                             stack.append((prod, stages))
#                         seen.add(prod)
#
#                     while stack:
#                         mt, st = stack.pop(0)
#                         st -= 1
#                         if mt in df:
#                             val = df[mt]
#                             for pro in val:
#                                 # update defauldict
#                                 pairs[mol_id].update({pro: abs(st - 10)})
#                                 if pro not in seen:
#                                     if st:
#                                         stack.append((pro, st))
#                                 seen.add(pro)
#
#     print('--- finish na ---')
#     print('--> usego', n)

# dump defauldict of pairs
# with open('pairs_from_reactant_to_product.pickle', 'wb') as bds:
#     pickle.dump(pairs, bds)
#
#
# with open('old_to_new_ids.pickle', 'rb') as f1:
#     old_to_new = pickle.load(f1)
#
# with open('array_of_similarity_A.pickle', 'rb') as f2:
#     big = pickle.load(f2)
#
# with open('list_of_A.pickle', 'rb') as we:
#     sigma = pickle.load(we)
#
#
# def diff(t, n):
#     upd = 0
#     diff = len(t) - n
#     if diff < 0:
#         upd = - diff
#     return upd
#
#
# def writer(m_A, t1, t2, num):
#     unic_tshki_1 = t1 - t2
#     print('started mols', len(unic_tshki_1), 'to', num)
#     unic_tshki_2 = t2 - t1
#     print('u_tshki_2--->', len(unic_tshki_2))
#     ts = list(unic_tshki_2)
#     shuffle(ts)
#     for t_2 in ts[:num + 1]:
#         with open(f'database/{t_2}.pickle', 'rb') as p2:
#             fp_2 = pickle.load(p2)
#         ind = 0
#         for t_1 in unic_tshki_1:
#             with open(f'database/{t_1}.pickle', 'rb') as p1:
#                 fp_1 = pickle.load(p1)
#             new_ind = evaluation(fp_2, fp_1)
#             if new_ind > ind:
#                 ind = new_ind
#         fw.write(ReactionContainer(reactants=[m_A],
#                                    products=[Molecule[t_2].structure], meta={'fer': ind / 2}))


# with RDFwrite('False_pairs.rdf') as fw:
#     for i, mol_id_1 in enumerate(sigma):
#         print(f'-- {i} powel --')
#         ini = time()
#         with db_session:
#             mol_A = Molecule[mol_id_1].structure
#             high = middle = low = 0
#             high_mol_id = middle_mol_id = low_mol_id = 0
#             for mol_id_2 in sigma:
#                 if mol_id_1 == mol_id_2:
#                     continue
#                 ind_tan = big[old_to_new[mol_id_1]][old_to_new[mol_id_2]] / 255
#                 if .9 <= ind_tan <= 1:
#                     if ind_tan > high:
#                         high_mol_id = mol_id_2
#                         high = ind_tan
#                 elif .4 <= ind_tan <= .5:
#                     if ind_tan > middle:
#                         middle_mol_id = mol_id_2
#                         middle = ind_tan
#                 elif 0 < ind_tan <= .1:
#                     if ind_tan > low:
#                         low_mol_id = mol_id_2
#                         low = ind_tan
#
#             tshki_1 = pairs[mol_id_1]
#             total = len(tshki_1)
#             try:
#                 tshki_2_h = pairs[high_mol_id]
#                 tshki_2_m = pairs[middle_mol_id]
#                 tshki_2_l = pairs[low_mol_id]
#                 hi = len(tshki_2_h)
#                 mi = len(tshki_2_m)
#                 lo = len(tshki_2_l)
#                 to = hi + mi + lo
#                 hi = hi / to
#                 mi = mi / to
#                 lo = lo / to
#             except Exception as e:
#                 print(e)
#             if high_mol_id:
#                 print('high', 100 * hi)
#                 writer(mol_A, tshki_1, tshki_2_h, ceil(total * hi))
#             if middle_mol_id:
#                 print('middle', 100 * mi)
#                 writer(mol_A, tshki_1, tshki_2_m, ceil(total * mi))
#             if low_mol_id:
#                 print('low', lo * 100)
#                 writer(mol_A, tshki_1, tshki_2_l, ceil(total * lo))
#             print('prowlo', time() - ini)
#             print(f'-- finish {i} --\n')

# sigma = set(pairs)
# sigma = list(sigma)
# print(type(sigma))
# print(sigma)
# print(len(sigma))
#
# with open('list_of_A.pickle', 'wb') as w0:
#     pickle.dump(sigma, w0)
#
# big = np.eye(len(sigma), dtype=np.uint8)
# old_to_new = {sigma[0]: 0}
# new_to_old = {0: sigma[0]}
# n = 0
# for i, id_1 in enumerate(sigma):
#     print(f'--- {i} powel ---')
#     fp_1 = fps[id_1]
#     for id_2 in sigma[i+1:]:
#         fp_2 = fps[id_2]
#         if not i:
#             n += 1
#             old_to_new[id_2] = n
#             new_to_old[n] = id_2
#         if not n % 100:
#             print(f'--- uwlo {n} ---')
#         ind = evaluation(fp_1, fp_2)
#         new_id_2 = old_to_new[id_2]
#         big[i][new_id_2] = big[new_id_2][i] = ind * 255
# print('finish')


# select all data of reactions from database and save pairs of tuples (ids reactant, ids product)
# reaction_id = 0
# yse = list()
# rs = dict()
# with db_session:
#     ale = select((x.reaction.id, x.molecule.id, x.is_product) for x in MoleculeReaction.entity)\
#         .order_by(1)
#     print(ale.count())
#     b = ceil(ale.count() / 1000)
#     new_t = []
#     for i in range(1, b + 1):
#         print(i * 1000, 'powel')
#         dt = ale.page(i, 1000)
#         for t in dt:
#             # create tuples for every reaction (t[0] = id of reaction)
#             if reaction_id == t[0]:
#                 new_t.append(t)
#             else:
#                 if new_t:
#                     rs[reaction_id] = new_t
#                 reaction_id = t[0]
#                 new_t = [t]
#
#     print('dict powel')
#     for v in rs.values():
#         reacts = []
#         prodts = []
#         [prodts.append(m[1]) if m[2] else reacts.append(m[1]) for m in v]
#         # create all pairs [id reactant, id product]
#         yse.extend(list(product(reacts, prodts)))
# with open('right_pairs.pickle', 'rb') as rp:
#     yse = pickle.load(rp)

# create one step reactions: key == reactant from patents, value == set of products from patents
# dd = defaultdict(set)
# for pair in yse:
#     dd[pair[0]].add(pair[1])
# with open('from_prod_to_react.pickle', 'wb') as er:
#     pickle.dump(dd, er)

# def trees(m):
#     if m in pairs:
#         new_tshki = pairs[m]
#         new_tshki -= tshki
#     else:
#         new_tshki = set()
#         m_db = Molecule.find_structure(m)
#         if m_db:
#             new_seen = {m_db}
#             new_stack = [(m_db, stages)]
#             new_added_reactions = set()
#             print('create_new_graph')
#             while new_stack:
#                 mt, st = new_stack.pop(0)
#                 reactions = mt.reactions_entities(pagesize=10000, product=False)
#                 st -= 1
#                 for r in reactions:
#                     if r.id in new_added_reactions:
#                         continue
#                     new_added_reactions.add(r.id)
#                     for mm in r.molecules:
#                         obj_mol = mm.molecule
#                         structure = obj_mol.structure
#                         if mm.is_product:
#                             pairs[m].add(structure)
#                             if structure not in tshki:
#                                 new_tshki.add(structure)
#                             if st:
#                                 if obj_mol not in new_seen:
#                                     new_seen.add(obj_mol)
#                                     new_stack.append((obj_mol, st))
#             print('new_done')
#         for t in tshki:
#             ind = 0
#             s = 0
#             for n_t in new_tshki:
#                 new_ind = evaluation(t, n_t)
#                 if new_ind > ind:
#                     ind = new_ind
#                     s = n_t
#             if s:
#                triples.add((m, t, ind / 2))


# with open('sigmaldrich.pickle', 'rb') as f:
#     gh = pickle.load(f)
#     with RDFwrite('true_pairs_1.rdf') as w:
#         for i, z_mol in enumerate(gh[86:], 86):
#             print(i, 'go')
#             tshki = set()
#             with db_session:
#                 print('search')
#                 mol_db = Molecule.find_structure(z_mol)
#                 if mol_db:
#                     print('create graph')
#                     seen = {mol_db}
#                     stages = 10
#                     stack = [(mol_db, stages)]
#                     added_reactions = set()
#                     while stack:
#                         mt, st = stack.pop(0)
#                         reactions = mt.reactions_entities(pagesize=10000, product=False)
#                         st -= 1
#                         for r in reactions:
#                             if r.id in added_reactions:
#                                 continue
#                             added_reactions.add(r.id)
#                             for m in r.molecules:
#                                 obj_mol = m.molecule
#                                 structure = obj_mol.structure
#                                 if m.is_product:
#                                     tshki.add(structure)
#                                     triples.add((z_mol, structure, 1))
#                                     w.write(ReactionContainer(reactants=[z_mol], products=[structure], meta={'fer': 1,
#                                                                                                         'id': r.id}))
#                                     if st:
#                                         if obj_mol not in seen:
#                                             seen.add(obj_mol)
#                                             stack.append((obj_mol, st))
#                     print('done')
#                 r = len(triples)
#                 for m in all_mols[i+1:]:
#                     ind_tan = evaluation(z_mol, m)
#                     if .9 <= ind_tan <= 1:
#                         trees(m)
#                     elif .4 <= ind_tan <= .5:
#                         trees(m)
#                     elif 0 < ind_tan <= .1:
#                         trees(m)
#         print(i, 'finish')
#         if i == 0:
#             break


# with open('zinc.pickle', 'rb') as f:
#     zinc = pickle.load(f)
#     n = 0
#     with db_session:
#         for chunk in chunks(range(0, 1600000), size=1000):
#             with db_session:
#                 for x in chunk:
#                     try:
#                         m = Molecule[x]
#                         if bytes(m.structure) in zinc:
#                             if not n % 1000:
#                                 print(n, 'done')
#                             try:
#                                 MoleculeProperties(user=User[1], structure=m, data={'zinc': 1})
#                             except:
#                                 pass
#                             n += 1
#                     except:
#                         pass

# n = 0
# for signat in mols:
#     with db_session:
#         mol = list(MoleculeStructure.select(lambda x: x.signature == signat))
#         if mol:
#             if not n % 100:
#                 print(n, 'done')
#             m = mol[0].molecule
#             try:
#                 MoleculeProperties(user=User[1], structure=m, data={'zinc': 1})
#             except:
#                 pass
#             n += 1


# nu = 0
# with RDFread('../../../tmp/final3.rdf', indexable=True) as f:
#     f.seek(3453000)
#     with db_session:
#         for r in f:
#             if not nu % 1000:
#                 print(f'---{nu}-done---')
#             nu += 1
#             try:
#                 r.standardize()
#                 r.aromatize()
#                 year = r.meta['file_name']
#                 try:
#                     year = int(year[:4])
#                 except:
#                     try:
#                         year = int(year[6:10])
#                     except:
#                         try:
#                             year = int(year[1:5])
#                         except:
#                             pass
#                 reagentes = ', '.join([str(x) for x in r.reagents])
#                 r.reagents.clear()
#                 metaa = {'source_id': r.meta['source_id'], 'text': r.meta['text'], 'reagents': reagentes, 'year': year}
#                 if not Reaction.structure_exists(r):
#                     ri = Reaction(r, User[1])
#                     if ri:
#                         if not ReactionConditions.exists(structure=ri, data=metaa):
#                             ReactionConditions(user=User[1], structure=ri, data=metaa)
#             except Exception as e:
#                 warning(e)
# drugs_in_patents_as_product = 447
# zinc = set()
# for x in range(1, 1872):
#     try:
#         with SDFread(f'zinc/{x}.sdf', indexable=True) as f:
#             print(f'do zinc/{x}.sdf')
#             num = 0
#             l = len(f)
#             print('teor_mols', l)
#             for z_mol in f:
#                 num += 1
#                 z_mol.standardize()
#                 z_mol.aromatize()
#                 z_mol.implicify_hydrogens()
#                 zinc.add(bytes(z_mol))
#                 with db_session:
#                     if Molecule.structure_exists(z_mol):
#                         print('exists')
#                         w.write(z_mol)
#                         continue
#                     mols = Molecule.find_similar(z_mol, page=1, pagesize=1)
#                     if mols and mols[0][-1] >= .95:
#                         print('similar')
#                         mol = [mols[0][0].structure]
#                         r = ReactionContainer(reactants=[z_mol], products=mol)
#                         we.write(r)
#             print('fact_mols', num, '\n')
#     except:
#         pass
# with open('zinc.pickle', 'wb') as z:
#     pickle.dump(zinc, z)
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
#         r.standardize()
#         r.aromatize()
#         metaa = r.meta
#         year = metaa['file_name']
#         try:
#             year = int(year[:4])
#         except:
#             try:
#                 year = int(year[6:10])
#             except:
#                 try:
#                     year = int(year[1:5])
#                 except:
#                     pass
#         reagentes = ', '.join([str(x) for x in r.reagents])
#         r.reagents = []
#         with db_session:
#             ri = Reaction(r, User[1])
#             description = '; '.join([metaa['source_id'], f'reagents: {reagentes}', metaa['text'], f'year: {year}'])
#             if not Reaction.structure_exists(r):
#                 ReactionConditions(user=User[1], structure=ri, data={'data': description})
#             if not ReactionConditions.exists(structure=ri, data=description):
#                 ReactionConditions(user=User[1], structure=ri, data={'data': description})
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


# if __name__ == '__main__':
#     i = 0
#     with RDFread('../../../tmp/final3.rdf') as reactions:
#         # inp = Queue()
#         # for _ in range(12):
#         #     Process(target=worker, args=(inp,)).start()
#         for r in reactions:
#             r.standardize()
#             r.aromatize()
#             year = r.meta['file_name']
#             try:
#                 year = int(year[:4])
#             except:
#                 try:
#                     year = int(year[6:10])
#                 except:
#                     try:
#                         year = int(year[1:5])
#                     except:
#                         pass
#             reagentes = ', '.join([str(x) for x in r.reagents])
#             r.reagents.clear()
#             metaa = {'source_id': r.meta['source_id'], 'text': r.meta['text'], 'reagents': reagentes, 'year': year}
#             with db_session:
#                 ri = Reaction(r, User[1])
#                 if not Reaction.structure_exists(r):
#                     ReactionConditions(user=User[1], structure=ri, data=metaa)
#                 if not ReactionConditions.exists(structure=ri, data=metaa):
#                     ReactionConditions(user=User[1], structure=ri, data=metaa)
#             if not i % 1000:
#                 print(f'--------{i} done--------')
        #     inp.put(x)
        #     i += 1
        # for _ in range(12):
        #     inp.put('STOP')
# from psycopg2.extensions import register_adapter
# from psycopg2.extras import Json
# from time import sleep
# from joblib import Parallel, delayed
# register_adapter(dict, Json)
#
#
# def worker(q):
#     i = 3453000
#     with RDFread('../../../tmp/final3.rdf', indexable=True) as f:
#         f.seek(3453000)
#         for re in f:
#             if q.qsize() > 1000:
#                 sleep(65)
#             q.put(re)
#             if not i % 1000:
#                 print(f'--------{i} done--------')
#             i += 1
#         for _ in range(12):
#             q.put('done')
#         warning('worker done')


# def exists(reaction, db):
#     with db_session:
#         return db.Reaction.structure_exists(reaction)
#
#
# def insert(reaction, meta, db):
#     with db_session:
#         u = db.User[1]
#         if not db.Reaction.structure_exists(reaction):
#             ri = db.Reaction(reaction, u)
#             ReactionConditions(user=u, structure=ri, data=meta)
#
#
# def add_meta(reaction, meta, db):
#     with db_session:
#         u = db.User[1]
#         ri = db.Reaction.find_structure(reaction)
#         if ri:
#             if not db.ReactionConditions.exists(structure=ri, data=meta):
#                 ReactionConditions(user=u, structure=ri, data=meta)
#         else:
#             insert(reaction, meta, db)
#
#
# def populate_reactions(r):
#     db = load_schema('all_patents', user='postgres', password='jyvt0n3', host='localhost', database='postgres')
#     try:
#         r.standardize()
#         r.aromatize()
#         year = r.meta['file_name']
#         try:
#             year = int(year[:4])
#         except:
#             try:
#                 year = int(year[6:10])
#             except:
#                 try:
#                     year = int(year[1:5])
#                 except:
#                     pass
#         reagentes = ', '.join([str(x) for x in r.reagents])
#         r.reagents.clear()
#         metaa = {'source_id': r.meta['source_id'], 'text': r.meta['text'], 'reagents': reagentes, 'year': year}
#         if not exists(r, db):
#             try:
#                 insert(r, metaa, db)
#             except:
#                 try:
#                     insert(r, metaa, db)
#                 except:
#                     add_meta(r, metaa, db)
#         else:
#             add_meta(r, metaa, db)
#         del r, reagentes, year, metaa
#         db.disconnect()
#     except Exception as e:
#         db.disconnect()
#         warning(e)
#
#
# def main():
#     q = Queue()
#     Process(target=worker, args=[q]).start()
#     warning('started')
#     pr = [Process(target=populate_reactions, args=[q]) for _ in range(12)]
#     for p in pr:
#         p.start()


# f = RDFread('../../../tmp/final3.rdf', indexable=True)
# f.seek(3453000)
# Parallel(n_jobs=3, verbose=1000, pre_dispatch='1.5*n_jobs')(delayed(populate_reactions)(r) for r in f)
