from CGRdb import load_schema
from CGRdb.database import Molecule
from CGRdb.search.fingerprints import FingerprintMolecule
from CGRtools.files import RDFwrite
from CGRtools.containers import ReactionContainer
import numpy as np
from math import ceil
from pony.orm import select, db_session
from pickle import dump, load
from itertools import islice, cycle, filterfalse
from random import randint
from time import time

load_schema('all_patents', )


def evaluation(query, res):
    """
    оценка нод.
    возвращает танимото для пары запрос-результат.
    """
    # query_fp, res_fp = (FingerprintMolecule.get_fingerprint(x) for x in [query, res])
    qc, rc, common = len(query), len(res), len(query.intersection(res))
    return common / (qc + rc - common)


with open('pairs_from_reactant_to_product.pickle', 'rb') as f:
    pairs = load(f)


# big = np.eye(len(sigma), dtype=np.uint8)
# old_to_new = {sigma[0]: 0}
# n = 0
# for i, id_1 in enumerate(sigma):
#     print(f'--- {i} powel ---')
#     with open(f'database/{id_1}.pickle', 'rb') as f1:
#         fp_1 = load(f1)
#     for id_2 in sigma[i+1:]:
#         with open(f'database/{id_2}.pickle', 'rb') as f2:
#             fp_2 = load(f2)
#         if not i:
#             n += 1
#             old_to_new[id_2] = n
#         if not n % 100:
#             print(f'--- uwlo {n} ---')
#         ind = evaluation(fp_1, fp_2)
#         new_id_2 = old_to_new[id_2]
#         big[i][new_id_2] = big[new_id_2][i] = ind * 255

with open('old_to_new_ids.pickle', 'rb') as f1:
    old_to_new = load(f1)
with open('array_of_similarity_A.pickle', 'rb') as f2:
    big = load(f2)
with open('list_of_A.pickle', 'rb') as we:
    sigma = load(we)
    print(len(sigma))
with open('naoborot.pickle', 'rb') as f3:
    new_dict = load(f3)


# def writer(t1, t2, num):
#     unic_tshki_1 = t1 - t2
#     print('started mols', len(unic_tshki_1), 'to', num)
#     unic_tshki_2 = t2 - t1
#     print('u_tshki_2--->', len(unic_tshki_2))
#     for t_1 in list(unic_tshki_1)[: num + 1]:
#         with open(f'database/{t_1}.pickle', 'rb') as p1:
#             fp_1 = load(p1)
#         ind = 0
#         for t_2 in unic_tshki_2:
#             with open(f'database/{t_2}.pickle', 'rb') as p2:
#                 fp_2 = load(p2)
#             new_ind = evaluation(fp_1, fp_2)
#             if new_ind > ind:
#                 ind = new_ind
#                 tsh_2 = t_2
#         if ind:
#             fw.write(ReactionContainer(reactants=[Molecule[mol_id_1].structure],
#                                        products=[Molecule[tsh_2].structure], meta={'fer': ind / 2}))
#
#         upd = 0
#         diff = len(unic_tshki_1) - num
#         if diff < 0:
#             upd = - diff
#         return upd


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
    for tu_sh, tu in iterables:
        print(tu_sh, tu)
        try:
            tu_sh = next(tu_sh)
            dt.append((tu_sh, tu))
            # print(tu_sh, tu)
            yield tu_sh, tu
        except StopIteration:
            pass
    exit()
    num_active = len(dt)
    iterables = cycle(dt)
    while num_active:
        try:
            for tu_sh, tu in iterables:
                yield next(tu_sh), tu
        except StopIteration:
            num_active -= 1
            dt.pop(0)
            iterables = cycle(dt)


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


with RDFwrite('False_pairs_.rdf') as fw:
    for i, mol_id_1 in enumerate(sigma):
        print(f'-- {i} powel --')
        new_id = old_to_new[mol_id_1]
        with db_session:
            tshki_1 = pairs[mol_id_1]
            total = len(tshki_1)
            print('total', total)
            t_sh01 = ((iter(pairs[new_dict[n]] - tshki_1), tshki_1 - pairs[new_dict[n]])
                      for n, x in enumerate(big[new_id]) if x <= 25.5 and new_id != n)
            t_sh45 = ((iter(pairs[new_dict[n]] - tshki_1), tshki_1 - pairs[new_dict[n]])
                      for n, x in enumerate(big[new_id]) if 102 <= x <= 127.5 and new_id != n)
            t_sh90 = ((iter(pairs[new_dict[n]] - tshki_1), tshki_1 - pairs[new_dict[n]])
                      for n, x in enumerate(big[new_id]) if x >= 229.5 and new_id != n)
            # print('ini', time())
            # rr1 = roundrobin1(t_sh01, t_sh45, t_sh90)
            # print('rr1', time())
            # rr2 = roundrobin2(rr1)
            # print('rr2', time())
            # uee = unique_everseen(rr2, lambda x: x[0])
            # print('uee', time())
            # exit()
            n = 0
            for ts, t in unique_everseen(roundrobin2(roundrobin1(t_sh01, t_sh45, t_sh90)), lambda x: x[0]):
                n += 1
                with open(f'database/{ts}.pickle', 'rb') as p1:
                    fp_1 = load(p1)
                ind = 0
                for ii in t:
                    with open(f'database/{ii}.pickle', 'rb') as p2:
                        fp_2 = load(p2)
                    new_ind = evaluation(fp_1, fp_2)
                    if new_ind > ind:
                        ind = new_ind
                if ind:
                    fw.write(ReactionContainer(reactants=[Molecule[mol_id_1].structure],
                                               products=[Molecule[ts].structure], meta={'fer': ind / 2}))
            print('usego', n)

            # try:
            #     tshki_2_h = pairs[high_mol_id]
            #     tshki_2_m = pairs[middle_mol_id]
            #     tshki_2_l = pairs[low_mol_id]
            #     hi = len(tshki_2_h)
            #     mi = len(tshki_2_m)
            #     lo = len(tshki_2_l)
            #     to = hi + mi + lo
            #     hi = hi / to
            #     mi = mi / to
            #     lo = lo / to
            # except Exception as e:
            #     print(e)
            # balance = 0
            # if high_mol_id:
            #     print('high', 100 * hi)
            #     balance = writer(tshki_1, tshki_2_h, ceil(total * hi))
            # if middle_mol_id:
            #     print('middle', 100 * mi)
            #     if balance:
            #         balance = writer(tshki_1, tshki_2_m, ceil(total * mi) + balance)
            #     else:
            #         balance = writer(tshki_1, tshki_2_m, ceil(total * mi))
            # if low_mol_id:
            #     print('low', lo * 100)
            #     if balance:
            #         writer(tshki_1, tshki_2_l, ceil(total * lo) + balance)
            #     else:
            #         writer(tshki_1, tshki_2_l, ceil(total * lo))
            # print(f'-- finish {i} --')
