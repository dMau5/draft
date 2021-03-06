from CGRdb import load_schema
from CGRdb.database import Molecule
from CGRdb.search.fingerprints import FingerprintMolecule
from CGRtools.files import RDFwrite
from CGRtools.containers import ReactionContainer
import numpy as np
from math import ceil
from multiprocessing import Queue, Process
from pony.orm import select, db_session
from pickle import dump, load
from itertools import islice, cycle, filterfalse
from random import randint
from time import time

load_schema('all_patents',)

with open('pairs_from_reactant_to_product.pickle', 'rb') as f:
    pairs = load(f)

trues = []
rank = 1000
n = 1
file = 10154
with open('True_pairs2.txt', 'w') as w:
    for reactant, values in pairs.items():
        with db_session:
            reactant = Molecule[reactant].structure
            print('tts -->', len(values))
            for product in Molecule.select(lambda x: x.id in values.keys()):
                stages = values[product.id]
                product = product.structure
                trues.append((reactant, product, True, stages))
                rank -= 1
                if not rank:
                    with open(f'True_pairs/{file}.pickle', 'wb') as ww:
                        dump(trues, ww)
                    rank = 1000
                    file += 1
                    trues = []
                w.write(str(reactant) + ', ')
                w.write(str(product) + ', ')
                w.write('True, ')
                w.write(str(stages) + '\n')
        print('done -->', n, ' file -->', file, ' rank -->', rank)
        n += 1


def evaluation(query, res):
    """
    оценка нод.
    возвращает танимото для пары запрос-результат.
    """
    # query_fp, res_fp = (FingerprintMolecule.get_fingerprint(x) for x in [query, res])
    qc, rc, common = len(query), len(res), len(query.intersection(res))
    return common / (qc + rc - common)


# with open('pairs_from_reactant_to_product.pickle', 'rb') as f:
#     pairs = load(f)


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

# with open('fingerprints.pickle', 'rb') as f0:
#     fps = load(f0)

# with open('old_to_new_ids.pickle', 'rb') as f1:
#     old_to_new = load(f1)
#
# with open('array_of_similarity_A.pickle', 'rb') as f2:
#     big = load(f2)
#
# with open('list_of_A.pickle', 'rb') as we:
#     sigma = load(we)
#
# with open('naoborot.pickle', 'rb') as f3:
#     new_dict = load(f3)


# def roundrobin1(*iterables):
#     num_active = len(iterables)
#     iterables = cycle(iterables)
#     while num_active:
#         try:
#             for nxt in iterables:
#                 yield next(nxt)
#         except StopIteration:
#             num_active -= 1
#             iterables = cycle(islice(iterables, num_active))
#
#
# def roundrobin2(iterables):
#     dt = []
#     for tu_sh, tu, A_sh in iterables:
#         try:
#             yield next(tu_sh), tu, A_sh
#             dt.append((tu_sh, tu, A_sh))
#         except StopIteration:
#             pass
#
#     num_active = len(dt)
#     iterables = cycle(dt)
#     while num_active:
#         try:
#             for tu_sh, tu, A_sh in iterables:
#                 yield next(tu_sh), tu, A_sh
#         except StopIteration:
#             num_active -= 1
#             iterables = cycle(islice(iterables, num_active))
#
#
# def unique_everseen(iterable, key=None):
#     "List unique elements, preserving order. Remember all elements ever seen."
#     # unique_everseen('AAAABBBCCDAABBB') --> A B C D
#     # unique_everseen('ABBCcAD', str.lower) --> A B C D
#     seen = set()
#     seen_add = seen.add
#     if key is None:
#         for element in filterfalse(seen.__contains__, iterable):
#             seen_add(element)
#             yield element
#     else:
#         for element in iterable:
#             k = key(element)
#             if k not in seen:
#                 seen_add(k)
#                 yield element
#
#
# def worker(p, o):
#     for m_A, m_ash, tsh in iter(p.get, 'STOP'):
#         fp_1 = fps[tsh]
#         tu = pairs[m_A] - pairs[m_ash]
#         # fp_1, molecule_1 = m_1['fp'], m_1['molecule']
#         index = 0
#         for ii in tu:
#             # with open(f'database/{ii}.pickle', 'rb') as p2:
#             #     fp_2 = load(p2)
#             #     fp_2 = m_2['fp']
#             fp_2 = fps[ii]
#             new_ind = evaluation(fp_1, fp_2)
#             if new_ind > index:
#                 index = new_ind
#         # with db_session:
#         #     o.put(ReactionContainer(reactants=[m_A],
#         #                             products=[Molecule[tsh].structure], meta={'fer': index / 2}))
#         o.put((m_A, tsh, index))
#
#
# if __name__ == '__main__':
#     with RDFwrite('False_pairs_.rdf') as fw:
#         inp = Queue()
#         out = Queue()
#         for _ in range(20):
#             Process(target=worker, args=(inp, out, )).start()
#
#         n = 23958
#         rank = 1000
#         respa = []
#         # holost = 15700
#         for i, mol_id_1 in enumerate(sigma[6905:], 6905):
#             # with open(f'database_1/{mol_id_1}.pickle', 'rb') as pp:
#             #     m = load(pp)
#             #     molecule = m['molecule']
#             print(f'-- {i} powel --')
#             # if i == 3:
#             #     exit()
#             new_id = old_to_new[mol_id_1]
#
#             tshki_1 = pairs[mol_id_1]
#             total = len(tshki_1)
#             print('total', total)
#             # ini = time()
#             t_sh01 = ((iter(pairs[new_dict[n]] - tshki_1), tshki_1 - pairs[new_dict[n]], new_dict[n])
#                       for n, x in enumerate(big[new_id]) if x <= 25.5 and new_id != n)
#             t_sh45 = ((iter(pairs[new_dict[n]] - tshki_1), tshki_1 - pairs[new_dict[n]], new_dict[n])
#                       for n, x in enumerate(big[new_id]) if 102 <= x <= 127.5 and new_id != n)
#             t_sh90 = ((iter(pairs[new_dict[n]] - tshki_1), tshki_1 - pairs[new_dict[n]], new_dict[n])
#                       for n, x in enumerate(big[new_id]) if x >= 229.5 and new_id != n)
#
#             for triple in islice(unique_everseen(roundrobin2(roundrobin1(t_sh01, t_sh45, t_sh90)), lambda x: x[0]),
#                                  total):
#                 t_sh, t, m_Ash = triple
#                 inp.put((mol_id_1, m_Ash, t_sh))
#             for _ in range(total):
#                 print(f'for {i} ||', 'inp-->', inp.qsize(), 'out-->', out.qsize(), 'rank-->', rank, 'file-->', n)
#                 respa.append(out.get())
#                 # holost -= 1
#                 # if holost <= 0:
#                 rank -= 1
#                 if not rank:
#                     with open(f'triples/{n}.pickle', 'wb') as wq:
#                         dump(respa, wq)
#                     respa = []
#                     n += 1
#                     rank = 1000
#             # print('prowlo', time() - ini, 'sec')
#
#         for _ in range(20):
#             inp.put('STOP')
