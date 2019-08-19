from CGRdb import load_schema
from CGRdb.database import Molecule
from CGRdb.search.fingerprints import FingerprintMolecule
from CGRtools.files import RDFwrite
from CGRtools.containers import ReactionContainer
import numpy as np
from pony.orm import select, db_session
from pickle import dump, load
from random import choices

load_schema('all_patents', )


def evaluation(query, res):
    """
    оценка нод.
    возвращает танимото для пары запрос-результат.
    """
    query_fp, res_fp = (FingerprintMolecule.get_fingerprint(x) for x in [query, res])
    qc, rc, common = len(query_fp), len(res_fp), len(query_fp.intersection(res_fp))
    return common / (qc + rc - common)


with open('pairs_from_reactant_to_product.pickle', 'rb') as f:
    pairs = load(f)


sigma = list(pairs.keys())
print(len(sigma))
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


def writer(t1, t2, tto):
    unic_tshki_1 = t1 - t2
    print(len(unic_tshki_1))
    unic_tshki_2 = t2 - t1
    print(len(unic_tshki_2), '\n')
    tshki = choices(list(unic_tshki_1), k=2)
    for t_1 in tshki:
        mol_t_1 = Molecule[t_1].structure
        ind = 0
        for t_2 in unic_tshki_2:
            mol_t_2 = Molecule[t_2].structure
            new_ind = evaluation(mol_t_1, mol_t_2)
            if new_ind > ind:
                ind = new_ind
                tsh_2 = t_2
        fw.write(ReactionContainer(reactants=[Molecule[mol_id_1].structure],
                                   products=[Molecule[tsh_2].structure], meta={'fer': tto / 2}))
        print('written')


with RDFwrite('False_pairs_.rdf') as fw:
    for i, mol_id_1 in enumerate(sigma):
        print(f'-- {i} powel --')
        with db_session:
            high = middle = low = 0
            high_mol_id = middle_mol_id = low_mol_id = 0
            for mol_id_2 in sigma:
                if mol_id_1 == mol_id_2:
                    continue
                ind_tan = big[old_to_new[mol_id_1]][old_to_new[mol_id_2]] / 255
                if .9 <= ind_tan <= 1:
                    if ind_tan > high:
                        high_mol_id = mol_id_2
                        high = ind_tan
                elif .4 <= ind_tan <= .5:
                    if ind_tan > middle:
                        middle_mol_id = mol_id_2
                        middle = ind_tan
                elif 0 < ind_tan <= .1:
                    if ind_tan > low:
                        low_mol_id = mol_id_2
                        low = ind_tan

            tshki_1 = pairs[mol_id_1]
            if high_mol_id:
                writer(tshki_1, pairs[high_mol_id], high)
            if middle_mol_id:
                writer(tshki_1, pairs[middle_mol_id], middle)
            if low_mol_id:
                writer(tshki_1, pairs[low_mol_id], low)
            print(f'-- finish {i} --')
