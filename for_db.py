from CGRdb import load_schema
from CGRdb.database import Molecule
from CGRdb.search.fingerprints import FingerprintMolecule
from CGRtools.files import RDFwrite
from CGRtools.containers import ReactionContainer
from pony.orm import select, db_session
from pickle import dump, load

load_schema('all_patents',)


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
with RDFwrite('False_pairs_.rdf') as fw:
    for i, mol_id_1 in enumerate(sigma):
        print(f'-- {i} powel --')
        with db_session:
            for mol_id_2 in sigma:
                if mol_id_1 == mol_id_2:
                    continue
                molecule_1 = Molecule[mol_id_1].structure
                molecule_2 = Molecule[mol_id_2].structure
                ind_tan = evaluation(molecule_1, molecule_2)
                if .9 <= ind_tan <= 1 or .4 <= ind_tan <= .5 or 0 < ind_tan <= .1:
                    tshki_1, tshki_2 = pairs[mol_id_1], pairs[mol_id_2]
                    unic_tshki_1 = tshki_1 - tshki_2
                    unic_tshki_2 = tshki_2 - tshki_1
                    for t_1 in unic_tshki_1:
                        mol_t_1 = Molecule[t_1].structure
                        ind = 0
                        for t_2 in unic_tshki_2:
                            mol_t_2 = Molecule[t_2].structure
                            new_ind = evaluation(mol_t_1, mol_t_2)
                            if new_ind > ind:
                                ind = new_ind
                                tsh_2 = t_2
                        fw.write(ReactionContainer(reactants=[molecule_1],
                                                   products=[Molecule[tsh_2].structure], meta={'fer': ind / 2}))
            print(f'-- finish {i} --')
