import pickle
from pony.orm import db_session
from CGRtools.containers import ReactionContainer
from CGRtools.files import SDFread, RDFread, SDFwrite, RDFwrite
from CGRtools.containers import MoleculeContainer
import networkx as nx
from multiprocessing import Process, Queue
from CGRdbUser import User
from CGRdb import load_schema, Molecule, Reaction
from CGRdb.database import MoleculeStructure
from CGRdbData import ReactionConditions, MoleculeProperties
from logging import warning, basicConfig, ERROR
from itertools import islice, chain
from time import sleep
from networkx import DiGraph

load_schema('bb')


def chunks(iterable, size=10):
    iterator = iter(iterable)
    for first in iterator:
        yield chain([first], islice(iterator, size - 1))


# with db_session:
#     all = set()
#     for molecule in Molecule.select():
#         m = molecule.structure
#         m.standardize()
#         m.implicify_hydrogens()
#         all.add(m)
# print(len(all))
# with open('building_blocks.pickle', 'wb') as w:
#     pickle.dump(all, w)
#
# pairs = []
# for z_mol in all:
#     mol_db = Molecule.find_structure(z_mol)
#     if mol_db:
#         print('create graph')
#         seen = {mol_db}
#         stages = 10
#         stack = [(mol_db, stages)]
#         added_reactions = set()
#         while stack:
#             mt, st = stack.pop(0)
#             reactions = mt.reactions_entities(pagesize=10000, product=False)
#             st -= 1
#             for r in reactions:
#                 if r.id in added_reactions:
#                     continue
#                 added_reactions.add(r.id)
#                 for m in r.molecules:
#                     obj_mol = m.molecule
#                     structure = obj_mol.structure
#                     if m.is_product:
#                         pairs.append((z_mol, structure))
#                         if st:
#                             if obj_mol not in seen:
#                                 seen.add(obj_mol)
#                                 stack.append((obj_mol, st))
#         print('done')

n = 0
for x in range(47):
    with SDFread(f'zinc/{x}.sdf') as f:
        with db_session:
            for m in f:
                try:
                    m.aromatize()
                except Exception as e:
                    print(e)
                    print(m)
                m.standardize()
                m.implicify_hydrogens()
                m_db = Molecule.find_structure(m)
                if m_db:
                    if not n % 1000:
                        print(n, 'done')
                    n += 1
                    try:
                        MoleculeProperties(user=User[1], structure=m_db, data={'zinc': 1})
                    except:
                        pass

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
# n = 0
# with SDFread('approved_drugs.sdf') as f:
#     with SDFread('drugs_in_patents_as_product.sdf') as w:
#         with db_session:
#             for drug in f:
#                 drug.standardize()
#                 drug.implicify_hydrogens()
#                 try:
#                     drug.aromatize()
#                 except:
#                     print(drug.meta)
#                 dr_db = Molecule.find_structure(drug)
#                 if dr_db:
#                     reactions = dr_db.reactions_entities(pagesize=100, product=True)
#                     if reactions:
#                         print(reactions)
#                         sleep(7)
#                         w.write(drug)
#                         n += 1
#                         print(n)
# basicConfig(level=ERROR)
# cycle_ring_error = 2669
# ideal_drugs = {173: 'antidepressant',
#                190: 'Antifibrinolytic hemostatic used in severe hemorrhage.',
#                274: 'A dopamine D2-receptor antagonist',
#                317: 'Allopurinol is a xanthine oxidase enzyme inhibitor that is considered to be one of the most '
#                     'effective drugs used to decrease urate levels and is frequently used in the treatment of chronic '
#                     'gout 6',
#                337: 'is a drug used to treat hypertension. Prazosin is marketed by Pfizer and was initially approved '
#                     'by the FDA in 1988 16. It belongs to the class of drugs known as alpha-1 antagonists Label, 8.',
#                379: 'A non-steroidal anti-inflammatory agent (anti-inflammatory agents, NON-steroidal) similar in mode '
#                     'of action to indomethacin.',
#                394: 'Trandolapril may be used to treat mild to moderate hypertension, to improve survival following '
#                     'myocardial infarction in clinically stable patients with left ventricular dysfunction, as an '
#                     'adjunct treatment for congestive heart failure, and to slow the rate of progression of renal '
#                     'disease in hypertensive individuals with diabetes mellitus and microalbuminuria or overt nephropathy.',
#                1117: 'Bevantolol is a beta-1 adrenoceptor antagonist that has been shown to be as effective as other beta '
#                      'blockers for the treatment of angina pectoris and hypertension.'}
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
