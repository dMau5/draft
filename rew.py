import pickle
from pony.orm import db_session
from CGRtools.containers import ReactionContainer
from CGRtools.files import SDFread, RDFread, SDFwrite, RDFwrite
from CGRtools.containers import MoleculeContainer
import networkx as nx
from multiprocessing import Process, Queue
from CGRdbUser import User
from CGRdb import load_schema, Molecule, Reaction
from CGRdbData import ReactionConditions


def worker(input_queue):
    for r in iter(input_queue.get, 'STOP'):
        r.standardize()
        r.aromatize()
        with db_session:
            if not Reaction.structure_exists(r):
                metaa = r.meta
                year = metaa['file_name']
                reagentes = ', '.join([str(x) for x in r.reagents])
                r.reagents = []
                ri = Reaction(r, User[1])
                description = '; '.join([metaa['source_id'],  f'reagents: {reagentes}', metaa['text'], f'year: {year}'])
                if not ReactionConditions.exists(structure=ri, data=description):
                    ReactionConditions(user=User[1], structure=ri, data={'data': description})


if __name__ == '__main__':
    i = 0
    with db_session:
        db = load_schema('sandbox')
        with RDFread('file_name.rdf') as reactions:
            inp = Queue()
            for _ in range(12):
                Process(target=worker, args=(inp,)).start()
            for x in reactions:
                if not i % 1000:
                    print(f'--------{i} done--------')
                inp.put(x)
                i += 1
            for _ in range(12):
                inp.put('STOP')


# def worker(q):
#     for r in RDFread('data/two_rules.rdf'):
#         q.put(r)
#     for _ in range(10):
#         q.put('done')
#
#
# def exists(reaction):
#     with db_session:
#         return db.Reaction.structure_exists(reaction)
#
#
# def insert(reaction, meta):
#     with db_session:
#         u = db.User[1]
#         if not db.Reaction.structure_exists(reaction):
#             ri = db.Reaction(reaction, u)
#             ReactionConditions(user=u, structure=ri, data=meta)
#
#
# def add_meta(reaction, meta):
#     with db_session:
#         u = db.User[1]
#         ri = db.Reaction.find_structure(reaction)
#         if not db.ReactionConditions.exists(structure=ri, data=meta):
#             ReactionConditions(user=u, structure=ri, data=meta)
#
#
# def populate_reactions(q):
#     for reaction in iter(q.get, 'done'):
#         meta = reaction.meta
#         if not exists(reaction):
#             try:
#                 insert(reaction, meta)
#             except:
#                 try:
#                     insert(reaction, meta)
#                 except:
#                     add_meta(reaction, meta)
#         else:
#             add_meta(reaction, meta)
#
#
# def main():
#     q = Queue()
#     Process(target=worker, args=[q]).start()
#     warning('started')
#     pr = [Process(target=populate_reactions, args=[q]) for _ in range(10)]
#     for p in pr:
#         p.start()
