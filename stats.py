from CGRtools.files import SDFread, SDFwrite
from CGRdb import Molecule, load_schema
from collections import defaultdict
from pony.orm import db_session
from pickle import load, dump

with open('zinc.pickle', 'rb') as z:
    zinc = load(z)

with SDFwrite('er_aromatize.sdf') as w:
    with SDFread('drugs500mw.sdf') as f:
        with db_session:
            for m in f:
                m.standardize()
                try:
                    m.aromatize()
                except:
                    w.write(m)
                molecule = Molecule.find_structure(m)
                stats = defaultdict(dict)
                if molecule:
                    reactions_pr = molecule.reactions_entities(pagesize=100, product=True)
                    if reactions_pr:
                        stats[molecule.structure]['is_product'] = [x.id for x in reactions_pr]
                    reactions_re = molecule.reactions_entities(pagesize=100)
                    if reactions_re:
                        stats[molecule.structure]['is_reactant'] = [x.id for x in reactions_re]
with open('stats.pickle', 'wb') as p:
    dump(stats, p)

            # patent = dict()
            # added_reactions = set()
            # seen = set()
            # stack = [(molecule, 100)]
            # while stack:
            #     mt, st = stack.pop(0)
            #     st -= 1
            #     reactions = mt.reactions_entities(pagesize=100, product=True)
            #     for r in reactions:
            #         if r.id not in added_reactions:
            #             data = list(r.metadata)[0].data
            #             year = [x for x in data['text'].split(' ') if len(x) == 4]
                        # for y in year:
                        #     try:
                        #         patent[r.id] = int(y)
                        #     except:
                        #         pass
                        # if st == 99:
                        #     single[mt.id].add(r.id)
                        # elif st == 98:
                        #     double[mt.id].add(r.id)
                        # elif st == 97:
                        #     triple[mt.id].add(r.id)
                        # g.add_node(n, data=f"{data['source_id']}, {data['text']}".replace('+\n', ''))
                        # g.add_edge(n, mt.structure)
                        # added_reactions.add(r.id)
                    # else:
                    #     continue
                    # for m in r.molecules:
                    #     obj_mol = m.molecule
                    #     structure = obj_mol.structure
                    #     if not m.is_product:
                    #         if st and obj_mol not in seen:
                    #             seen.add(obj_mol)
                    #             stack.append((obj_mol, st))
                            # g.add_edge(structure, n)
                        # else:
                            # g.add_edge(n, structure)
                        # if bytes(structure) in zinc:
                        #     g.nodes[structure]['zinc'] = 1
                    # n += 1
