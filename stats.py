from CGRtools.files import SDFread, SDFwrite, RDFread, RDFwrite
from CGRdb import Molecule, load_schema, Reaction
from networkx import DiGraph
from collections import defaultdict
from pony.orm import db_session
from pickle import load, dump
from watch import paths_of_synthesis_for_target_molecule
from watch import visualization

load_schema('all_patents', user='postgres', password='jyvt0n3', host='localhost', database='postgres')
with open('zinc.pickle', 'rb') as z:
    zinc = load(z)


drugs_in_reactions = []
with SDFread('drugs_in_patents_as_product.sdf') as d:
    # with SDFwrite('trusted_drugs.sdf') as w:
        with db_session:
            for drug in d:
                drug.aromatize()
                drug.standardize()
                drug.implicify_hydrogens()
                drug_in_db = Molecule.find_structure(drug)
                seen = {drug_in_db}
                stack = [(drug_in_db, 32)]
                added_reactions = set()
                g = DiGraph()
                n = 0
                paths = 0
                flag = False
                subgr = [drug]
                while stack:
                    mt, st = stack.pop(0)
                    st -= 1
                    reactions = mt.reactions_entities(pagesize=100, product=True)
                    for r in reactions:
                        if r.id in added_reactions:
                            continue
                        data = [x.data['source_id'] for x in list(r.metadata)]
                        g.add_node(n, data='; '.join(data))
                        g.add_edge(n, mt.structure)
                        added_reactions.add(r.id)
                        components = r.molecules
                        reactants = [x.molecule for x in components if not x.is_product]
                        if all(bytes(x.structure) in zinc for x in reactants):
                            subgr.append(n)
                            flag = True
                        for m in components:
                            obj_mol = m.molecule
                            structure = obj_mol.structure
                            if flag:
                                subgr.append(structure)
                            if obj_mol in reactants:
                                if st and obj_mol not in seen:
                                    seen.add(obj_mol)
                                    if bytes(mt.structure) not in zinc:
                                        stack.append((obj_mol, st))
                                g.add_edge(structure, n)
                            else:
                                g.add_edge(n, structure)
                        n += 1
                sub = g.subgraph(subgr)


                # while stack:
                #     mt, st = stack.pop(0)
                #     st -= 1
                #     reactions = mt.reactions_entities(pagesize=20, product=True)
                #     for r in reactions:
                #         if r.id in added_reactions:
                #             continue
                #         added_reactions.add(r.id)
                #         mooleecuulees = [x.molecule for x in r.molecules]
                #         reactants = [x.structure for x in mooleecuulees if not x.is_product]
                #         if all(bytes(x) in zinc for x in reactants):
                #             data = [x.data['source_id'] for x in list(r.metadata)]
                #             g.add_node(n, data=data)
                #             g.add_edge(n, mt.structure)
                #             for m in r.molecules:
                #                 obj_mol = m.molecule
                #                 structure = obj_mol.structure
                #                 if not m.is_product:
                #                     g.add_edge(structure, n)
                #                 else:
                #                     g.add_edge(n, structure)
                #         else:
                #             for m2 in r.molecules:
                #                 if not m2.is_product:
                #                     obj_mol = m2.molecule
                #                     if st and obj_mol not in seen:
                #                         seen.add(obj_mol)
                #                         stack.append((obj_mol, st))
                #         if st == 49:
                #             paths += 1
                #         n += 1
                #     w.write(m)
                #     t = (m, 50 - st, paths, g)
                #     drugs_in_reactions.append(t)
# with open('all_stars.pickle', 'wb') as c:
#     dump(drugs_in_reactions, c)

