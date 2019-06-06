from CGRtools.files import SDFread, SDFwrite, RDFread, RDFwrite
from CGRdb import Molecule, load_schema, Reaction
from networkx import DiGraph
from collections import defaultdict
from pony.orm import db_session
from pickle import load, dump
from watch import visualization

load_schema('all_patents', )
with open('zinc.pickle', 'rb') as z:
    zinc = load(z)

ter = 0
drugs_in_reactions = []
with SDFread('drugs_in_patents_as_product.sdf') as d:
    with SDFwrite('trusted_drugs.sdf') as w:
        with db_session:
            for drug in d:
                drug.aromatize()
                drug.standardize()
                drug.implicify_hydrogens()
                print('start ==', drug)
                drug_in_db = Molecule.find_structure(drug)
                print('create graph')
                seen = {drug_in_db}
                stages = 32
                stack = [(drug_in_db, stages)]
                added_reactions = set()
                g = DiGraph()
                n = 0
                while stack:
                    mt, st = stack.pop(0)
                    st -= 1
                    reactions = mt.reactions_entities(pagesize=100, product=True)
                    for r in reactions:
                        if r.id in added_reactions:
                            continue
                        data = list(r.metadata)[0].data
                        g.add_node(n, data=f"{data['source_id']}, {data['text']}".replace('+\n', ''))
                        added_reactions.add(r.id)
                        for m in r.molecules:
                            obj_mol = m.molecule
                            structure = obj_mol.structure
                            atoms = [atom[-1].element for atom in structure.atoms()]
                            if not m.is_product:
                                if st and obj_mol not in seen:
                                    seen.add(obj_mol)
                                    if mt == drug_in_db:
                                        stack.append((obj_mol, st))
                                elif bytes(structure) not in zinc and bytes(structure) not in zinc \
                                        and not (atoms.count('C') <= 2 or len(atoms) <= 6):
                                        stack.append((obj_mol, st))
                                g.add_edge(structure, n)
                            else:
                                g.add_edge(n, structure)
                            if bytes(structure) in zinc or atoms.count('C') <= 2 or len(atoms) <= 6:
                                g.nodes[structure]['zinc'] = 1
                        n += 1
                print('done')

                try:
                    synky = list(g._succ[drug])
                    if synky:
                        g.remove_nodes_from(synky)

                    st = []
                    for x in g.nodes():
                        if not isinstance(x, int) and x != drug:
                            atoms = [atom[-1].element for atom in x.atoms()]
                            if bytes(x) not in zinc and not (atoms.count('C') <= 2 or len(atoms) <= 6):
                                st.append(x)

                    while True:
                        for n, structure in enumerate(st):
                            if not g._pred[structure]:
                                st.pop(n)
                                synky = g._succ[structure]
                                g.remove_node(structure)
                                g.remove_nodes_from(synky)
                                break
                        else:
                            break

                    seen = {drug}
                    queue = [drug]
                    while queue:
                        start = queue.pop(0)
                        for i in g._pred[start].keys() - seen:
                            queue.append(i)
                            seen.add(i)
                        if isinstance(start, int):
                            for i in g._succ[start]:
                                seen.add(i)

                    if len(seen) > 1:
                        g = g.subgraph(seen)
                        g = DiGraph(g)
                        w.write(drug)
                        t = (drug_in_db.id, g, drug.meta)
                        drugs_in_reactions.append(t)
                        # t = visualization(g, drug)
                        # r = 6
                    else:
                        print('not found\n')
                except:
                    print(Exception(KeyError), drug_in_db.id)

                ter += 1
                print(ter)


with open('stats.pickle', 'wb') as p:
    dump(drugs_in_reactions, p)

    # with db_session:
    #     for x in drugs_in_reactions:
    #         idd, g, meta = x
    #         molecule = Molecule[idd]
    #         t = visualization(g, molecule.structure)
    #         r = 6
    # seen = {drug_in_db}
    # stages = 32
    # stack = [(drug_in_db, stages, None)]
    # added_reactions = set()
    # g = DiGraph()
    # n = 0
    # paths = defaultdict(set)
    # while stack:
    #     mt, st, path = stack.pop(0)
    #     st -= 1
    #     if mt:
    #         reactions = mt.reactions_entities(pagesize=100, product=True)
    #         if bytes(mt.structure) not in zinc:
    #             if not st or not reactions or any(x.id in added_reactions for x in reactions):
    #                 common = set()
    #                 other_del = set()
    #                 for k, x in paths.items():
    #                     if k != path:
    #                         common.update(x)
    #                     if mt.structure in x:
    #                         other_del.update(x)
    #                 deel = paths[path]
    #                 deel.update(other_del)
    #                 deel = deel.difference(common)
    #                 g.remove_nodes_from(deel)
    #                 stack = [x for x in stack if x[2] != path]
    #                 continue
    #         for r in reactions:
    #             if r.id in added_reactions:
    #                 continue
    #             data = list(r.metadata)[0].data
    #             g.add_node(n, data=f"{data['source_id']}, {data['text']}".replace('+\n', ''))
    #             g.add_edge(n, mt.structure)
    #             added_reactions.add(r.id)
    #             for m in r.molecules:
    #                 obj_mol = m.molecule
    #                 structure = obj_mol.structure
    #                 if not m.is_product:
    #                     if st and obj_mol not in seen:
    #                         seen.add(obj_mol)
    #                         if mt == drug_in_db:
    #                             stack.append((obj_mol, st, n))
    #                         elif bytes(mt.structure) not in zinc and bytes(structure) not in zinc:
    #                             stack.append((obj_mol, st, path))
    #                     g.add_edge(structure, n)
    #                 else:
    #                     g.add_edge(n, structure)
    #                 if mt == drug_in_db:
    #                     paths[n].add(structure)
    #                     paths[n].add(n)
    #                 else:
    #                     paths[path].add(structure)
    #                     paths[path].add(n)
    #                 if bytes(structure) in zinc:
    #                     g.nodes[structure]['zinc'] = 1
    #             n += 1
    # ter += 1
    # print(ter)
    # if len(g) > 1 and g.edges():
    #     t = (drug_in_db.id, g, drug.meta)
    #     y = visualization(g, drug)
    #     drugs_in_reactions.append(t)
    #     w.write(drug)
    # else:
    #     print(drug, 'paths not found')
    #

    # with db_session:
    #     for x in drugs_in_reactions:
    #         idd, g, meta = x
    #         drug = Molecule[idd]
    #         t = visualization(g, drug.structure)
    #         b = 6
    #
    #
    #             while stack:
    #                 mt, st = stack.pop(0)
    #                 st -= 1
    #                 reactions = mt.reactions_entities(pagesize=20, product=True)
    #                 for r in reactions:
    #                     if r.id in added_reactions:
    #                         continue
    #                     added_reactions.add(r.id)
    #                     mooleecuulees = [x.molecule for x in r.molecules]
    #                     reactants = [x.structure for x in mooleecuulees if not x.is_product]
    #                     if all(bytes(x) in zinc for x in reactants):
    #                         data = [x.data['source_id'] for x in list(r.metadata)]
    #                         g.add_node(n, data=data)
    #                         g.add_edge(n, mt.structure)
    #                         for m in r.molecules:
    #                             obj_mol = m.molecule
    #                             structure = obj_mol.structure
    #                             if not m.is_product:
    #                                 g.add_edge(structure, n)
    #                             else:
    #                                 g.add_edge(n, structure)
    #                     else:
    #                         for m2 in r.molecules:
    #                             if not m2.is_product:
    #                                 obj_mol = m2.molecule
    #                                 if st and obj_mol not in seen:
    #                                     seen.add(obj_mol)
    #                                     stack.append((obj_mol, st))
    #                     if st == 49:
    #                         paths += 1
    #                     n += 1
    #                 w.write(m)
    #                 t = (m, 50 - st, paths, g)
    #                 drugs_in_reactions.append(t)
# with open('all_stars.pickle', 'wb') as c:
#     dump(drugs_in_reactions, c)

