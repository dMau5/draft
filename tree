from pickle import load, dump
from pony.orm import db_session, select
from CGRdb import load_schema, Molecule, Reaction


n = 1
# ондостадийные реакции (ид реактанта ключ, значение сет тюплов ид продуктов и ид реакции, которые могут быть получены 
# из реактанта в данной реакции)
with open('new_from_prod_to_react.pickle', 'rb') as f:
    reac_prod = load(f)
# ид реакций, которые дала Адклия, с них я продолжаю цепочку, то есть продукты этих реакций являются началом цепочки до
# 4 стадий
with open('reactions.pickle', 'rb') as f:
    dt = load(f)
print(len(dt))
gl_seen = set()
for i, r in enumerate(dt):
    print(i, 'powel')
    with db_session:
        reaction = Reaction[r]
        if reaction:
            data = []
            flag = False
            for m in reaction._molecules:
                obj_mol = m.molecule
                mol = obj_mol.structure
                m_id = obj_mol.id
                # проверяю если это не продукт или атомов в молекуле меньше 2 или эта ид уже рассматривалась
                if not m.is_product or len(mol._atoms) < 2 or m_id in gl_seen:
                    continue
                print('mol', obj_mol)
                gl_seen.add(m_id)
                # бфс для молекулы
                if m_id in reac_prod:
                    targets = []
                    stages = 4
                    n += 1
                    stack = []
                    seen = set()
                    stages -= 1
                    poi = abs(stages - 5)
                    for (product, reaction_id) in reac_prod[m_id]:
                        if product not in seen:
                            stack.append((product, stages))
                        seen.add(product)

                    while stack:
                        mt, st = stack.pop(0)
                        st -= 1
                        poii = abs(st - 5)
                        # если молекулы нет в словаре или стадии исчерпаны, добавляю молекулу и стадию в лист
                        if mt in reac_prod:
                            for (duct, tion_id) in reac_prod[mt]:
                                if st:
                                    if duct not in seen:
                                        stack.append((duct, st))
                                else:
                                    targets.append((mt, 5))
                                seen.add(duct)
                        else:
                            targets.append((mt, poii - 1))
                    # сохраняю полученные молекулы
                    with open(f'{m_id}.pickle', 'wb') as w1:
                        dump(targets, w1)

    print('usego', n)
