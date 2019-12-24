from multiprocessing import Pool, Queue, Process
from ordered_set import OrderedSet
from CGRtools.files import RDFRead, RDFWrite
from CGRtools.containers import ReactionContainer, MoleculeContainer
from logging import basicConfig, INFO, info


def group_ions(reaction: ReactionContainer):
    """
    Ungroup molecules recorded as ions, regroup ions. Returns a tuple with the corresponding ReactionContainer and
    return code as int (0 - nothing was changed, 1 - ions were regrouped, 2 - ions are unbalanced).
    :param reaction: current reaction
    :return: tuple[ReactionContainer, int]
    """
    reaction_parts = []
    return_codes = []
    test_index = 0
    for molecules in (reaction.reactants, reaction.reagents, reaction.products):
        divided_molecules = [x for m in molecules for x in m.split('.')]
        test_index += 1

        if len(divided_molecules) == 0:
            reaction_parts.append(())
            continue
        elif len(divided_molecules) == 1 and int(divided_molecules[0]) == 0:
            return_codes.append(0)
            reaction_parts.append(molecules)
            continue
        elif len(divided_molecules) == 1:
            return_codes.append(2)
            reaction_parts.append(molecules)
            continue

        new_molecules = []
        cations, anions = [], []
        total_charge = 0
        for molecule in divided_molecules:
            mol_charge = int(molecule)
            total_charge += mol_charge
            if mol_charge == 0:
                new_molecules.append(molecule)
            elif mol_charge > 0:
                cations.append((mol_charge, molecule))
            else:
                anions.append((mol_charge, molecule))

        if len(cations) == 0 and len(anions) == 0:
            return_codes.append(0)
            reaction_parts.append(tuple(new_molecules))
            continue
        elif total_charge != 0:
            return_codes.append(2)
            reaction_parts.append(tuple(divided_molecules))
            continue
        else:
            error = False
            for cation_charge, cation in cations:
                salt = MoleculeContainer()
                salt = salt.union(cation)
                total_charge = cation_charge
                number_of_iterations = 1000
                while total_charge != 0 and number_of_iterations > 0 and len(anions):
                    anion_charge, anion = anions.pop(0)
                    salt = salt.union(anion)
                    total_charge += anion_charge
                    number_of_iterations -= 1
                if number_of_iterations == 0:
                    raise ValueError
                if len(anions) == 0 and total_charge != 0:
                    error = True
                    break
                new_molecules.append(salt)
            if error:
                return_codes.append(2)
                reaction_parts.append(tuple(new_molecules))
            else:
                return_codes.append(1)
                reaction_parts.append(tuple(new_molecules))

    return ReactionContainer(reactants=reaction_parts[0], reagents=reaction_parts[1], products=reaction_parts[2],
                             meta=reaction.meta), max(return_codes)


def remove_unchanged_parts(reaction: ReactionContainer) -> ReactionContainer:
    """
    Ungroup molecules, remove unchanged parts from reactants and products.
    :param reaction: current reaction
    :return: ReactionContainer
    """
    reactants = [m for m in reaction.reactants]
    products = [m for m in reaction.products]
    new_reactants, new_reagents, new_products = reactants, [m for m in reaction.reagents], products
    for reactant in reactants:
        if reactant in products:
            new_reagents.append(reactant)
            new_reactants.remove(reactant)
            new_products.remove(reactant)
    return ReactionContainer(reactants=tuple(new_reactants), reagents=tuple(new_reagents),
                             products=tuple(new_products), meta=reaction.meta)


def standardize(reaction: ReactionContainer):
    try:
        reaction.standardize()
        reaction.thiele()
        reaction.kekule()
        reaction.implicify_hydrogens()
        reaction.fix_mapping()
    except:
        return

    mistakes = []
    for molecule in reaction.molecules():
        for n, atom in molecule.atoms():
            try:
                if atom.is_radical:
                    return
            except:
                return
            try:
                if atom.isotope:
                    return
            except:
                return
        try:
            mistakes.extend(molecule.check_valence())
        except:
            return
    if mistakes:
        reaction.meta['mistake'] = f'Valence mistake in {set(mistakes)}'
        return

    try:
        reaction, return_code = group_ions(reaction)
        if return_code == 2:
            return
    except:
        return

    try:
        reaction = remove_unchanged_parts(reaction)
        if not reaction.reactants:
            return
        if not reaction.products:
            return
    except:
        return

    return reaction


def worker(ii):
    for i in iter(ii.get, 'stop'):
        data = OrderedSet()
        with RDFRead(f'/media/dinar/ec5b12ac-2741-412d-96ce-16b7e50e8072/reaxys/rdf_one_step/{i}.rdf', ignore=True) \
                as file:
            for r in file:
                standardized_reaction = standardize(r)
                if standardized_reaction:
                    data.add(standardized_reaction)
        print("{0} reactions passed..".format(len(data)))
        with RDFWrite(f'/media/data/arkadii_standardize/reactions/{i}.rdf') as out:
            for rr in data:
                out.write(rr)
        print(f' ---- {i} done')


if __name__ == '__main__':
    inp = Queue()
    for x in range(12):
        Process(target=worker, args=(inp, )).start()

    for x in range(4586, 33687):
        inp.put(x)

    for _ in range(12):
        inp.put('stop')
