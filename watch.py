from base64 import encodebytes
from CGRdb import load_schema, Molecule
from CGRdbUser import User
from json import dumps
from logging import info
from networkx import DiGraph
from pickle import load
from pony.orm import db_session

load_schema('patents', )
with open('zinc.pickle', 'rb') as z:
    zinc = load(z)


def b64(m):
    pass


def visualization(G, target):
    pass


def search_synth_paths(target, stages=3, first_number_of_paths=3):
    target.standardize()
    target.aromatize()
    if not Molecule.structure_exists(target):
        info('find_similar')
        similar_molecules = Molecule.find_similar(target, page=1, pagesize=1)
        info('done')
        for pair in similar_molecules:
            molecule = pair[0].structure
            if Molecule.find_reactions_by_product(molecule):
                info(f'Target = {molecule}, Index Tanimoto = {pair[-1]}')
                break
            else:
                info('Similar molecules not found as product')
                raise Exception

    sign = set(target)
    stack = [(target, stages)]
    while stack:
        reactions = Molecule.find_reactions_by_product(target)
