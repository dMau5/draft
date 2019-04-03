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
    return encodebytes(m.depict().encode()).decode().replace('\n', '')


def visualization(G, target):
    template = 'data:image/svg+xml;base64,'
    info('go to visualization')
    data = {}
    _all = {target: 1}
    n = 2
    nodes = [{'id': 1, 'image': template + b64(target), 'shape': 'image', 'size': 100,
              "shapeProperties": {"useBorderWithImage": True}, "color": {"border": "red", "background": "white"}}]
    with db_session:
        for node in G.nodes():
            if node == target:
                continue
            if not isinstance(node, int):
                _all.update({node: n})
                color_border = "white"
                if 'zinc' in G.nodes[node]:
                    color_border = "green"
                item = {'id': n, 'image': template + b64(node), 'shape': 'image', 'size': 35,
                        "shapeProperties": {"useBorderWithImage": True},
                        "color": {"border": color_border, "background": "white"}}
            else:
                descr = ', '.join(G.node[node]['data'])
                _all.update({node: descr})
                item = {'id': descr, 'label': f'{n}',
                        'color': {'border': 'black',
                                  'background': 'white'},
                        'shape': 'box',
                        'widthConstraint': {'maximum': 300}, 'font': {'size': 20}}
            nodes.append(item)
            n += 1
        data['nodes'] = nodes
        edges = []
        for e in G.edges():
            item = {'width': 2}
            source = e[0]
            item['from'] = _all[source]
            item['to'] = _all[e[-1]]
            if not isinstance(source, int):
                item['color'] = {'color': 'blue'}
            else:
                item['color'] = {'color': 'red'}
            edges.append(item)
        data['edges'] = edges
        _json = dumps(data)
    return "<html>\n" \
           "<head>\n" \
           "    <meta charset='utf-8'>\n" \
           "    <title>Hierarchical Layout without Physics</title>\n" \
           "    <script type='text/javascript' src='static/vis/dist/vis.js'></script>\n" \
           "    <link href='static/vis/dist/vis-network.min.css' rel='stylesheet' type='text/css' />\n" \
           "    <style type='text/css'>\n" \
           "        #network{\n" \
           "            width: 80%;\n" \
           "            height: 100%;\n" \
           "        }\n" \
           "    </style>\n" \
           "    <style type='text/css'>\n" \
           "    #selection{\n" \
           "        position: absolute;\n" \
           "        right: 0;\n" \
           "        padding: 10% 5px;\n" \
           "        border-left: 1px solid blue;\n" \
           "        width: 19%;\n" \
           "        height: 100%;\n" \
           "    }\n" \
           "    </style>\n" \
           "    <style>\n" \
           "    button {\n" \
           "    position: static;\n" \
           "    border: 1px black solid;\n" \
           "    border-radius: 3px;\n" \
           "    padding: 5px;\n" \
           "    background: white;\n" \
           "    }\n" \
           "    button:hover {\n" \
           "    background: silver;\n" \
           "    }\n" \
           "    </style>\n" \
           "</head>\n" \
           "<body>\n" \
           "<button type='button' onclick=\"window.location.href = '/index'\">Go back</button>\n" \
           "<div id='selection'></div>\n" \
           "<div id='network'></div>\n" \
           "<script>\n" \
           "    var nodesArray, nodes, edgesArray, edges, network;\n" \
           f"    var json = {_json};\n" \
           "    nodesArray =  json.nodes;\n" \
           "    nodes = new vis.DataSet(nodesArray);\n" \
           "    edgesArray = json.edges;\n" \
           "    edges = new vis.DataSet(edgesArray);\n" \
           "    var data = {\n" \
           "        nodes: nodes,\n" \
           "        edges: edges\n" \
           "    };\n" \
           "    var container = document.getElementById('network');\n" \
           "    var options = {\n" \
           "        layout: {\n" \
           "            hierarchical: {\n" \
           "                direction: 'LR',\n" \
           "                sortMethod: 'directed',\n" \
           "                levelSeparation: 530,\n" \
           "                nodeSpacing: 150,\n" \
           "                parentCentralization: true,\n" \
           "            }\n" \
           "        },\n" \
           "        interaction: {dragNodes :true},\n" \
           "        physics: {\n" \
           "            enabled: true\n" \
           "        },\n" \
           "    };\n" \
           "    network = new vis.Network(container, data, options);\n" \
           "    network.on('select', function (params) {\n" \
           "    document.getElementById('selection').innerHTML = params.nodes;\n" \
           " });\n" \
           "</script>\n" \
           "</body>" \
           "</html>"


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

    sign = set()
    stack = [(target, stages)]

