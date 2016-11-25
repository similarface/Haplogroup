__author__ = 'similarface'
import sys
import json
import UserDict


import sys
import json
import UserDict

class Node(object):
    def __init__(self, id, parent, name,rsid):
        self.id = id
        self.parent = parent
        self.children = []
        self.name = name
        self.rsid=rsid


class NodeDict(UserDict.UserDict):
    def addNodes(self, nodes):
        """ Add every node as a child to its parent by doing two passes."""
        for i in (1, 2):
            for node in nodes:
                self.data[node.id] = node
                if node.parent in self.data.keys():
                    if node.parent != "none" and node not in self.data[node.parent].children:
                        self.data[node.parent].children.append(node)

class NodeJSONEncoder(json.JSONEncoder):
    def default(self, node):
        if type(node) == Node:
            return {"iid":node.id,"pid":node.parent, "name":node.name,"rsid":node.rsid, "children":node.children}
        raise TypeError("{} is not an instance of Node".format(node))

if __name__ == "__main__":
    nodes = []
    filename="/Users/similarface/Documents/wocao.txt"
    with open(filename) as f:
        for row in f.readlines()[1:]:
            #nid, parent, name = row.split()
            id,pid,name,rsid=row.split('\t')
            #nodes.append(Node(nid, parent, name))
            nodes.append(Node(id,pid,name,rsid.strip().split(',')))

    nodeDict = NodeDict()
    nodeDict.addNodes(nodes)

    rootNodes = [node for id, node in nodeDict.items() if node.parent =="1" ]
    for rootNode in rootNodes:
        print NodeJSONEncoder().encode(rootNode)