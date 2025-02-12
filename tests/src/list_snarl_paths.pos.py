import bdsg
import argparse


class Path:
    def __init__(self):
        self.nodes = []
        self.orients = []

    def addNode(self, node, orient):
        self.nodes.append(node)
        self.orients.append(orient)

    def addNodeHandle(self, node_h, stree):
        node_s = stree.net_handle_as_string(node_h)
        # trivial chain?
        if stree.is_trivial_chain(node_h):
            node_s = node_s.replace(' pretending to be a chain '
                                    'in a simple snarl', '')
        # parse node info
        node_s = node_s.replace('node ', '')
        node_o = '>'
        if 'rev' in node_s:
            node_o = '<'
        node_s = node_s.replace('rev', '').replace('fd', '')
        self.nodes.append(node_s)
        self.orients.append(node_o)

    def print(self):
        out_path = []
        for ii in range(len(self.nodes)):
            out_path.append(self.orients[ii] + str(self.nodes[ii]))
        return (''.join(out_path))

    def flip(self):
        self.nodes.reverse()
        self.orients.reverse()
        for ii in range(len(self.orients)):
            if self.nodes[ii] == '*':
                continue
            if self.orients[ii] == '>':
                self.orients[ii] = '<'
            else:
                self.orients[ii] = '>'

    def size(self):
        return (len(self.nodes))

    def nreversed(self):
        return (sum(['<' == orient for orient in self.orients]))


parser = argparse.ArgumentParser('List path through the netgraph of each snarl'
                                 ' in a pangenome')
parser.add_argument('-p', help='the input pangenome .pg file', required=True)
parser.add_argument('-d', help='the input distance index .dist file',
                    required=True)
parser.add_argument('-o', help='the output TSV file', required=True)
args = parser.parse_args()

# args = parser.parse_args(['-p', 'pg.pg', '-d', 'pg.full.dist', '-o',
#                           'pg.snarl_netgraph.paths.tsv'])

# load graph and snarl tree
pg = bdsg.bdsg.PackedGraph()
pg.deserialize(args.p)
ppo = bdsg.bdsg.PackedPositionOverlay(pg)
stree = bdsg.bdsg.SnarlDistanceIndex()
stree.deserialize(args.d)

# list snarls
# init with the child (only one ideally) of the root
root = stree.get_root()

# to save a list of (snarl handle, path name, position)
snarls = []

# to save the position of snarl/chain/etc using their string form as ID
snarls_pos = {}

# reference paths (TODO add to input parameters)
ref_paths = ['ref']


# given a node handle (dist index) return a position on a reference path
def get_node_position(node):
    node_h = stree.get_handle(node, pg)
    ret_pos = []

    def step_callback(step_handle):
        path_handle = pg.get_path_handle_of_step(step_handle)
        path_name = pg.get_path_name(path_handle)
        if path_name in ref_paths:
            position = ppo.get_position_of_step(step_handle)
            ret_pos.append(path_name)
            ret_pos.append(position)
            return False
        return True

    pg.for_each_step_on_handle(node_h, step_callback)
    return (ret_pos)


# given a handle for a node, snarl, chain, whatever, return a "start" position
def get_net_start_position(net):
    # if node, look at its position
    if stree.is_node(net):
        return (get_node_position(net))
    # otherwise check boundaries
    bnode1 = stree.get_bound(net, True, False)
    bnode1_p = get_node_position(bnode1)
    bnode2 = stree.get_bound(net, False, False)
    bnode2_p = get_node_position(bnode2)
    # if one of the bondaries is not on a reference path, return it
    if len(bnode1_p) == 0:
        return (bnode1_p)
    if len(bnode2_p) == 0:
        return (bnode2_p)
    assert bnode1_p[0] == bnode1_p[0], 'boundary nodes on different reference paths'
    # as "start" position, let's return the smallest position
    if bnode1_p[1] < bnode2_p[1]:
        return (bnode1_p)
    else:
        return (bnode2_p)


def save_snarl_tree_node(net):
    # check if we can find a position for this snarl on the reference path
    snarl_pos = get_net_start_position(net)
    # if we couldn't find a position, use the parent's that we should have
    # found and saved earlier
    if len(snarl_pos) == 0:
        par_net = stree.get_parent(net)
        snarl_pos = snarls_pos[stree.net_handle_as_string(par_net)]
    # save this position
    snarls_pos[stree.net_handle_as_string(net)] = snarl_pos
    # save snarl
    if stree.is_snarl(net):
        snarls.append([net, snarl_pos[0], snarl_pos[1]])
    # explore children
    if not stree.is_node(net) and not stree.is_sentinel(net):
        stree.for_each_child(net, save_snarl_tree_node)
    return (True)


stree.for_each_child(root, save_snarl_tree_node)
print('{} snarls found'.format(len(snarls)))

outf = open(args.o, 'wt')
outf.write('snarl\tpaths\tseqname\tpos\n')
npaths = 0

# for each snarl, lists paths and write to output TSV
for snarl_path_pos in snarls:
    snarl = snarl_path_pos[0]
    # create a snarl ID as LEFT_RIGTH bondary nodes
    sstart = stree.get_bound(snarl, False, True)
    sstart = stree.get_node_from_sentinel(sstart)
    send = stree.get_bound(snarl, True, True)
    send = stree.get_node_from_sentinel(send)
    snarl_id = '{}_{}'.format(stree.node_id(sstart),
                              stree.node_id(send))
    # init unfinished paths to the first boundary node
    paths = [[stree.get_bound(snarl, False, True)]]
    finished_paths = []
    # while there are still unifinshed paths, pop one path, follow net edges
    # to get the next step from the last node in the path, and then either
    # remember the completed path or remember to continue traversing
    while len(paths) > 0:
        path = paths.pop()
        # helper function to add the next child to the path we're building and
        # either add it to the list of completed paths or the list of paths
        # to continue building
        def add_to_path(next_child):
            if stree.is_sentinel(next_child):
                # If this is the bound of the snarl then we're done
                # Because we only traverse in the netgraph, it can only be the
                # bound of the parent snarl
                # Explicitly make a deep copy of the path because
                # idk how to do it in python
                finished_paths.append([])
                for net in path:
                    finished_paths[-1].append(net)
                finished_paths[-1].append(next_child)
            else:
                # If we reached sibling child of the snarl,
                # then continue the traversal
                # Do another copy into the other list of paths
                paths.append([])
                for net in path:
                    paths[-1].append(net)
                paths[-1].append(next_child)
            return True
        # run add_to_path for everything one step out
        # from the last thing in the path
        stree.follow_net_edges(path[-1], pg, False, lambda n: add_to_path(n))

    # prepare path list to output
    pretty_paths = []
    seq_net = []
    for path in finished_paths:
        ppath = Path()
        for net in path:
            if stree.is_sentinel(net):
                net = stree.get_node_from_sentinel(net)

            if stree.is_node(net):
                ppath.addNodeHandle(net, stree)
                node_start_id = stree.node_id(net)
                node_handle = pg.get_handle(node_start_id)
                seq_net.append("node " + pg.get_sequence(node_handle))

            if stree.is_trivial_chain(net):
                ppath.addNodeHandle(net, stree)
                stn_start = stree.get_bound(net, False, True)
                node_start_id = stree.node_id(stn_start)
                net_trivial_chain = pg.get_handle(node_start_id)
                seq_net.append("trivial_chain " + pg.get_sequence(net_trivial_chain))

            elif stree.is_chain(net):
                if stree.starts_at_start(net):
                    nodl = stree.get_bound(net, False, True)
                    nodr = stree.get_bound(net, True, False)
                else:
                    nodl = stree.get_bound(net, True, True)
                    nodr = stree.get_bound(net, False, False)
                nodl_s = stree.net_handle_as_string(nodl)
                nodr_s = stree.net_handle_as_string(nodr)
                ppath.addNodeHandle(nodl, stree)
                ppath.addNode('*', '>')
                ppath.addNodeHandle(nodr, stree)
                seq_net.append('chain _')

        # check if path is mostly traversing nodes in reverse orientation
        if ppath.nreversed() > ppath.size() / 2:
            ppath.flip()
            seq_net.reverse()
        pretty_paths.append(ppath.print())
    # write the paths out
    outf.write('{}\t{}\t{}\t{}\t{}\n'.format(snarl_id,
                                         ','.join(pretty_paths),
                                         ','.join(seq_net),
                                         snarl_path_pos[1],
                                         snarl_path_pos[2]))
    npaths += len(pretty_paths)

outf.close()

print('{} paths written in {}.'.format(npaths, args.o))
