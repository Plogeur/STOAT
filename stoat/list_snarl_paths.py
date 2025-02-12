import bdsg # type: ignore
import argparse
from stoat import utils
import time 
import os 

# class to help make paths from BDSG objects
# and deal with orientation, flipping, etc
class Path:
    def __init__(self):
        self.nodes = []
        self.orients = []

    def addNode(self, node, orient):
        # here we know the actual node id and orientation
        self.nodes.append(node)
        self.orients.append(orient)

    def addNodeHandle(self, node_h, stree):
        # we have a BDSG handle and need to extract node id and orientation
        # couldn't find a better way than parsing the
        # string representation of the node...
        node_s = stree.net_handle_as_string(node_h)
        # trivial chain
        if stree.is_trivial_chain(node_h):
            node_s = node_s.replace(' pretending to be a chain', '')
            node_s = node_s.replace(' in a simple snarl', '')

        # parse node info
        node_s = node_s.replace('node ', '')
        node_o = '>'
        if 'rev' in node_s:
            node_o = '<'
        node_s = node_s.replace('rev', '').replace('fd', '')
        # add node to path
        self.nodes.append(node_s)
        self.orients.append(node_o)

    def print(self):
        # write the string representation of the path
        # e.g. ">I>J<K>L" or ">I>J>*>M>N"
        out_path = []
        for ii in range(len(self.nodes)):
            out_path.append(self.orients[ii] + str(self.nodes[ii]))
        return (''.join(out_path))

    def flip(self):
        # some paths are traversing the nodes entirely or mostly in reverse
        # for convenience we might can to flip them
        self.nodes.reverse()
        self.orients.reverse()
        for ii in range(len(self.orients)):
            if self.nodes[ii] == '*':
                # don't flip the orientation of the * because
                # it should be ">" by construction
                continue
            if self.orients[ii] == '>':
                self.orients[ii] = '<'
            else:
                self.orients[ii] = '>'

    def size(self):
        return (len(self.nodes))

    def nreversed(self):
        # counts how many nodes are traversed in reverse
        return (sum(['<' == orient for orient in self.orients]))

def calcul_type_variant(list_list_length_paths) :
    """ 
    Calcul type variant of a tested snarl
    """
    list_type_variant = []
    for path_lengths in list_list_length_paths :
        
        # Case snarl in snarl or 4 node snarl
        if len(path_lengths) > 3 or path_lengths[1] == '-1' :
            list_type_variant.append("COMPLEX")

        # Case simple path len 3
        elif len(path_lengths) == 3 :
            list_type_variant.append("SNP" if path_lengths[1] == "1" else "INS")

        # length < 3 / Deletion
        else :
            list_type_variant.append("DEL")

    return list_type_variant

def check_threshold(proportion) :
    proportion = float(proportion)
    if proportion <= 0 :
        raise argparse.ArgumentTypeError("Proportion value must be >0.")

    return proportion

def find_snarl_id(stree, snarl) :
    # create a snarl ID as LEFT_RIGTH bondary nodes
    sstart = stree.get_bound(snarl, False, True)
    sstart = stree.get_node_from_sentinel(sstart)
    send = stree.get_bound(snarl, True, True)
    send = stree.get_node_from_sentinel(send)
    snarl_id = '{}_{}'.format(stree.node_id(send), stree.node_id(sstart))

    return snarl_id

def follow_edges(stree, finished_paths, path, paths, pg) :
    def add_to_path(next_child) :

        if stree.is_sentinel(next_child):
            # If this is the bound of the snarl then we're done
            # Because we only traverse in the netgraph, it can only be the
            # bound of the parent snarl
            finished_paths.append([])
            for net in path:
                finished_paths[-1].append(net)
            finished_paths[-1].append(next_child)
        else :
            for i in path : 
                # Case where we find a loop 
                if stree.net_handle_as_string(i) == stree.net_handle_as_string(next_child) :
                    return False
            paths.append([])
            for net in path:
                paths[-1].append(net)
            paths[-1].append(next_child)
        return True

    # from the last thing in the path
    stree.follow_net_edges(path[-1], pg, False, lambda n: add_to_path(n))
 
def save_snarls(stree, root, pg, ref_paths, ppo) :
    snarls = []
    snarls_pos = {}

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
    
    stree.for_each_child(root, save_snarl_tree_node)
    return snarls

def parse_graph_tree(pg_file, dist_file) :

    # load graph and snarl tree
    pg = bdsg.bdsg.PackedGraph()
    pg.deserialize(pg_file)
    stree = bdsg.bdsg.SnarlDistanceIndex()
    stree.deserialize(dist_file)
    pp_overlay = bdsg.bdsg.PackedPositionOverlay(pg)

    # list all snarls in pangenome
    # init with the child (only one ideally) of the root
    root = stree.get_root()
    return stree, pg, root, pp_overlay
   
def fill_pretty_paths(stree, pg, finished_paths) :
    pretty_paths = []
    length_net_paths = []

    for path in finished_paths:
        ppath = Path()
        length_net = []

        for net in path :
            if stree.is_sentinel(net) :
                net = stree.get_node_from_sentinel(net)

            # case node : get the node length
            if stree.is_node(net) :
                ppath.addNodeHandle(net, stree)
                length_net.append(str(stree.node_length(net)))

            # case trivial_chain : get the first node length
            if stree.is_trivial_chain(net) :
                ppath.addNodeHandle(net, stree)
                stn_start = stree.get_bound(net, False, True)
                node_start_id = stree.node_id(stn_start)
                net_trivial_chain = pg.get_handle(node_start_id)
                length_net.append(str(pg.get_length(net_trivial_chain)))

            elif stree.is_chain(net) :
                # if it's a chain, we need to write someting like ">Nl>*>Nr"
                if stree.starts_at_start(net):
                    nodl = stree.get_bound(net, False, True)
                    nodr = stree.get_bound(net, True, False)
                else:
                    nodl = stree.get_bound(net, True, True)
                    nodr = stree.get_bound(net, False, False)
                ppath.addNodeHandle(nodl, stree)
                ppath.addNode('*', '>')
                ppath.addNodeHandle(nodr, stree)
                length_net.append("-1")

        # check if path is mostly traversing nodes in reverse orientation
        if ppath.nreversed() > ppath.size() / 2 :
            ppath.flip()
        pretty_paths.append(ppath.print()) 
        length_net_paths.append(length_net)

    type_variants = calcul_type_variant(length_net_paths)
    assert len(type_variants) == len(pretty_paths)
    return pretty_paths, type_variants

def loop_over_snarls_write(stree, snarls, pg, output_file, output_snarl_not_analyse, children_treshold=50, bool_return=True) :

    with open(output_file, 'w') as out_snarl, open(output_snarl_not_analyse, 'w') as out_fail:
        out_snarl.write('snarl\tpaths\ttype\tchr\tpos\n')
        out_fail.write('snarl\treason\n')

        snarl_paths = []
        paths_number_analysis = 0
        time_threshold = 2 # 2s max per snarl analysis

        children = [0]
        def count_children(net):
            children[0] += 1
            return (True)

        # for each snarl, lists paths through the netgraph and write to output TSV
        for snarl_path_pos in snarls:
            
            snarl = snarl_path_pos[0]
            snarl_time = time.time()
            snarl_id = find_snarl_id(stree, snarl)
            not_break = True
            children = [0]

            stree.for_each_child(snarl, count_children)
            if children[0] > children_treshold :
                out_fail.write('{}\t{}\n'.format(snarl_id, "too_many_children"))
                continue

            # we'll traverse the netgraph starting at the left boundary
            # init unfinished paths to the first boundary node
            paths = [[stree.get_bound(snarl, False, True)]]
            finished_paths = []
            while len(paths) > 0 :
                path = paths.pop()
        
                if time.time() - snarl_time > time_threshold :
                    out_fail.write('{}\t{}\n'.format(snarl_id, "time_calculation_out"))
                    not_break = False
                    break

                follow_edges(stree, finished_paths, path, paths, pg)

            if not_break :

                # prepare path list to output and write each path directly to the file
                pretty_paths, type_variants = fill_pretty_paths(stree, pg, finished_paths)
                out_snarl.write('{}\t{}\t{}\t{}\t{}\n'.format(snarl_id, ','.join(pretty_paths), ','.join(type_variants), snarl_path_pos[1], snarl_path_pos[2]))

                if bool_return :
                    snarl_paths.append((snarl_id, pretty_paths, ','.join(type_variants), snarl_path_pos[1], snarl_path_pos[2]))
                
                paths_number_analysis += len(pretty_paths)

    return snarl_paths, paths_number_analysis

if __name__ == "__main__" :

    parser = argparse.ArgumentParser('List path through the netgraph of each snarl in a pangenome')
    parser.add_argument('-p', type=utils.check_file, help='The input pangenome .pg file', required=True)
    parser.add_argument('-d', type=utils.check_file, help='The input distance index .dist file', required=True)
    parser.add_argument('-c', "--chr", type=utils.check_file, help='The input reference chr file', required=False)
    parser.add_argument("-t", type=check_threshold, help='Children threshold', required=False)
    parser.add_argument('-o', help='output dir', type=str, required=False)
    args = parser.parse_args()

    reference = utils.parse_chr_reference(args.chr) if args.chr else "ref"
    output_dir = args.o or "output"
    os.makedirs(output_dir, exist_ok=True)
    output = os.path.join(output_dir, "list_snarl_paths.tsv")
    output_snarl_not_analyse = os.path.join(output_dir, "snarl_not_analyse.tsv")
    stree, pg, root, pp_overlay = parse_graph_tree(args.p, args.d)

    snarls = save_snarls(stree, root, pg, reference, pp_overlay)
    print(f"Total of snarls found : {len(snarls)}")
    print("Saving snarl path decomposition...")

    threshold = args.t if args.t else 10
    _, paths_number_analysis = loop_over_snarls_write(stree, snarls, pg, output, output_snarl_not_analyse, threshold, False)
    print(f"Total of paths analyse : {paths_number_analysis}")

    # python3 stoat/list_snarl_paths.py -p /home/mbagarre/Bureau/droso_data/fly/fly.pg -d /home/mbagarre/Bureau/droso_data/fly/fly.dist -o output/test
    # vg find -x ../snarl_data/fly.gbz -r 5176878:5176884 -c 10 | vg view -dp - | dot -Tsvg -o ../snarl_data/subgraph.svg

    # binary 
    # python3 stoat/list_snarl_paths.py -p tests/simulation/binary_data/pg.pg -d tests/simulation/binary_data/pg.dist -c tests/simulation/binary_data/pg.chromosome -o output/test
    
    # quantitative
    # python3 stoat/list_snarl_paths.py -p tests/simulation/quantitative_data/pg.pg -d tests/simulation/quantitative_data/pg.dist -c tests/simulation/quantitative_data/pg.chromosome -o output/test
