import bdsg # type: ignore
import argparse
from stoat import utils
import re
from collections import defaultdict
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

def split_paths(path) :
    return re.findall(r'\d+', path)

def length_node(pg, node_id) :
    return pg.get_length(node_id)

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
            list_type_variant.append("SNP" if path_lengths[1] == 1 else "INS")

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
    stree.follow_net_edges(path[-1], pg, False, add_to_path)

def save_snarls(stree, root) :

    # list storing the snarl objects
    snarls = []

    def save_snarl_tree_node(net):
        if stree.is_snarl(net):
            snarls.append(net)

        if not stree.is_node(net) and not stree.is_sentinel(net):
            stree.for_each_child(net, save_snarl_tree_node)
        return (True)
    
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

def find_node_position_and_chromosome(pg, pp_overlay, reference, handle_t):
    positions = ["-1"]
    chromosomes = ["-1"]
    
    def step_callback(step_handle):
        path_handle = pg.get_path_handle_of_step(step_handle)
        path_name = pg.get_path_name(path_handle)
        print("path_name : ", path_name)

        while path_name not in reference :
            node_handle_t = pg.get_handle_of_step(step_handle)
            pg.for_each_step_on_handle(node_handle_t, step_callback)
        
        chromosomes.append(path_name)
        position = pp_overlay.get_position_of_step(step_handle)
        positions.append(position)
        return True
    
    pg.for_each_step_on_handle(handle_t, step_callback)
    return chromosomes[-1], positions[-1]
   
def fill_pretty_paths(stree, pg, pp_overlay, reference, finished_paths) :
    pretty_paths = []
    length_net_paths = []
    chromosomes = []
    positions = []

    for path in finished_paths:
        ppath = Path()
        length_net = []

        # stree : bdsg.bdsg.SnarlDistanceIndex
        # net : bdsg.handlegraph.net_handle_t
        for net in path :
            if stree.is_sentinel(net) :
                net = stree.get_node_from_sentinel(net)
 
            # case node : get the node length
            if stree.is_node(net) :
                ppath.addNodeHandle(net, stree)
                length_net.append(str(stree.node_length(net)))
                node_handle_t = stree.get_handle(net, pg)
                chromosome, position = find_node_position_and_chromosome(pg, pp_overlay, reference, node_handle_t)
                chromosomes.append(chromosome)
                positions.append(position)

            # case trivial_chain : get the first node length
            elif stree.is_trivial_chain(net) :
                ppath.addNodeHandle(net, stree)
                stn_start = stree.get_bound(net, False, True)
                node_start_id = stree.node_id(stn_start)
                net_trivial_chain = pg.get_handle(node_start_id)
                length_net.append(str(pg.get_length(net_trivial_chain)))

            elif stree.is_chain(net) :
                # if it's a chain, we need to write someting like ">Nl>*>Nr"
                nodl = stree.get_bound(net, False, True)
                nodr = stree.get_bound(net, True, False)
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
    return pretty_paths, type_variants, chromosomes[0], positions[0]

def write_header_output(output_file) :
    with open(output_file, 'w') as outf:
        outf.write('snarl\tpaths\ttype\tchr\tpos\n')

def write_output(output_file, snarl_id, pretty_paths, type_variants, chromosome, position) :
    with open(output_file, 'a') as outf:
        outf.write('{}\t{}\t{}\t{}\t{}\n'.format(snarl_id, ','.join(pretty_paths), ','.join(type_variants), chromosome, position))

def write_header_output_not_analyse(output_file) :
    with open(output_file, 'w') as outf:
        outf.write('snarl\treason\n')

def write_output_not_analyse(output_file, snarl_id, reason) :
    with open(output_file, 'a') as outf:
        outf.write('{}\t{}\n'.format(snarl_id, reason))

def loop_over_snarls_write(stree, snarls, pg, pp_overlay, reference, output_file, output_snarl_not_analyse, children_treshold=50, bool_return=True) :

    write_header_output(output_file)
    write_header_output_not_analyse(output_snarl_not_analyse)

    snarl_paths = defaultdict(list)
    paths_number_analysis = 0
    time_threshold = 2 # 2s max per snarl analysis

    children = [0]
    def count_children(net):
        children[0] += 1
        return (True)

    # for each snarl, lists paths through the netgraph and write to output TSV
    for snarl in snarls:
    
        snarl_time = time.time()
        snarl_id = find_snarl_id(stree, snarl)
        not_break = True
        children = [0]

        stree.for_each_child(snarl, count_children)
        if children[0] > children_treshold :
            write_output_not_analyse(output_snarl_not_analyse, snarl_id, "too_many_children")
            continue

        # we'll traverse the netgraph starting at the left boundary
        # init unfinished paths to the first boundary node
        paths = [[stree.get_bound(snarl, False, True)]]
        finished_paths = []
        while len(paths) > 0 :
            path = paths.pop()
    
            if time.time() - snarl_time > time_threshold :
                write_output_not_analyse(output_snarl_not_analyse, snarl_id, "time_calculation_out")
                not_break = False
                break

            follow_edges(stree, finished_paths, path, paths, pg)

        if not_break :
            # prepare path list to output and write each path directly to the file
            pretty_paths, type_variants, chromosome, position = fill_pretty_paths(stree, pg, pp_overlay, reference, finished_paths)
            write_output(output_file, snarl_id, pretty_paths, type_variants, chromosome, position)

            if bool_return :
                snarl_paths[snarl_id].extend(pretty_paths)
            
            paths_number_analysis += len(pretty_paths)

    return snarl_paths, paths_number_analysis

if __name__ == "__main__" :

    parser = argparse.ArgumentParser('List path through the netgraph of each snarl in a pangenome')
    parser.add_argument('-p', type=utils.check_file, help='The input pangenome .pg file', required=True)
    parser.add_argument('-d', type=utils.check_file, help='The input distance index .dist file', required=True)
    parser.add_argument('-c', "--chr", type=utils.check_file, help='The input reference chr file', required=True)
    parser.add_argument("-t", type=check_threshold, help='Children threshold', required=False)
    parser.add_argument('-o', help='output file', type=str, required=False)
    args = parser.parse_args()

    reference = utils.parse_chr_reference(args.chr)
    output_dir = args.o or "output"    
    os.makedirs(output_dir, exist_ok=True)
    output = os.path.join(output_dir, "list_snarl_paths.tsv")
    output_snarl_not_analyse = os.path.join(output_dir, "snarl_not_analyse.tsv")

    stree, pg, root, pp_overlay = parse_graph_tree(args.p, args.d)
    snarls = save_snarls(stree, root)
    print(f"Total of snarls found : {len(snarls)}")
    print("Saving snarl path decomposition...")

    threshold = args.t if args.t else 10
    _, paths_number_analysis = loop_over_snarls_write(stree, snarls, pg, pp_overlay, reference, output, output_snarl_not_analyse, threshold, False)
    print(f"Total of paths analyse : {paths_number_analysis}")

    # python3 stoat/list_snarl_paths.py -p /home/mbagarre/Bureau/droso_data/fly/fly.pg -d /home/mbagarre/Bureau/droso_data/fly/fly.dist -o output/test/test_list_snarl.tsv
    # vg find -x ../snarl_data/fly.gbz -r 5176878:5176884 -c 10 | vg view -dp - | dot -Tsvg -o ../snarl_data/subgraph.svg

    # python3 stoat/list_snarl_paths.py -p tests/simulation/binary_data/pg.full.pg -d tests/simulation/binary_data/pg.dist -c tests/simulation/binary_data/pg.chromosome -o output/test/test_list_snarl.tsv
