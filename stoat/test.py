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

def find_node_position_and_chromosome(pg, pp_overlay, handle_t):
    positions = []
    chromosomes = []

    def step_callback(step_handle):
        path_handle = pg.get_path_handle_of_step(step_handle)
        path_name = pg.get_path_name(path_handle)
        if path_name not in chromosomes:
            chromosomes.append(path_name)
        position = pp_overlay.get_position_of_step(step_handle)
        positions.append(position)
        return True

    pg.for_each_step_on_handle(handle_t, step_callback)
    try :
        return chromosomes[-1], positions[-1]
    except :
        return "-1","-1"

def fill_pretty_paths(stree, pg, pp_overlay, finished_paths) :
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
                chromosome, position = find_node_position_and_chromosome(pg, pp_overlay, node_handle_t)
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
