from bdsg.bdsg import HashGraph

def create_simple_snp_graph(filename="simple_snp.hg"):
    gr = HashGraph()
    seqs = ["TATA", "CCGT", "C", "G", "GATAA", "T"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[4], nodes[5])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 4, 5]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("badpath")
    for idx in [0, 1, 3, 5]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)

def create_simple_loop_graph(filename="loop.hg"):
    gr = HashGraph()
    seqs = ["A", "C", "T", "G"]
    nodes = [gr.create_handle(s) for s in seqs]

    # A -> C -> T -> G -> C (loop)
    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[3], nodes[1])

    path = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3, 1, 2]:
        gr.append_step(path, nodes[idx])

    gr.serialize(filename)

def create_fork_graph(filename="fork.hg"):
    gr = HashGraph()
    seqs = ["AAA", "CCC", "GGG", "TTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    # AAA -> {CCC, GGG}
    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[0], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[3])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 3]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("fork_path2")
    for idx in [0, 2, 3]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)

def create_linear_path(filename="linear.hg"):
    gr = HashGraph()
    seqs = ["AAA", "T", "GGG", "C"]
    nodes = [gr.create_handle(s) for s in seqs]

    for i in range(len(nodes) - 1):
        gr.create_edge(nodes[i], nodes[i + 1])

    path = gr.create_path_handle("ref")
    for node in nodes:
        gr.append_step(path, node)

    gr.serialize(filename)

def create_bubble_graph(filename="bubble.hg"):
    gr = HashGraph()
    seqs = ["A", "G", "C", "T", "A"]
    nodes = [gr.create_handle(s) for s in seqs]

    # Bubble: A -> G/C -> T -> A
    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[0], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[3], nodes[4])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 3, 4]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("bubble_alt")
    for idx in [0, 2, 3, 4]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)

def create_4th_graph(filename="4th.hg"):
    gr = HashGraph()
    seqs = ["AAA", "GA", "CG", "TCC", "CT", "AA"]
    nodes = [gr.create_handle(s) for s in seqs]

    # 4TH :   AAA -> GA -> CG -> TCC
    #                GA -> AA
    #         AAA -> CT -> AA -> TCC
    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[0], nodes[4])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[5])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[4], nodes[5])
    gr.create_edge(nodes[5], nodes[3])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("bubble_alt")
    for idx in [0, 4, 5, 3]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)

def create_4th_modify_graph(filename="4th_modify.hg"):
    gr = HashGraph()
    seqs = ["TTT", "AAA", "GT", "CATG", "TCC", "CTTTT", "A", "TA"]
    nodes = [gr.create_handle(s) for s in seqs]

    # 4TH :   AAA -> GA -> CG -> TCC
    #                GA -> AA
    #         AAA -> CT -> AA -> TCC
    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[5])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[2], nodes[6])
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[5], nodes[6])
    gr.create_edge(nodes[6], nodes[4])
    gr.create_edge(nodes[4], nodes[7])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3, 4, 7]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("bubble_alt")
    for idx in [0, 1, 5, 6, 4, 7]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)
    
# Example call to generate all
if __name__ == "__main__":
    create_simple_snp_graph()
    create_simple_loop_graph()
    create_fork_graph()
    create_linear_path()
    create_bubble_graph()
    create_4th_graph()
    create_4th_modify_graph()

# vg convert simple_snp.hg > simple_snp.pg
# vg index simple_snp.pg -j simple_snp.dist
# vg convert bubble.hg > bubble.pg
# vg index bubble.pg -j bubble.dist
# vg convert fork.hg > fork.pg
# vg index fork.pg -j fork.dist
# vg convert linear.hg > linear.pg
# vg index linear.pg -j linear.dist
# vg convert loop.hg > loop.pg
# vg index loop.pg -j loop.dist
# vg convert 4th.hg > 4th.pg
# vg index 4th.pg -j 4th.dist

# vg find -x 4th_modify.pg -r 1:6 -c 10 | vg view -dp - | dot -Tsvg -o 4th_modify.svg
