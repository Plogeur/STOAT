from bdsg.bdsg import HashGraph

def create_simple_snp_graph(filename="simple_snp.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "C", "G", "AAAA", "TTTT"]
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

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 3, 5]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)

def create_insert_snp_graph(filename="insert_snp.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "C", "GTT", "AAAA", "TTTT"]
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

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 3, 4, 5]:
        gr.append_step(path2, nodes[idx])

    gr.serialize(filename)

def create_deletion_snp_graph(filename="deletion_snp.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "C", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[3], nodes[4])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3, 4]:
        gr.append_step(path1, nodes[idx])

    gr.serialize(filename)

def create_insert_deletion_graph(filename="insert_deletion.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "GTT", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[2], nodes[3])
    gr.create_edge(nodes[3], nodes[4])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3, 4]:
        gr.append_step(path1, nodes[idx])

    gr.serialize(filename)

def create_simple_loop_graph(filename="loop.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "A", "C", "T", "G", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[2], nodes[1])
    gr.create_edge(nodes[1], nodes[3])
    gr.create_edge(nodes[3], nodes[4])
    gr.create_edge(nodes[4], nodes[5])

    path = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 3, 5]:
        gr.append_step(path, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 2, 5]:
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

def create_snp_and_nested_snp_graph(filename="snp_and_nested_snp.hg"):
    gr = HashGraph()
    seqs = ["TTTT", "AAAA", "T", "G", "A", "C", "T", "AAAA", "TTTT"]
    nodes = [gr.create_handle(s) for s in seqs]

    gr.create_edge(nodes[0], nodes[1])
    gr.create_edge(nodes[1], nodes[2])
    gr.create_edge(nodes[2], nodes[3]) # NESTED
    gr.create_edge(nodes[2], nodes[4])
    gr.create_edge(nodes[3], nodes[5]) # NESTED
    gr.create_edge(nodes[3], nodes[6]) # NESTED
    gr.create_edge(nodes[5], nodes[7]) # NESTED
    gr.create_edge(nodes[6], nodes[7]) # NESTED
    gr.create_edge(nodes[7], nodes[8]) # NESTED
    gr.create_edge(nodes[4], nodes[8])
    gr.create_edge(nodes[8], nodes[9])

    path1 = gr.create_path_handle("ref")
    for idx in [0, 1, 2, 4, 8, 9]:
        gr.append_step(path1, nodes[idx])

    path2 = gr.create_path_handle("alt")
    for idx in [0, 1, 2, 3, 5, 6, 7, 8, 9]:
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
    create_insert_snp_graph()
    create_deletion_snp_graph()
    create_insert_deletion_graph()
    create_snp_and_nested_snp_graph()
    create_simple_loop_graph()
    create_linear_path()
    create_4th_modify_graph()

# vg convert simple_snp.hg > simple_snp.pg
# vg index simple_snp.pg -j simple_snp.dist
# vg find -x 4th_modify.pg -r 1:6 -c 10 | vg view -dp - | dot -Tsvg -o 4th_modify.svg
