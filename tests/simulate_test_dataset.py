import random
import argparse
# -*- coding: utf-8 -*-

# Written by Jean Monlong (jean.monlong@inserm.fr)
# Simulate a simple pangenome graph with SNPs and indels. Some indels can have
# nested SNPs too. See below for the parameters used.

# Briefly, the pangenome is first build, one variant site at a time. Then, it
# is traversed by one (random) path that will be used as reference path. A
# first GFA file can be written at this point. At each beginning "fork" of a
# variant site, we had set random probability to follow one edge or the other.
# To create, difference between two groups of samples, we then select some
# variants and slightly shift those probabilities for the second group. Then
# we simulate paths for each haplotype of each sample and each group following
# those probabilities. A GFA file including all of those haplotype paths can
# also be outputted (to be used later to simulate reads from each sample)

# Limitation: simple pangenome with only "biallelic" snarls (two out-edges max)
# and only one level of nesting.

def randSeq(length):
    """Simulate a random nucleotide sequence"""
    return ''.join(random.choices(_nuc, k=length))

class Graph:
    def __init__(self, phenotype_type_='binary'):
        # phenotype type
        self.phenotype_type = phenotype_type_
        # map node name to sequence
        self.nodes = {}
        # map predecessor nodes to their successor nodes
        self.edges = {}
        # map path names to list of nodes
        self.paths = {}
        # map a snarl start to frequencies
        self.snarls_freq = {}
        # phenotype for each sample (sample : phenotype)
        self.phenotypes = []
        # keep track of a few things
        self.next_node_id = 1
        self.ngroups = 2

    def quantitativePhenotype(self, samp):
        """Simulate a quantitative phenotype for a sample"""
        samp_t = nsamp//2
        for samp in range(0, nsamp):
            if samp < samp_t:
                self.phenotypes.append(random.uniform(0.0, 1.0))
            else:
                self.phenotypes.append(random.uniform(0.0, -1.0))

    def binaryPhenotype(self, samp):
        """Simulate a binary phenotype for a sample"""
        samp_t = nsamp//2
        for samp in range(0, nsamp):
            if samp < samp_t:
                self.phenotypes.append(0)
            else:
                self.phenotypes.append(1)
    
    def addNode(self, pred_nodes=[], min_size=50, max_size=300):
        # create a node and connect it to specified predecessors
        seq = randSeq(random.randint(min_size, max_size))
        # add node
        new_node = self.next_node_id
        self.nodes[new_node] = seq
        # add edges
        for pnod in pred_nodes:
            self.addEdge(pnod, new_node)
        # update next node ID
        self.next_node_id += 1
        return (new_node)

    def addEdge(self, pnode, snode):
        if pnode not in self.edges:
            self.edges[pnode] = {}
        self.edges[pnode][snode] = True

    def addSNP(self, pred_node):
        # pick two different alleles for the SNP
        als = random.sample(_nuc, 2)
        # create nodes and edges to predecessors
        node_a = self.next_node_id
        self.nodes[node_a] = als[0]
        self.next_node_id += 1
        node_b = self.next_node_id
        self.nodes[node_b] = als[1]
        self.next_node_id += 1
        # add edges
        self.addEdge(pred_node, node_a)
        self.addEdge(pred_node, node_b)
        # create a successor node
        suc_node = self.addNode([node_a, node_b])
        # init snarl frequency (same in all groups)
        freqs = {}
        freq_init = random.random()
        freqs[node_a] = [freq_init] * self.ngroups
        freqs[node_b] = [1 - freq_init] * self.ngroups
        self.snarls_freq[pred_node] = freqs
        return (suc_node)

    def addIndel(self, pred_node, snp_prob=.5):
        # create a middle node
        mid_node = self.addNode([pred_node])
        final_mid_node = mid_node
        # create a nested SNP sometimes
        while random.random() < snp_prob:
            final_mid_node = self.addSNP(final_mid_node)
        # create a successor node
        suc_node = self.addNode([final_mid_node])
        # add edge from predecessor directly to successor
        self.addEdge(pred_node, suc_node)
        # init snarl frequency (same in all groups)
        freqs = {}
        freq_init = random.random()
        freqs[mid_node] = [freq_init] * self.ngroups
        freqs[suc_node] = [1 - freq_init] * self.ngroups
        self.snarls_freq[pred_node] = freqs
        return (suc_node)

    def addPath(self, path_name='ref', group=0, samp=0):
        # start at node 1
        path = [1]
        while path[-1] in self.edges:
            # look for successor nodes
            snodes = self.edges[path[-1]]
            snodes = list(snodes.keys())
            # if there are some, pick one randomly
            if len(snodes) > 0:
                if path[-1] in self.snarls_freq:
                    freqs = self.snarls_freq[path[-1]]
                    rr = random.random()
                    
                    # pick next node based on phenotype
                    if self.phenotype_type == 'binary':
                        tot_freq = 0
                        # choose path base on group
                        for snode in freqs:
                            tot_freq += freqs[snode][group]
                            if rr < tot_freq:
                                break
                    else :
                        sample_pheno = (self.phenotypes[samp] / 2) + 0.5  # [0,1]
                        for snode in freqs:
                            if rr < sample_pheno:
                                break

                    path.append(snode)
                else:
                    path.append(random.sample(snodes, 1)[0])
        self.paths[path_name] = path

    def createMarkerFreq(self, prop_markers=.1, min_dev=.01, max_dev=.5):
        # loop over snarls and change a group's frequency sometimes
        for pnode in self.snarls_freq:
            # skip if not selected as marker
            if random.random() > prop_markers:
                continue
            # pick a group to deviate
            gp = random.randint(1, self.ngroups) - 1
            # pick a deviation size
            dev = random.uniform(min_dev, max_dev)
            # apply deviation
            dev_added = False
            for snode in self.snarls_freq[pnode]:
                if self.snarls_freq[pnode][snode][gp] < 1 - dev and not dev_added:
                    self.snarls_freq[pnode][snode][gp] += dev
                    dev_added = True
                else:
                    self.snarls_freq[pnode][snode][gp] -= dev

    def writeGfa(self, out_fn):
        outf = open(out_fn, 'wt')
        outf.write('H\tVN:Z:1.1\tRS:Z:ref\n')
        for nod in self.nodes:
            outf.write('S\t{}\t{}\n'.format(nod, self.nodes[nod]))
        for pnode in self.edges:
            for snode in self.edges[pnode]:
                outf.write('L\t{}\t+\t{}\t+\t0M\n'.format(pnode, snode))
        for pathn in self.paths:
            path = ','.join([str(nn) + '+' for nn in self.paths[pathn]])
            outf.write('P\t{}\t{}\t*\n'.format(pathn, path))

    def writeSnarlsFreq(self, out_fn):
        outf = open(out_fn, 'wt')
        outf.write('start_node\tnext_node\tgroup\tfreq\n')
        oft = '{}\t{}\t{}\t{}\n'
        for pnode in self.snarls_freq:
            for snode in self.snarls_freq[pnode]:
                freqs = self.snarls_freq[pnode][snode]
                for group in range(len(freqs)):
                    freq = round(freqs[group], 3)
                    outf.write(oft.format(pnode, snode, group, freq))

    def writePhenotype(self, out_ph):
        outf = open(out_ph, 'wt')
        outf.write('FID\tIID\tPHENO\n')
        for idx, samp in enumerate(self.phenotypes):
            # write sample FID, ID and phenotype
            outf.write(f"samp_{idx}\tsamp_{idx}\t{samp}\n")

if "__main__" == __name__ :

    # parse command line arguments
    parser = argparse.ArgumentParser(description='Simulate a simple pangenome graph with SNPs and indels.')
    parser.add_argument('-v', '--nvar', type=int, default=1000, help='Number of top-level variants (default: 1000)')
    parser.add_argument('-n', '--nsamp', type=int, default=200, help='Number of samples (default: 100)')
    parser.add_argument('--snp_prop', type=float, default=0.7, help='Proportion of top-level variants that are SNPs (default: 0.7)')
    parser.add_argument('--nested_prop', type=float, default=0.8, help='Proportion of indels that we want to try to add nested SNPs into (default: 0.8)')
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-b', '--binary', action='store_true', help='Binary phenotypre for each sample')
    group.add_argument('-q', '--quantitative', action='store_true', help='Quantitative phenotype for each sample')
    args = parser.parse_args()

    # number of top-level variants
    nvar = 1000 if args.nvar is None else args.nvar
    # proportion of top-level variants that are SNPs also probability to add or keep adding SNPs when we want to add nested SNPs
    snp_prop = .7 if args.snp_prop is None else args.snp_prop
    # proportion of indels that we want to try to add nested SNPs into
    nested_prop = .8 if args.nested_prop is None else args.nested_prop
    # number of samples
    nsamp = 200 if args.nsamp is None else args.nsamp
    # random seed for reproducibility
    seed = 42
    random.seed(seed)

    _nuc = ["A", "T", "C", "G"]

    # write phenotypes
    if args.binary:
        pg = Graph("binary")
        pg.binaryPhenotype(nsamp)
        pg.writePhenotype('pg.phenotypes.binary.tsv')

    elif args.quantitative:
        pg = Graph("quantitative")
        pg.quantitativePhenotype(nsamp)
        pg.writePhenotype('pg.phenotypes.quantitative.tsv')

    # first node larger than read length
    pnod = pg.addNode(min_size=300, max_size=500)

    # add each variant
    for ii in range(nvar):
        if random.random() < snp_prop:
            # add a SNP
            pnod = pg.addSNP(pnod)
        else:
            # add an indel
            if random.random() < nested_prop:
                # with potentially some nested SNPs
                pnod = pg.addIndel(pnod, snp_prop)
            else:
                # no nested variants
                pnod = pg.addIndel(pnod, 0)

    # add last node larger than read length
    pnod = pg.addNode([pnod], min_size=300, max_size=500)

    # traverse the graph to create a reference path
    pg.addPath()

    # write simple GFA with just the pangenome and reference path
    pg.writeGfa('pg.gfa')

    # Shift the frequencies/probabilities at a subset of variant sites.
    pg.createMarkerFreq(prop_markers=.1, min_dev=.01, max_dev=.5)
    # write the "truth" frequencies
    pg.writeSnarlsFreq('pg.snarls.freq.tsv')

    # traverse the graph to create samples for each group
    # two haplotypes for each sample
    # one sample per group
    for samp in range(nsamp):
        pg.addPath('samp_g0_' + str(samp) + '_h0', group=0, samp=samp)
        pg.addPath('samp_g0_' + str(samp) + '_h1', group=0, samp=samp)
        pg.addPath('samp_g1_' + str(samp) + '_h0', group=1, samp=samp)
        pg.addPath('samp_g1_' + str(samp) + '_h1', group=1, samp=samp)

    # write GFA containing all those haplotype paths
    pg.writeGfa('pg.full.gfa')

# simulate_test_dataset.py
