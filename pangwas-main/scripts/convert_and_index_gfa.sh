set -ex

# This is used to convert gfa to gbz, set the reference, and get the distance index
# Setting the reference may not work well if the sample names of paths aren't specified
# See https://github.com/vgteam/vg/wiki/Changing-References 

GRAPH_BASE=$1
OUT_GRAPH_BASE=$2
REF_SAMPLE=$3

if [[ ! -f ${GRAPH_BASE}.gbz ]]
then
    # Convert gfa to gbz
    vg gbwt --gbz-format -g ${OUT_GRAPH_BASE}.temp.gbz -G ${GRAPH_BASE}.gfa
    vg gbwt --set-tag "reference_samples=${REF_SAMPLE}" --gbz-format -g ${OUT_GRAPH_BASE}.gbz -Z ${OUT_GRAPH_BASE}.temp.gbz
    rm ${OUT_GRAPH_BASE}.temp.gbz
fi

if [[ ! -f ${OUT_GRAPH_BASE}.dist ]]
then
    # Distance index chopped gbz
    vg index -j ${OUT_GRAPH_BASE}.dist ${OUT_GRAPH_BASE}.gbz
fi


# This converts to hg format, which may not work well for gbz's with a lot of paths

# If there is no reference path, try to promote the given reference sample to a reference

# If there are no reference paths
#if [[ $(vg paths -L -R -x ${OUT_GRAPH_BASE}.gbz | wc -l) -eq 0 ]]
#then
#    # If the reference we want is not a haplotype
#    if [[ $(vg paths -L -H -x ${OUT_GRAPH_BASE}.gbz | grep ${REF_SAMPLE} | wc -l) -eq 0 ]]
#    then
#        # Convert the generic path with the reference locus to a haplotype path with the reference as its sample name
#        vg convert --hap-locus ${REF_SAMPLE} --new-sample ${REF_SAMPLE} -a ${OUT_GRAPH_BASE}.hg > ${OUT_GRAPH_BASE}.hap.hg
#        mv ${OUT_GRAPH_BASE}.hap.hg ${OUT_GRAPH_BASE}.hg
#    fi
#    # Convert the haplotype path with the reference as its sample name into a reference
#    vg convert --ref-sample ${REF_SAMPLE} -a ${OUT_GRAPH_BASE}.hg > ${OUT_GRAPH_BASE}.ref.hg
#    mv ${OUT_GRAPH_BASE}.ref.hg ${OUT_GRAPH_BASE}.hg
#fi

