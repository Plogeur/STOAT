# Report comparison stoat vs plink

## Intro
Pour tester les différentes  outils nous réalisons 2 type de simulation : binaire et quantitative.

La simulation binaire contient 100 échantillons dont chacun d’eux présentent 1000 variations, chaque échantillon est simulé pour appartenir à un dès 2 groupes. Les probabilité de passage d’un nœud à un autre dans le graph est connues pour chacun des groupes.

La simulation quantitative contient 200 échantillons pour 1000 variantions, chaque échantillon est simulé avec un phénotype quantitative normalizer en 0, tout comme la simulation binaire les échantillons sont soumis à une appartenance à 2 groupes mais cette dépendante est moduler ici la valeur du phénotype de chaque échantillons. Les probabilités de passage d’un nœud à un autre dans le graph est connues pour chacun des groupes.

Fichiers présents : 
- VCFs (sortie du pipeline de calling)
- Probabilité de passage des noeuds
- Information sur le graph pangénomique (.dist, .pg, .gfa)
- Phenotype

## Méthodes
Dans un premier temps nous voulons tester les différents outils sans un preprocessing des données au préalable, pour cela nous mergons les VCFs avec la commande bcftools merge, puis nous lançons les différents outils (stoat/plink) sur le VCF mergé obtenue.  

```bash
bcftools merge $dir_vcf/*.vcf.gz --threads 6 -m none -Oz -o $output/merged_vcf.vcf
```

Puis nous avons décidé de réaliser une décomposition des VCFs avec un pipeline en snakemake afin d’enlever les snarl imbriquer qui ne serait pas bien analyser par plink. Dans cette étapes les VCFs seronts normalizer puis merger.

```bash
snakemake --cores 8 --snakefile Snakefile --configfile config/config.yaml
```

## Result
### stoat
#### Binaire 




#### Quantitatif 





### Plink 
#### Binaire 
#### Quantitatif 

Interprétation
