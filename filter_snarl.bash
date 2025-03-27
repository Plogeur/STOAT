#!/bin/bash

# Vérification des arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <snarl_paths_file> <reference_file>"
    exit 1
fi

snarl_paths_file="$1"
reference_file="$2"

# Extraction de la colonne SNARL du fichier de référence
awk 'NR>1 {print $2}' "$reference_file" | sort > reference_snarls.tmp

echo "chr	pos	snarl	paths	type" > filtered_snarl_paths.txt

# Filtrage du fichier snarl_paths
awk 'NR==1 || FNR==NR { ref[$1]; next } $3 in ref' reference_snarls.tmp "$snarl_paths_file" >> filtered_snarl_paths.txt

# Nettoyage du fichier temporaire
rm reference_snarls.tmp

echo "Fichier filtré sauvegardé sous filtered_snarl_paths.txt"