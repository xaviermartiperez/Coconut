#!/bin/bash
#SBATCH --job-name=merge_Ref
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm-%j.out 
#SBATCH --partition=normal
#SBATCH --mail-user=ainara.jimenez@upf.edu
#SBATCH --mail-type=ALL

GWAS=/homes/users/ajimenez/scratch/gwas_tfm/GCST90092809_buildGRCh37.tsv.gz
GWAS_desc=/homes/users/ajimenez/scratch/gwas_tfm/GCST90092809_buildGRCh37.tsv
REF_ID=/homes/users/ajimenez/scratch/gwas_tfm/id_ref_GCST90092809.tsv
output=/homes/users/ajimenez/scratch/gwas_tfm/GCST90092809_buildGRCh37_mod.tsv

gunzip -c "$GWAS" > $GWAS_desc

#awk 'NR==FNR{a[$1]=$2; next} ($2 in a){print $0, a[$2]}' $REF_ID $GWAS_desc > $output
awk -v OFS="\t" '
    BEGIN {
        # Leer el archivo de referencia en un array asociativo
        while ((getline < "'$REF_ID'") > 0) {
            chr_map[$3] = $1  # Columna 1 -> chr
            pos_map[$3] = $2  # Columna 2 -> pos
            ref_map[$3] = $4  # Columna 4 -> ref
        }
    }
    NR == 1 {
        # Imprimir el encabezado original y agregar "chr", "pos", y "ref"
        print $0, "chr", "pos", "ref"
    }
    NR > 1 {
        # Verificar si el rsID (columna 3 del GWAS) está en el mapa
        if ($1 in chr_map) {
            print $0, chr_map[$1], pos_map[$1], ref_map[$1]
        } else {
            # Colocar "NA" si no hay mapeo
            print $0, "NA", "NA", "NA"
        }
    }
' $GWAS_desc > $output

rm $GWAS_desc
#rm $REF_ID