#!/bin/bash
#SBATCH --job-name=extract_ref
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm-%j.out 
#SBATCH --partition=haswell
#SBATCH --mail-user=ainara.jimenez@upf.edu
#SBATCH --mail-type=ALL

##This script intends to add the reference allele, chromosome and position column to GWAS summary statistics
##when it is needed.

#Load necessary packages
module load BCFtools 

#Possible names for REF, CHR and POS columns in header
chr_pattern="^(chr|chromosome|chrom|snp-chr|chr_ID|#chrom)$"
pos_pattern="^(pos|base_pair_location|snppos|snp_pos|bp|position|chr_pos|pos\\(b37\\))$"
ref_pattern="^(ref|reference_allele|ref_allele)$"

#Define objects
##GWAS list
gwas_list=$1 # Añadir la lista de GWAS (hemos de comprobar la columna ref_allele)

##Directories
gwas_dir=$2
vcf_dir=$3




#Define objects
##Define GWAS file
GWAS=/homes/users/ajimenez/scratch/gwas_tfm/GCST90012622_buildGRCh37.tsv.gz

##Extract the base name (without path and extension)
GWAS_BASE=$(basename "$GWAS" .tsv.gz)

##Define the column numbers where rsIDs are stored
ID_COLUMN=2

##Define an intermediate file to store the raw rsIDs
INTER=/homes/users/ajimenez/scratch/gwas_tfm/unique_rsIDs.txt

##Define file that stores rsIDs and reference allele columns
BASE_OUT=/homes/users/ajimenez/scratch/gwas_tfm/id_ref_${GWAS_BASE}

##Define output
output=/homes/users/ajimenez/scratch/gwas_tfm/${GWAS_BASE}_mod.tsv

#Define temp file
TEMP=/homes/users/ajimenez/gwas_tfm/gcat/temp.tsv

#Create an array of chromosomes (1 to 22)
CHROMOSOMES=( $(seq 1 22) )

##Loop to do the whole process per chromosome
for CHR in "${CHROMOSOMES[@]}"; do
  #Name the input file split by chromosome
  IN_CHR=${GWAS}_chr"$CHR".tsv
  #Create the new input file with just 1 chromosome
  gunzip -c "$GWAS" | awk -F '\t' -v chromosome="$CHR" '$1 == chromosome {print > output_file}' output_file="$IN_CHR"
  
  ##Make a file just with the list of unique rsIDs and query with bcftools
  #Select just the column with the rsID
  cut -f $ID_COLUMN $IN_CHR | \
  #Make just a unique list of rsID
  sort -u > $INTER
  #Name the appropriate dbSNP dataset to query from
  dbSNP=/gpfs/projects/lab_dcomas/1000genomes_phase3_dcomas/vcf/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

  #Name the output file
  OUT=${BASE_OUT}_chr${CHR}.tsv
  #Query SNP database for rsIDs, retrieving REF allele
  bcftools query -f '%ID\t%REF\n' -i 'ID=@/homes/users/ajimenez/scratch/gwas_tfm/unique_rsIDs.txt' $dbSNP > $OUT
  
  #Add header to the output file
  sed -i '1i rsID\tref' $OUT
done

##Merge all the datasets into one
#Name the final output file
REF_ID=${BASE_OUT}.tsv

#Set a condition for the header situation
HEADER_WRITTEN=false
# Create an empty file
> "$REF_ID"

# Iterate over the file names and concatenate them
for ((i=1; i<=22; i++))
do
    INPUT_FILE=${BASE_OUT}_chr${i}.tsv
    # Check if it's the first file and write the header
    if [[ "$HEADER_WRITTEN" = false ]]; then
        cat "$INPUT_FILE" >> "$REF_ID"
        HEADER_WRITTEN=true
    else
        # Skip the header and append the contents
        tail -n +2 "$INPUT_FILE" >> "$REF_ID"
    fi
done

##Clean-up
rm $INTER
rm ${BASE_OUT}_chr*
rm ${IN}_chr*

#Descompress GWAS file
GWAS_desc=/homes/users/ajimenez/scratch/gwas_tfm/${GWAS_BASE}.tsv
gunzip -c "$GWAS" > "$GWAS_desc"

#Add REF allele information to GWAS data
awk -v OFS="\t" '
    BEGIN {
        #Read REF allele mappings into an associative array
        while ((getline < "'$REF_ID'") > 0) {
            variant_map[$1] = $2
        }
    }
    NR == 1 {
        #Print original header and add "ref" column
        print $0, "ref"
    }
    NR > 1 {
        #Print the row with REF if mapping exists, otherwise "NA"
        if ($2 in variant_map) {
            print $0, variant_map[$2]
        } else {
            print $0, "NA"  # Colocar "NA" si no hay mapeo
        }
    }
' $GWAS_desc > $output

##Clean-up
rm $GWAS_desc
#rm $REF_ID




#!/bin/bash

# Definir patrones de búsqueda
chr_pattern="^(chr|chromosome|chrom|snp-chr|chr_ID|#chrom)$"
pos_pattern="^(pos|base_pair_location|snppos|snp_pos|bp|position|chr_pos|pos\\(b37\\))$"
ref_pattern="^(ref|reference_allele|ref_allele)$"

# Directorios donde se encuentran los archivos
gwas_dir="/homes/users/ajimenez/scratch/gwas_tfm"
vcf_dir="/gpfs/projects/lab_dcomas/1000genomes_phase3_dcomas/vcf"

# Archivo con la lista de GWAS y parámetros
input_file="gwas_list.txt"

# Leer cada línea del archivo GWAS list
while IFS=$'\t' read -r gwas_file param; do
    [[ -z "$gwas_file" || "$gwas_file" == "GWAS" ]] && continue
    
    echo "Procesando: $gwas_file con parámetro $param"

    gwas_path="$gwas_dir/$gwas_file"
    gwas_base=$(basename "$gwas_file" .tsv.gz)
    
    # Detectar si el archivo está comprimido
    is_compressed=0
    if [[ "$gwas_file" == *.gz ]]; then
        is_compressed=1
    fi

    # Detectar extensión del archivo
    file_extension="${gwas_file##*.}"
    delimiter="\t"  # Predeterminado para .vcf y .txt

    if [[ "$file_extension" == "csv" || "$file_extension" == "tsv" ]]; then
        delimiter=","
    fi

    # Leer el encabezado del archivo
    if [[ "$is_compressed" -eq 1 ]]; then
        header=$(zcat "$gwas_path" | head -n 1)
    else
        header=$(head -n 1 "$gwas_path")
    fi

    # Buscar si existen las columnas CHR, POS y REF
    chr_match=$(echo "$header" | tr "$delimiter" '\n' | grep -iE "$chr_pattern")
    pos_match=$(echo "$header" | tr "$delimiter" '\n' | grep -iE "$pos_pattern")
    ref_match=$(echo "$header" | tr "$delimiter" '\n' | grep -iE "$ref_pattern")

    # Determinar si falta CHR, POS o REF
    missing_chr_pos=$([[ -z "$chr_match" || -z "$pos_match" ]] && echo "YES" || echo "NO")
    missing_ref=$([[ -z "$ref_match" ]] && echo "YES" || echo "NO")

    # Si falta CHR/POS o el parámetro es NULL, integrar desde VCF
    if [[ "$param" == "NULL" || "$missing_chr_pos" == "YES" ]]; then
        echo "Faltan CHR/POS o el parámetro es NULL. Integrando desde VCF..."

        # Procesar archivo GWAS por cromosoma
        for CHR in $(seq 1 22); do
            # Nombre de archivo de entrada dividido por cromosomas
            IN_CHR="${gwas_base}_chr${CHR}.tsv"
            OUT_CHR="${gwas_base}_chr${CHR}_modified.tsv"
            
            # Crear un archivo con solo una cromosoma
            if [[ "$is_compressed" -eq 1 ]]; then
                zcat "$gwas_path" | awk -F "$delimiter" -v chromosome="$CHR" '$1 == chromosome {print > "'$IN_CHR'"}'
            else
                awk -F "$delimiter" -v chromosome="$CHR" '$1 == chromosome {print > "'$IN_CHR'"}' "$gwas_path"
            fi

            # Extraer una lista única de rsIDs
            cut -f 2 "$IN_CHR" | sort -u > "${gwas_base}_rsIDs_chr${CHR}.txt"
            
            # Definir el archivo VCF correspondiente
            dbSNP="$vcf_dir/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

            # Consultar el archivo VCF para obtener los alelos de referencia
            bcftools query -f '%ID\t%REF\n' -i "ID=@${gwas_base}_rsIDs_chr${CHR}.txt" "$dbSNP" > "${gwas_base}_ref_chr${CHR}.tsv"

            # Añadir el encabezado con la columna "ref"
            sed -i '1i rsID\tref' "${gwas_base}_ref_chr${CHR}.tsv"

            # Integrar el alelo de referencia en el archivo GWAS
            awk -v OFS="$delimiter" '
                BEGIN {
                    # Leer mapeo de REF en un array asociativo
                    while ((getline < "'${gwas_base}_ref_chr${CHR}.tsv'") > 0) {
                        variant_map[$1] = $2
                    }
                }
                NR == 1 {
                    # Imprimir el encabezado original y agregar "ref"
                    print $0, "ref"
                }
                NR > 1 {
                    # Imprimir fila con REF si existe el mapeo, sino "NA"
                    if ($2 in variant_map) {
                        print $0, variant_map[$2]
                    } else {
                        print $0, "NA"
                    }
                }
            ' "$IN_CHR" > "$OUT_CHR"

            # Eliminar archivos temporales
            rm "$IN_CHR" "${gwas_base}_rsIDs_chr${CHR}.txt" "${gwas_base}_ref_chr${CHR}.tsv"
        done
        
        # Combinar todos los archivos de cromosomas en un solo archivo
        cat "${gwas_base}_chr"*"_modified.tsv" > "${gwas_base}_modified_combined.tsv"
        
        echo "Integración completada para $gwas_file: ${gwas_base}_modified_combined.tsv"
    else
        echo "No es necesario modificar $gwas_file"
    fi

done < "$input_file"
