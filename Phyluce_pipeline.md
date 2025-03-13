```markdown
# UCE Data Pipeline for SI Hydra by Miles Zhang

**Documentation:** [Phyluce Documentation](https://phyluce.readthedocs.io/en/latest/)

**General Note:** If you make an error and create an output, the file will quit and not run. So delete the empty directory and try again.

## Illumiprocessor (Quality Control + Trimming)

- **Make config file** (do not use textedit because of invisible characters causing error)
- **Or use dos2unix** to change your file
- **Conda config -- channels** to see priority
- **Conda config -- remove channels 'x'** to remove
- **Make sure for R1 and R2 tags have {} in front**
- **Run Illumiprocessor:**
  ```bash
  illumiprocessor --input ./raw-fastq --output ./clean-fastq --config ./conf_file --cores 16
  ```

- **Quality Control:**
  ```bash
  # move to the directory holding our cleaned reads
  cd clean-fastq/
  # run this script against all directories of reads
  for i in *; do
      phyluce_assembly_get_fastq_lengths --input $i/split-adapter-quality-trimmed/ --csv;
  done
  ```

## Assembly (ABySS or SPAdes)

- **Subsample if using SPAdes** (as it is very slow and memory intensive)
  - **Calculate # of reads after trimming:**
    ```bash
    for i in clean-reads/*; do echo $i; gunzip -c $i/split-adapter-quality-trimmed/*-READ1.fastq.gz | wc -l | awk '{print $1/4}'; done
    ```
  - **Make a list of samples to downsize** (anything > 1.5 mil if you want 3 mil reads in total, or >2 mil if you want 4 mil reads in total)
  - **Run the assembly:**
    ```bash
    # make a directory for log files
    mkdir log
    # run the assembly
    phyluce_assembly_assemblo_spades \
    --conf assembly.conf \
    --output spades-assemblies \
    --cores 12
    ```

- **Smithsonian Hydra code:**
  ```bash
  nohup ./2.spades_submission.sh ../clean-fastq2-merge ../spades_assemblies2-merge ./spades.job > spades_submission_nohup.out&
  ```

- **Quality Control:**
  ```bash
  for i in ./*.fasta; do
      phyluce_assembly_get_fasta_lengths --input $i --csv;
  done
  ```

## Match to UCE Probe

- **Make sure to download probe set** (hymenopteraV2P.fasta)
- **Run:**
  ```bash
  phyluce_assembly_match_contigs_to_probes \
    --contigs whatever-assemblies/contigs \
    --probes hymenopteraV2P.fasta \
    --output uce-search-results
  ```

## Extracting UCE Loci

- **Config file for taxon set** (name of all species)
- **Run:**
  ```bash
  phyluce_assembly_get_match_counts \
    --locus-db uce-search-results/probe.matches.sqlite \
    --taxon-list-config config/taxon-set.conf \
    --taxon-group 'xxx' \
    --incomplete-matrix \
    --output taxon-sets/xxx-incomplete.conf
  ```

## Get FASTA

- **Run:**
  ```bash
  phyluce_assembly_get_fastas_from_match_counts \
    --contigs ../spades-assemblies/contigs/ \
    --locus-db ../uce-search-results/probe.matches.sqlite \
    --match-count-output xxx-incomplete.conf \
    --output xxx-incomplete.fasta \
    --incomplete-matrix xxx-incomplete.incomplete \
    --log-path log
  ```

## Explode Monolithic FASTA

- **Phasing data?** If so do it here (Andermann et al 2018, Tutorial #2)
- **Run:**
  ```bash
  phyluce_assembly_explode_get_fastas_file --input incomplete_matrix.fasta --output exploded-fastas --by-taxon
  ```
- **Summary stats:**
  ```bash
  for i in exploded-fastas/*.fasta; do
      phyluce_assembly_get_fasta_lengths --input $i --csv;
  done >> Platy_UCE_stats.csv
  ```

## Alignment (MAFFT implemented in Phyluce)

- **Internal trimming, Gblock is needed** (change b1-b4 settings from default)
  ```bash
  cd uce-tutorial/taxon-sets/all
  phyluce_align_seqcap_align \
    --input incomplete_matrix.fasta \
    --output internal \
    --taxa xxx \
    --aligner mafft \
    --cores 16 \
    --incomplete-matrix \
    --output-format fasta \
    --no-trim \
    --log-path log
  ```

- **Run Gblocks:**
  ```bash
  phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments internal \
    --output internal-gblocks \
    --b1 0.5 \
    --b2 0.5 \
    --b3 12 \
    --b4 7 \
    --cores 16 \
    --log log
  ```

- **Align summary data:**
  ```bash
  phyluce_align_get_align_summary_data \
    --alignments internal \
    --input-format fasta \
    --cores 24 \
    --log-path log
  ```

- **Remove locus names/alignment cleaning:**
  ```bash
  phyluce_align_remove_locus_name_from_files \
    --alignments internal-gblocks \
    --output internal-gblocks-clean \
    --cores 24 \
    --log-path log
  ```

## Generate Data Matrices

- **Run:**
  ```bash
  phyluce_align_get_only_loci_with_min_taxa \
    --alignments internal-gblocks-clean \
    --taxa xxx \
    --percent 0.x \
    --output internal-gblocks-clean-x0p \
    --cores 24 \
    --log-path log
  ```

## RAxML

- **Run:**
  ```bash
  phyluce_align_concatenate_alignments \
    --alignments internal-gblocks-clean-x0p \
    --output internal-gblocks-clean-x0p-phylip \
    --phylip \
    --log-path log
  ```

## NEXUS (works for SWSC-EN)

- **Run:**
  ```bash
  phyluce_align_concatenate_alignments \
    --alignments internal-gblocks-clean-x0p \
    --output internal-gblocks-clean-x0p-nexus \
    --nexus \
    --log-path log
  ```

## Extract COI

- **Run:**
  ```bash
  phyluce_assembly_extract_contigs_to_barcodes \
    --contigs ../../spades-results/spades-assemblies/contigs/ \
    --config coi.fas \
    --output barcodes \
    --log-path log
  ```

## PartitionFinder2

- **SWSC-EN:**
  - Use `phyluce_align_format_nexus_files_for_raxml --output nexus` to generate input concatenated nexus file
  - Modify charsets so it shows '#:UCE loci #' instead of the default 'UCE loci #: character lengths'
  - Run program using Anaconda python:
    ```bash
    python /path/to/SWSCEN.py /path/to/nexus
    ```
- **Use output to run PF2**

## IQ-TREE for Concatenated Data

- **Unpartitioned, partitioned by locus, and partitioned by SWSCEN:**
  ```bash
  module load iq-tree
  iqtree -s input.phylip -spp partition.txt -alrt 1000 -bb 1000 -bnni
  ```

- **Concordance Factor -czb to collapse low nodal support for gene trees:**
  ```bash
  iqtree -s ALN_FILE -p PARTITION_FILE --prefix concat -bb 1000 -nt
  iqtree -s ALN_FILE -S PARTITION_FILE --prefix loci -nt
  iqtree -t concat.treefile --gcf loci.treefile -s ALN_FILE --scf 100 --prefix concord -nt 8
  ```

## ASTRAL for Species Tree

- **Generate gene trees with IQ-TREE:**
  ```bash
  module load iq-tree
  iqtree -s input.phylip -bb 1000
  ```

- **Generate species tree with ASTRAL with gene trees as input:**
  ```bash
  java -jar astral.5.7.3.jar -i in.tree -o out.tre
  ```

## BEAST for Divergence Dating

- **Install SA package in your own account** (log in through VPN if off campus):
  ```bash
  Module load gui/2
  gui start --module beast/XXX -e beauti
  Install SA
  ```

## Phasing

- **Run:**
  ```bash
  phyluce_align_explode_alignments \
    --alignments edge \
    --input-format fasta \
    --output edge-exploded \
    --by-taxon
  ```

- **Make phasing.conf file:**
  ```ini
  [references]
  mus_musculus:/Users/bcf/tmp/phyluce/mafft-nexus-edge-trimmed-exploded/mus_musculus.fasta

  [individuals]
  mus_musculus:/Users/bcf/tmp/phyluce/clean-fastq/mus_musculus/split-adapter-quality-trimmed/

  [flowcell]
  alligator_mississippiensis:D1HTMACXX
  ```

- **Realign using BWA:**
  ```bash
  phyluce_snp_bwa_multiple_align \
    --config phasing.conf \
    --output multialign-bams \
    --cores 16 \
    --log-path log \
    --mem
  ```

- **Phase UCEs:**
  ```bash
  phyluce_snp_phase_uces \
    --config phasing.conf \
    --bams multialign-bams \
    --output multialign-bams-phased-reads
  ```

- **Align phased sequences:**
  ```bash
  phyluce_align_seqcap_align \
    --input joined_allele_sequences_all_samples.fasta \
    --output phased-edge-no-ambi \
    --taxa xxx \
    --aligner mafft \
    --cores 16 \
    --incomplete-matrix \
    --log-path log
  ```

- **Run Gblocks on phased data:**
  ```bash
  phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments phased-internal \
    --output phased-internal-gblocks \
    --b1 0.5 \
    --b2 0.5 \
    --b3 12 \
    --b4 7 \
    --cores 16 \
    --log log
  ```

- **Remove locus names in the files:**
  ```bash
  phyluce_align_remove_locus_name_from_files \
    --alignments phased-edge-no-ambi \
    --output phased-edge-no-ambi-clean \
    --cores 16 \
    --log-path log
  ```

- **Align summary data:**
  ```bash
  phyluce_align_get_align_summary_data \
    --alignments phased-edge-no-ambi \
    --cores 16 \
    --log-path log
  ```

- **Generate data matrices:**
  ```bash
  phyluce_align_get_only_loci_with_min_taxa \
    --alignments phased-edge-no-ambi-clean \
    --taxa xxx \
    --percent 0.x \
    --output phased-edge-no-ambi-clean-x0p \
    --cores 16 \
    --log-path log
  ```

- **Concatenate alignments:**
  ```bash
  phyluce_align_concatenate_alignments \
    --alignments phased-edge-no-ambi-clean-x0p \
    --output phased-edge-no-ambi-clean-x0p-phylip \
    --phylip \
    --log-path log
  ```

- **Use phrynomics to change nexus file to binary:**
  ```bash
  add FORMAT DATATYPE=integerdata symbols="012" gap=- MISSING=?;
  ```

## Protein Coding

- **Create protein database:**
  ```bash
  uce_to_protein.py blastdb -i proteins.fasta (download one from here)
  ```

- **Explode unaligned fasta by loci:**
  ```bash
  phyluce_assembly_explode_get_fastas_file --input incomplete_matrix.fasta --output unaligned
  ```

- **Remove locus names in the files:**
  ```bash
  for file in *.fasta; do
      sed -i 's/uce-[0-9]\{1,5\}_\([^|]*\) |uce-[0-9]\{1,5\}/\1/' "$file"
  done
  ```

- **Query BLAST search:**
  ```bash
  python uce_to_protein.py queryblast -i unaligned/*.fasta -c 16
  ```

- **Parse Blast Results:**
  ```bash
  Python uce_to_protein.py parse -i fasta_to_xml.conf -o my_uce_hits.sqlite
  ```

- **Write FASTA:**
  ```bash
  Python uce_to_protein.py queryprot -d my_uce_hits.sqlite -c taxon_config.conf -g taxon_group
  ```

- **Generate data matrices:**
  ```bash
  phyluce_align_get_only_loci_with_min_taxa \
    --alignments protein-aligned-trimal \
    --taxa xxx \
    --percent 0.x \
    --output protein-aligned-trimal-x0p \
    --input-format fasta \
    --cores 24 \
    --log-path log
  ```

- **RAxML:**
  ```bash
  phyluce_align_concatenate_alignments \
    --alignments protein-aligned-trimal-x0p \
    --output protein-aligned-trimal-x0p-phylip \
    --phylip \
    --input-format fasta \
    --log-path log
  ```
```