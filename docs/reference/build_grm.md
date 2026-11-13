# Build a genomic relationship matrix

Reads numeric, HapMap, or VCF genotypes; applies sample and marker QC;
optionally imputes missing calls; builds and checks a genomic
relationship matrix; and writes auditable result files.

## Usage

``` r
build_grm(
  geno_file,
  output_dir,
  input_format = "auto",
  run_label = "GRM",
  marker_col = "SNP",
  chrom_col = "CHR",
  position_col = "POS",
  ref_col = "REF",
  alt_col = "ALT",
  meta_cols = c("CHR", "POS", "REF", "ALT"),
  marker_weight_col = NULL,
  marker_weights_file = NULL,
  marker_weights_id_col = "MarkerID",
  marker_weights_value_col = "Weight",
  multiallelic_action = "remove",
  grm_method = "VanRaden",
  qc_call_rate = 0.95,
  qc_maf = 0,
  qc_sweep_sample = 1,
  qc_heterozygosity = 1,
  imputation_method = "wright",
  qc_seed = 123,
  rrblup_em_tol = 0.02,
  rrblup_cores = 1,
  beagle_jar = NULL,
  java_path = "java",
  beagle_memory = "4g",
  beagle_threads = 1,
  beagle_seed = 123,
  beagle_ne = NULL,
  beagle_ref = NULL,
  beagle_map = NULL,
  beagle_extra_args = character(0),
  write_imputed_genotypes = FALSE,
  n_clusters = NULL,
  cluster_n_pcs = 2,
  cluster_seed = 123,
  overwrite_outputs = FALSE
)
```

## Arguments

- geno_file:

  Path to the genotype file.

- output_dir:

  Directory for result files.

- input_format:

  One of "auto", "numeric", "hapmap", or "vcf".

- run_label:

  Prefix used for every output file.

- marker_col, chrom_col, position_col, ref_col, alt_col:

  Column names for numeric input.

- meta_cols:

  Additional non-genotype columns in numeric input.

- marker_weight_col:

  Optional inline marker-weight column for numeric input.

- marker_weights_file:

  Optional external marker-weight table.

- marker_weights_id_col, marker_weights_value_col:

  Columns in the external marker-weight table.

- multiallelic_action:

  Either "remove" or "error".

- grm_method:

  Either "VanRaden" or "Su".

- qc_call_rate:

  Minimum marker call rate.

- qc_maf:

  Minimum marker minor allele frequency.

- qc_sweep_sample:

  Maximum missing fraction per genotype.

- qc_heterozygosity:

  Maximum observed marker heterozygosity.

- imputation_method:

  One of "none", "wright", "random_hwe", "rrblup_em", or "beagle". With
  "none", missing calls are passed to AGHmatrix and do not cause GSPip
  to stop.

- qc_seed:

  Seed for stochastic imputation.

- rrblup_em_tol, rrblup_cores:

  rrBLUP EM controls.

- beagle_jar, java_path, beagle_memory, beagle_threads, beagle_seed:

  Beagle executable controls.

- beagle_ne, beagle_ref, beagle_map, beagle_extra_args:

  Optional Beagle population, reference, map, and additional arguments.

- write_imputed_genotypes:

  Write the final marker matrix when \`TRUE\`.

- n_clusters:

  Optional k-means cluster count.

- cluster_n_pcs:

  Number of leading GRM PCs used for clustering.

- cluster_seed:

  Clustering seed.

- overwrite_outputs:

  Permit replacement of files with the same run label.

## Value

Invisibly, a `gspip_grm` list containing the GRM, retained marker matrix
and metadata, QC tables, optional clusters, and output paths.
