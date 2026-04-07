#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//-- Help Message ---------------------------------------------------------------

def helpMessage() {
  log.info """
  ===============================
  BUSCO Phylogenomics Pipeline
  Nextflow pipeline for BUSCO-based multigene phylogenomics
  Version: 0.1.2
  Author : Akito Shima (ASUQ)
  Email: akito.shima@oist.jp
  ===============================
  Usage: nextflow run main.nf [parameters]

  Required parameters:
    --sample      Path to sample.csv (header: sample,fasta)
    --lineage     BUSCO lineage dataset (e.g. mycoplasmatota_odb12)

  Optional parameters:
    --help        Show this help message
    --outdir      Output directory (default: ./output)
    --fraction    Comma-separated completeness fractions (default: 0.8,0.9,1.0)
    --busco_opts  Extra BUSCO flags (default: "")
    --mafft_opts  MAFFT options (default: --globalpair --maxiterate 1000)
    --trimal_opts trimAl options (default: -automated1)
    --amas_opts   AMAS concat options (default: --in-format fasta --data-type aa --part-format nexus --cores \$task.cpus)
    --iqtree_opts IQ-TREE options (default: -B 1000 -alrt 1000 -m MFP+MERGE -T \$task.cpus)
  """.stripIndent()
}


def missingParametersError() {
    log.error "Missing input parameters"
    helpMessage()
    error "Please provide all required parameters: --sample and --lineage"
}


//-- Processes -----------------------------------------------------------------

// Download BUSCO dataset for offline use
process download_busco_dataset {
    label 'process_single'

    publishDir "${params.outdir}/busco_downloads", mode: 'copy'

    input:
    val lineage

    output:
    path "busco_downloads", emit: lineage_dir

    script:
    """
    busco --download "${lineage}"
    """

    stub:
    """
    echo "Stub process for downloading BUSCO dataset: ${lineage}"
    mkdir -p "busco_downloads/lineages/${lineage}"
    """
}

// Run BUSCO for each sample in offline mode
process busco {
    label 'process_high_memory'
    tag   "${sample}"

    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    tuple val(sample), path(fasta)
    path lineage_dir
    val busco_opts

    output:
    path "${sample}", emit: busco_dir

    script:
    """
    busco --in "${fasta}" \
          --lineage_dataset "${params.lineage}" \
          --out "${sample}/busco_output" \
          --mode genome \
          --cpu ${task.cpus} \
          --offline \
          ${busco_opts}
    """

    stub:
    """
    echo "Stub process for BUSCO: ${sample}"
    mkdir -p "${sample}/busco_output"
    """
}

// Collect per-gene FASTA files from BUSCO outputs
process collect_and_select_genes {
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sample_busco, stageAs: 'busco/*'
    val fractions

    output:
    path 'seqs', emit: seqs_dir
    path "frac*pct_results", emit: frac_results

    script:
    """
    "${projectDir}/bin/busco_multigene_tree.py" collect \
      --input_dir 'busco' --out_dir '.' --cores ${task.cpus}

    "${projectDir}/bin/busco_multigene_tree.py" select \
      --input_dir 'busco' --out_dir '.' \
      --fraction ${fractions} --cores ${task.cpus}
    """

    stub:
    """
    echo "Stub process for collecting and selecting BUSCO genes"
    mkdir -p seqs/raw
    touch seqs/raw/geneA.faa
    """
}

// Align and trim genes
process align_genes {
    label 'process_high'
    tag   "${gene}"

    publishDir "${params.outdir}/seqs/aligned", pattern: '*_aligned.faa', mode: 'copy'
    publishDir "${params.outdir}/seqs/trimmed", pattern: '*_trimmed.faa', mode: 'copy'

    input:
    tuple val(gene), path(faa)

    output:
    path "${gene}_aligned.faa", emit: aligned
    path "${gene}_trimmed.faa", emit: trimmed

    script:
    def mafft_opts  = (params.mafft_opts  ?: '')
    def trimal_opts = (params.trimal_opts ?: '')
    """
    mafft ${mafft_opts} --thread ${task.cpus} --threadtb ${task.cpus} \
      --threadit ${task.cpus} "${faa}" > "${gene}_aligned.faa"

    trimal ${trimal_opts} \
      -in  "${gene}_aligned.faa" \
      -out "${gene}_trimmed.faa"
    """

    stub:
    """
    echo "Stub process for aligning and trimming gene: ${gene}"
    touch "${gene}_aligned.faa"
    touch "${gene}_trimmed.faa"
    """
}

// Concatenate alignments and infer phylogenetic trees
process infer_trees {
    label 'process_medium'
    tag   "frac${pct}pct"

    publishDir "${params.outdir}/frac${pct}pct_results", mode: 'copy'

    input:
    tuple val(pct), path(gene_list)
    path  trimmed_all                 // all *_trimmed.faa staged as inputs

    output:
    path "concat.faa", optional: true
    path "partitions.nex", optional: true
    path "frac*", optional: true

    script:
    def amas_opts   = (params.amas_opts   ?: '')
    def iqtree_opts = (params.iqtree_opts ?: '')
    """
    readarray -t GENES < <(tail -n +3 "${gene_list}" | sed -e 's/\r\$//' -e '/^[[:space:]]*\$/d')

    if [ "\${#GENES[@]}" -eq 0 ]; then
      echo "[WARNING] No genes found in the fraction file: ${gene_list}" >&2
      exit 0
    fi

    TRIMMED_FILES=""
    for gene in "\${GENES[@]}"; do
      f="\${gene}_trimmed.faa"
      [[ -f "\$f" ]] || { echo "[ERROR] Missing trimmed file: \$f" >&2; exit 2; }
      TRIMMED_FILES+="\$f "
    done

    AMAS.py concat ${amas_opts} \
      --cores ${task.cpus} \
      --concat-out  concat.faa \
      --concat-part partitions.nex \
      --in-files \${TRIMMED_FILES}

    iqtree ${iqtree_opts} \
      -T   ${task.cpus} \
      -s   concat.faa \
      -p   partitions.nex \
      -pre "frac${pct}pct"
    """

    stub:
    """
    echo "Stub process for inferring trees for fraction ${pct}%"
    mkdir -p frac${pct}pct_results
    touch "frac${pct}pct_results/concat.faa"
    touch "frac${pct}pct_results/partitions.nex"
    touch "frac${pct}pct_results/frac${pct}pct.best_model.nex"
    touch "frac${pct}pct_results/frac${pct}pct.best_scheme"
    touch "frac${pct}pct_results/frac${pct}pct.best_scheme.nex"
    touch "frac${pct}pct_results/frac${pct}pct.bionj"
    touch "frac${pct}pct_results/frac${pct}pct.ckp.gz"
    touch "frac${pct}pct_results/frac${pct}pct.contree"
    touch "frac${pct}pct_results/frac${pct}pct_genes.txt"
    touch "frac${pct}pct_results/frac${pct}pct.iqtree"
    touch "frac${pct}pct_results/frac${pct}pct.log"
    touch "frac${pct}pct_results/frac${pct}pct.mldist"
    touch "frac${pct}pct_results/frac${pct}pct.model.gz"
    touch "frac${pct}pct_results/frac${pct}pct.splits.nex"
    touch "frac${pct}pct_results/frac${pct}pct.treefile"
    """
}


//-- Workflow ------------------------------------------------------------------
workflow {
  // Parameter parsing
  if (params.help) {
    helpMessage()
    exit 0
  }

  if (params.sample == null || params.lineage == null) {
    missingParametersError()
    exit 1
  }

  // Read fractions (e.g. "0.8,0.9,1.0"), pick smallest for alignment
  frac_list = (params.fraction as String).tokenize(',').collect { it as Double }.sort()
  smallest_pct = (frac_list[0] * 100) as int

  // Channel setup
  fasta_ch = Channel.fromPath(params.sample, checkIfExists: true)
                    .splitCsv(strip: true, header: true)

  // Download BUSCO lineage dataset
  busco_db = download_busco_dataset(params.lineage)

  // Run BUSCO for each sample
  busco_results = busco(fasta_ch, busco_db, params.busco_opts)
                      .collect()

  // Collect and select genes from BUSCO results
  busco_genes = collect_and_select_genes(busco_results, params.fraction)

  // Build gene channel for the smallest fraction only -> per-gene alignment input
  min_frac_gene_ch = busco_genes.frac_results
                                .flatten()
                                .filter { it.name == "frac${smallest_pct}pct_results" }
                                .map { dir -> file("${dir}/frac${smallest_pct}pct_genes.txt") }
                                .splitText()
                                .map { it.trim() }                   // remove newline (\n)
                                .filter {
                                  it &&                              // drop blanks
                                  !it.startsWith('Number of genes considered') &&
                                  !it.startsWith('Analyzed genes')   // drop the 2 header lines
                                }


  // raw dir path from 'seqs' output
  raw_dir_ch = busco_genes.seqs_dir
                          .map { dir -> file("${dir}/raw") }

  // Create (gene, fasta) tuples by combining gene list + raw dir
  gene_ch = min_frac_gene_ch.combine(raw_dir_ch)
                            .map { gene, rawDir -> tuple(gene, file("${rawDir}/${gene}.faa")) }

  // Align & trim each gene
  aligned = align_genes(gene_ch)

  // Gather all trimmed files to feed each infer task
  trimmed_all_ch = aligned.trimmed.collect()

  // Build (pct_int, gene_list_file) tuples for ALL fractions found
  frac_gene_file_ch = busco_genes.frac_results
                                .flatten()
                                .map { dir ->
                                    def m = (dir.name =~ /frac(\d+)pct_results/)
                                    assert m : "Unexpected directory name: ${dir.name}"
                                    def pct = (m[0][1] as int)
                                    tuple(pct, file("${dir}/frac${pct}pct_genes.txt"))
                                }

  // Run one infer job per integer percent in parallel
  infer_trees(frac_gene_file_ch, trimmed_all_ch)
}
