#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.types = true

params {
    ref: Path
    ref_index: Path
    sample_sheet: Path
    data_dir: Path
    out_prefix: String
    // from sample sheet
    n_samples: Integer
    // set higher for more parallelism across regions
    min_windows: Integer = 1
    primary_assembly: Boolean = true
}

workflow {

    main:
    // collate by size. for some reason round() doesn't return an integer
    size = Math.round(Math.sqrt(params.n_samples)) as Integer
    size = Math.max(size, 2)
    n_batches = Math.ceilDiv(params.n_samples, size)
    log.info("Merging samples in batches of ${size}")

    vcfs = channel.of(params.sample_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def rowMap = row as Map<String, String>
            tuple(
                rowMap['File ID'],
                params.data_dir / rowMap['File ID'] / rowMap['File Name'],
                params.data_dir / rowMap['File ID'] / (rowMap['File Name'] + '.tbi'),
            )
        }
        .take(params.n_samples)

    def prepped: Channel<Tuple<Path, Path>>
    prepped = PREP(vcfs)

    batched = prepped.collate(size, size, true) as Channel<List<Tuple<Path, Path>>>
    batched_uuid = batched.map { ls ->
        tuple(UUID.randomUUID().toString(), ls.collect { tup -> tup[0] }, ls.collect { tup -> tup[1] })
    }

    regions = WINDOW_GENOME(channel.of(params.ref_index), params.min_windows, params.primary_assembly)
        .splitCsv(header: false, sep: '\t')
        .map { row ->
            // chrom, start, end, name, score, strand
            def rowTup = row as Tuple<String, String, String, String, String, String>
            // make 1-based closed intervals
            "${rowTup[0]}:${rowTup[1].toInteger() + 1}-${rowTup[2]}"
        }

    def regions_batched_bcfs: Channel<Tuple<String, String, List<Path>, List<Path>>>
    regions_batched_bcfs = regions.combine(batched_uuid)

    // merge within regions
    def first_merge: Channel<Tuple<String, List<Path>, List<Path>>>
    first_merge = FIRST_MERGE(regions_batched_bcfs).groupTuple(size: n_batches)
    // join on region and merge
    final_merge = FINAL_MERGE(first_merge).toList() as Value<List<Tuple<Path, Path>>>
    transposed_merge = final_merge.map { ls ->
        tuple(ls.collect { p -> p[0] }, ls.collect { p -> p[1] })
    }
    (bcf, index) = CONCAT(transposed_merge, params.out_prefix)

    publish:
    bcf: Path = bcf
    index: Path = index
}

output {
    bcf: Path {
        mode 'copy'
        overwrite true
    }
    index: Path {
        mode 'copy'
        overwrite true
    }
}

process PREP {
    cpus 1
    memory 4.GB

    input:
    (id, vcf, _index): Tuple<String, Path, Path>

    script:
    """
    # Reheader, filter, and normalize
    echo "TUMOR ${id}" > renamer.txt
    
    bcftools view -O u -s TUMOR -f '.,PASS' -i 'GT="alt"' ${vcf} \\
    | bcftools norm -O u -a --atom-overlaps . -m -both --multi-overlaps . -f ${params.ref} \\
    | bcftools annotate -O u -c FORMAT/VAF:=FORMAT/AF \\
    | bcftools annotate -O u -x INFO,FILTER,^FMT/VAF,FMT/GT \\
    | bcftools reheader -s renamer.txt \\
    | bcftools view -W -O b -o ${id}.bcf
    """

    output:
    tuple(file("${id}.bcf"), file("${id}.bcf.csi"))
}

process WINDOW_GENOME {
    cpus 1
    memory 4.GB

    input:
    ref_index: Path
    min_windows: Integer
    primary_assembly: Boolean

    script:
    primary_assembly_flag = primary_assembly ? '--primary-assembly' : ''
    """
    fai_to_regions.py --fai ${ref_index} --min-windows ${min_windows} \\
        ${primary_assembly_flag} --out regions.bed
    """

    output:
    file('regions.bed')
}

process FIRST_MERGE {
    cpus 2 * task.attempt
    memory 16.GB * task.attempt
    errorStrategy task.exitStatus in 137..140 ? 'retry' : 'terminate'
    maxRetries 4
    maxForks 32

    input:
    (region, sample_group_id, bcfs, _indices): Tuple<String, String, List<Path>, List<Path>>

    script:
    bcfs = bcfs.toSorted { p -> p.name }
    if (bcfs.size() == 1) {
        // can't just symlink because we need to subset to the region
        """
        bcftools view --threads ${task.cpus} -r ${region} -O b -o ${sample_group_id}.bcf -W ${bcfs[0]}
        """
    }
    else {
        """
        bcftools merge --threads ${task.cpus} -r ${region} -m none -0 -O b -o ${sample_group_id}.bcf -W ${bcfs.join(' ')}
        """
    }

    output:
    tuple(region, file("${sample_group_id}.bcf"), file("${sample_group_id}.bcf.csi"))
}

process FINAL_MERGE {
    cpus 4 * task.attempt
    memory 32.GB * task.attempt
    errorStrategy task.exitStatus in 137..140 ? 'retry' : 'terminate'
    maxRetries 4
    maxForks 16

    input:
    (region, bcfs, indices): Tuple<String, List<Path>, List<Path>>

    script:
    region_id = region.replaceAll('[:-]', '_')
    bcfs = bcfs.toSorted { p -> p.name }
    if (bcfs.size() == 1) {
        """
        ln -s ${bcfs[0]} ${region_id}.bcf
        ln -s ${indices[0]} ${region_id}.bcf.csi
        """
    }
    else {
        """
        bcftools merge --threads ${task.cpus} -m none -0 -O b -o ${region_id}.bcf -W ${bcfs.join(' ')}
        """
    }

    output:
    tuple(file("${region_id}.bcf"), file("${region_id}.bcf.csi"))
}

process CONCAT {
    cpus 4 * task.attempt
    memory 64.GB * task.attempt
    errorStrategy task.exitStatus in 137..140 ? 'retry' : 'terminate'
    maxRetries 4

    input:
    (bcfs, _indices): Tuple<List<Path>, List<Path>>
    out_prefix: String

    script:
    """
    mapfile -t bcfs < <(echo ${bcfs.join(' ')} | tr ' ' '\\n' | sort -V)
    bcftools concat --threads ${task.cpus} --naive -O b -o ${out_prefix}.bcf \${bcfs[@]}
    bcftools index --threads ${task.cpus} ${out_prefix}.bcf
    """

    output:
    bcf: Path = file("${out_prefix}.bcf")
    index: Path = file("${out_prefix}.bcf.csi")
}
