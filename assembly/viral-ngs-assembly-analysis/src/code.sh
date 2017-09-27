#!/bin/bash

main() {
    set -e -x -o pipefail

    if [ -z "$name" ]; then
        name="${assembly_prefix%.refined.refined}"
    fi

    pids=()
    dx cat "$resources" | pigz -dc | tar x -C / & pids+=($!)
    dx download "$assembly" -o assembly.fasta & pids+=($!)
    dx download "$reads" -o reads.bam & pids+=($!)
    mkdir gatk/
    dx cat "$gatk_tarball" | tar jx -C gatk/
    for pid in "${pids[@]}"; do wait $pid || exit $?; done

    if [ "$novocraft_license" != "" ]; then
        dx cat "$novocraft_license" > novoalign.lic
    fi

    # index assembly
    viral-ngs bash -c "read_utils.py index_fasta_picard /user-data/assembly.fasta &&
                       read_utils.py index_fasta_samtools /user-data/assembly.fasta &&
                       novoindex /user-data/assembly.nix /user-data/assembly.fasta"

    # align reads, dedup, realign, filter
    viral-ngs read_utils.py align_and_fix /user-data/reads.bam /user-data/assembly.fasta \
        --outBamAll /user-data/all.bam --outBamFiltered /user-data/mapped.bam \
        --GATK_PATH /user-data/gatk \
        --aligner_options "$aligner_options" --NOVOALIGN_LICENSE_PATH /user-data/novoalign.lic
    samtools index mapped.bam

    # collect some statistics
    assembly_length=$(tail -n +1 assembly.fasta | tr -d '\n' | wc -c)
    alignment_read_count=$(samtools view -c mapped.bam)
    reads_paired_count=$(samtools flagstat all.bam | grep properly | awk '{print $1}')
    alignment_base_count=$(samtools view mapped.bam | cut -f10 | tr -d '\n' | wc -c)
    mean_coverage_depth=$(( alignment_base_count / assembly_length ))
    samtools flagstat  all.bam > stats.txt

    # only plot coverage if input bam has reads
    if [ $alignment_read_count -gt 0 ]; then
      viral-ngs reports.py plot_coverage /user-data/mapped.bam /user-data/coverage_plot.pdf --plotFormat pdf --plotWidth 1100 --plotHeight 850 --plotDPI 100
    else
      echo "No reads present in mapped bam file; skipping plot generation."
    fi

    # Continue gathering statistics
    genomecov=$(bedtools genomecov -ibam mapped.bam | dx upload -o "${name}.genomecov.txt" --brief -)

    # upload outputs
    dx-jobutil-add-output assembly_length $assembly_length
    dx-jobutil-add-output reads_paired_count $reads_paired_count
    dx-jobutil-add-output alignment_read_count $alignment_read_count
    dx-jobutil-add-output alignment_base_count $alignment_base_count
    dx-jobutil-add-output mean_coverage_depth $mean_coverage_depth
    dxid="$(dx upload all.bam --destination "${name}.all.bam" --brief)"
    dx-jobutil-add-output all_reads --class=file "$dxid"
    dxd="$(dx upload stats.txt --destination "${name}.flagstat.txt" --brief)"
    dx-jobutil-add-output bam_stat --class=file "$dxid"
    dxid="$(dx upload mapped.bam --destination "${name}.mapped.bam" --brief)"
    dx-jobutil-add-output assembly_read_alignments --class=file "$dxid"
    dxid="$(dx upload mapped.bam.bai --destination "${name}.mapped.bam.bai" --brief)" 
    dx-jobutil-add-output assembly_read_index --class=file "$dxid"
    dx-jobutil-add-output alignment_genomecov "$genomecov"
    dxid="$(dx upload assembly.fasta --destination "${name}.fasta" --brief)"
    dx-jobutil-add-output final_assembly --class=file "$dxid"
    if [ $alignment_read_count -gt 0 ]; then
      dxid="$(dx upload coverage_plot.pdf --destination "${name}.coverage_plot.pdf" --brief)"
      dx-jobutil-add-output coverage_plot --class=file "$dxid"
    else
      echo "No reads present in mapped bam file; skipping plot upload."
    fi
}
