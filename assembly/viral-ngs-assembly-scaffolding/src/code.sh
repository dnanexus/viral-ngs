#!/bin/bash

main() {
    set -e -x -o pipefail

    pids=()
    dx cat "$resources" | pigz -dc | tar x -C / & pids+=($!)
    dx download "$trinity_contigs" -o trinity_contigs.fasta & pids+=($!)
    dx download "$reference_genome" -o reference_genome.fasta & pids+=($!)
    dx download "$trinity_reads" -o reads.bam
    for pid in "${pids[@]}"; do wait $pid || exit $?; done

    # run assembly.py order_and_orient to scaffold the contigs
    viral-ngs assembly.py order_and_orient \
        /user-data/trinity_contigs.fasta /user-data/reference_genome.fasta /user-data/intermediate_scaffold.fasta

    if [ -z "$name" ]; then
        name=${trinity_contigs_prefix%_1}
    fi

    # run assembly.py impute_from_reference to check assembly quality and clean the contigs
    exit_code=0
    touch impute.stderr.log
    viral-ngs assembly.py impute_from_reference \
        /user-data/intermediate_scaffold.fasta /user-data/reference_genome.fasta /user-data/scaffold.fasta
        --newName "${name}" --replaceLength "$replace_length" \
        --minLengthFraction "$min_length_fraction" --minUnambig "$min_unambig" \
        --aligner "$aligner" 2> >(tee impute.stderr.log >&2) || exit_code=$?

    if [ "$exit_code" -ne "0" ]; then
        if grep PoorAssemblyError impute.stderr.log ; then
            dx-jobutil-report-error "The assembly failed quality thresholds (length fraction >= ${min_length_fraction}, non-N proportion >= ${min_unambig})" AppError
        else
            dx-jobutil-report-error "Please check the job log" AppInternalError
        fi
        exit $exit_code
    fi

    test -s scaffold.fasta

    # upload outputs
    dxid=$(dx upload scaffold.fasta --destination "${name}.scaffold.fasta" --brief)
    dx-jobutil-add-output modified_scaffold --class=file "$dxid"
    dxid=$(dx upload intermediate_scaffold.fasta --destination "${name}.mummer.fasta" --brief)
    dx-jobutil-add-output intermediate_scaffold --class=file "$dxid"
}

first_fasta_header() {
    head -n 1 "$1" | tr -d ">\n"
}
