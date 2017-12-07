#!/bin/bash

# (https://www.gnu.org/software/bash/manual/html_node/The-Set-Builtin.html)
#   - `-e`: Exit immediately if any pipeline returns nonzero status.
#   - `-o pipefail`: Return value of pipeline is value of last command to exit
#     with nonzero status, or zero if all commands succeeded.
#   - `-x`: Print trace of many of commands below; useful for debugging.
set -e -x -o pipefail

# Log system resource statistics every 1m.
dstat -cmdn 60 &

function main() {
  dx cat "$resources" | pigz -dc | tar x -C /

  # Fetch and decompress Kraken & Krona databases.
  extract_db "$kraken_db" "$kraken_db_name" "$kraken_db_prefix"
  extract_db "$krona_taxonomy_db" "$krona_taxonomy_db_name" "$krona_taxonomy_db_prefix"

  # Process input samples
  export -f process_bam
  export SHELL=/bin/bash
  parallel --delay 1 -P 8 -t process_bam "$kraken_db_prefix" "$krona_taxonomy_db_prefix" ::: "${mappings[@]}"

  # upload outputs
  dx-upload-all-outputs --parallel
}

function extract_db() {
  db_id="$1"
  db_name="$2"
  db_prefix="$3"

  mkdir "./$db_prefix"

  # Db can be compressed by gzip or LZ4; choose decompressor based on file
  # extension.
  decompressor="pigz -dc"
  if [[ "$db_name" == *.lz4 ]]; then
    decompressor="lz4 -d"
  fi

  dx cat "$db_id" | $decompressor | tar -C "./$db_prefix" -xvf -

  if [ $(find "./$db_prefix" -type f | cut -d / -f 2,3 | sort | uniq | wc -l) -eq "1" ]; then
    # If tarball has top-level dir, then move contents of that dir up one dir.
    mv ./${db_prefix}/$(find "./$db_prefix" -type f | cut -d / -f 3 | uniq)/* ./${db_prefix}
    find "./$db_prefix" -type f
  fi

  du -sh "./$db_prefix"
}

function process_bam() {
  set -e -x -o pipefail

  kraken_db_prefix="$1"
  krona_taxonomy_db_prefix="$2"

  # stage input BAM
  bam_id=$(dx-jobutil-parse-link --no-project "$3")
  bam_name=$(dx describe --name "$3")

  mkdir -p "scratch/${bam_id}/krona"
  dx download -o "scratch/${bam_id}/${bam_name}" "$bam_id"

  # folder structure for multi-lane outputs uses lane metadata recorded
  # in BAM property at the end of demux
  lane=$(dx describe --json "$bam_id" | jq -r .properties.lane)
  if [ "$lane" == "null" ]; then
      output_root_dir="out/outputs/"
  else
      output_root_dir="out/outputs/lane_$lane/"
  fi

  output_filename_prefix="${bam_name%.bam}"
  output_filename_prefix="${output_filename_prefix%.cleaned}"
  output_root_dir="${output_root_dir}${output_filename_prefix}"
  mkdir -p "$output_root_dir"

  # Use Kraken to classify taxonomic profile of sample.
  viral-ngs metagenomics.py kraken \
                  "/user-data/scratch/${bam_id}/${bam_name}" \
                  "/user-data/${kraken_db_prefix}" \
                  --outReads "/user-data/${output_root_dir}/${output_filename_prefix}.kraken-classified.txt.gz" \
                  --outReport "/user-data/${output_root_dir}/${output_filename_prefix}.kraken-report.txt" \
                  --threads 4

  # Use Krona to visualize taxonomic profiling output from Kraken.
  viral-ngs metagenomics.py krona \
                  "/user-data/${output_root_dir}/${output_filename_prefix}.kraken-classified.txt.gz" \
                  "/user-data/${krona_taxonomy_db_prefix}" \
                  "/user-data/scratch/${bam_id}/krona/${output_filename_prefix}.krona-report.html" \
                  --noRank

  # Standalone html file output
  cp "scratch/${bam_id}/krona/${output_filename_prefix}.krona-report.html" "${output_root_dir}/"

  # Tar all html and attached js files for easy download
  tar cf "${output_root_dir}/${output_filename_prefix}.krona-report.tar" -C "scratch/$bam_id/krona" .

  # cleanup
  rm -rf "scratch/${bam_id}"
}
