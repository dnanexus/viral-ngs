#!/bin/bash

main() {
    set -e -x -o pipefail

    # pull the viral-ngs docker image
    dx-docker pull broadinstitute/viral-ngs$viral_ngs_version
    ls -lhR /tmp/dx-docker-cache/
    find /tmp/dx-docker-cache -type f > /tmp/resources-manifest.txt
    # generate a script /usr/local/bin/viral-ngs to invoke the viral-ngs docker image
    cat > /usr/local/bin/viral-ngs <<EOF
#!/bin/bash
set -x
dx-docker run -v \$(pwd):/user-data --entrypoint /bin/bash broadinstitute/viral-ngs$viral_ngs_version -c "source /opt/viral-ngs/easy-deploy-viral-ngs.sh load && \$@"
EOF
    chmod +x /usr/local/bin/viral-ngs
    cat /usr/local/bin/viral-ngs
    echo /usr/local/bin/viral-ngs >> /tmp/resources-manifest.txt

    # upload a tarball with the new files
    rinsed_version=$(echo "$viral_ngs_version" | tr -d ":@")
    resources=`tar -c -v -z -T /tmp/resources-manifest.txt | \
               dx upload --brief -o "viral-ngs-${rinsed_version}.resources.tar.gz" -`

    dx-jobutil-add-output resources "$resources" --class=file
}
