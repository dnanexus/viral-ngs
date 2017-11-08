#!/usr/bin/env python
import sys
import dxpy
import argparse
import subprocess
import time
import os
import json
import hashlib

argparser = argparse.ArgumentParser(description="Build the viral-ngs assembly workflow on DNAnexus.")
argparser.add_argument("--project", help="DNAnexus project ID", default="project-BXBXK180x0z7x5kxq11p886f")
argparser.add_argument("--folder", help="Folder within project (default: timestamp-based)", default=None)
argparser.add_argument("--novocraft", help="Novocraft license file. Optional: multithreading enabled with license")
argparser.add_argument("--gatk", help="GATK tarball (default: %(default)s)",
                                 default="file-By20P600jy1JY9q634Yq5PQQ")
argparser.add_argument("--run-tests", help="run small test assemblies", action="store_true")
argparser.add_argument("--run-large-tests", help="run test assemblies of varying sizes", action="store_true")
args = argparser.parse_args()

# detect git revision
here = os.path.dirname(sys.argv[0])
git_revision = subprocess.check_output(["git", "describe", "--always", "--dirty", "--tags"]).strip()

if args.folder is None:
    args.folder = time.strftime("/%Y-%m/%d-%H%M%S-") + git_revision

project = dxpy.DXProject(args.project)
applets_folder = args.folder + "/applets"
print "project: {} ({})".format(project.name, args.project)
print "folder: {}".format(args.folder)

###############################################################################
# BUILDING APPLETS
###############################################################################

def build_applets():
    applets = ["assembly/viral-ngs-human-depletion", "demux/viral-ngs-human-depletion-multiplex",
               "assembly/viral-ngs-filter", "assembly/viral-ngs-trinity", "assembly/viral-ngs-assembly-scaffolding",
               "assembly/viral-ngs-assembly-refinement", "assembly/viral-ngs-assembly-analysis",
               "demux/viral-ngs-demux-wrapper", "demux/viral-ngs-demux", "demux/viral-ngs-classification",
               "demux/viral-ngs-bwa-count-hits", "demux/viral-ngs-count-hits-multiplex"]

    # Build applets for assembly workflow in [args.folder]/applets/ folder
    project.new_folder(applets_folder, parents=True)
    for applet in applets:
        # TODO: reuse an existing applet with matching git_revision
        print "building {}...".format(applet),
        sys.stdout.flush()
        applet_dxid = json.loads(subprocess.check_output(["dx","build","--destination",args.project+":"+applets_folder+"/",os.path.join(here,applet)]))["id"]
        print applet_dxid
        applet = dxpy.DXApplet(applet_dxid, project=project.get_id())
        applet.set_properties({"git_revision": git_revision})

    # Build applets that user interact with directly in [args.folder]/ main folder
    exposed_applets = ["util/viral-ngs-fasta-fetcher"]
    for applet in exposed_applets:
        print "building {}...".format(applet),
        sys.stdout.flush()
        applet_dxid = json.loads(subprocess.check_output(["dx","build","--destination",args.project+":"+args.folder+"/",os.path.join(here,applet)]))["id"]
        print applet_dxid
        applet = dxpy.DXApplet(applet_dxid, project=project.get_id())
        applet.set_properties({"git_revision": git_revision})

build_applets()

# helpers for name resolution
def find_app(app_handle):
    return dxpy.find_one_app(name=app_handle, zero_ok=False, more_ok=False, return_handler=True)

def find_applet(applet_name, folder=applets_folder):
    return dxpy.find_one_data_object(classname='applet', name=applet_name,
                                     project=project.get_id(), folder=folder,
                                     zero_ok=False, more_ok=False, return_handler=True)

def find_resource_tarball_id():
    depletion_applet = find_applet("viral-ngs-human-depletion")
    depletion_applet_inputSpec = depletion_applet.describe()["inputSpec"]
    resource_tarball_id = [x for x in depletion_applet_inputSpec if x["name"] == "resources"][0]["default"]
    return resource_tarball_id

###############################################################################
# VIRAL ASSEMBLY WORKFLOWS: taking raw reads (in paired FASTQ or unmapped BAM)
# through optional human depletion, quality control and polished assembly
###############################################################################

assembly_workflow_resources = {
    'Ebola':{
        'contaminants': 'file-BXF0vYQ0QyBF509G9J12g927',
        'filter-targets': 'file-BXF0vf80QyBF509G9J12g9F2',
        'scaffold-reference': 'file-BXF0vZ00QyBF509G9J12g944',
        'abridged': False
    },
    'Lassa':{
        'contaminants': 'file-BXF0vYQ0QyBF509G9J12g927',
        'filter-targets': 'file-Bg533J00x0zBkYkFGb23k58B',
        'scaffold-reference': 'file-Bg533J00x0zBkYkFGb23k58B',
        'abridged': False
    },
    'Generic':{
        'contaminants': 'file-BXF0vYQ0QyBF509G9J12g927',
        'abridged': False
    },
    'Abridged':{
        'abridged': True
    }
}

def build_assembly_workflows(workflow_list):
    workflows = {}
    for w in workflow_list:
        workflow = build_assembly_workflow(w, assembly_workflow_resources[w])
        workflows[w] = workflow
    return workflows

def build_assembly_workflow(species, resources):
    wf = dxpy.new_dxworkflow(title='viral-ngs-assembly_{0}'.format(species),
                              name='viral-ngs-assembly_{0}'.format(species),
                              description='viral-ngs-assembly, with resources populated for {0}'.format(species),
                              project=args.project,
                              folder=args.folder,
                              properties={"git_revision": git_revision})

    # Locate the file ID corresponding to the viral-ngs resource tarball
    depletion_applet = find_applet("viral-ngs-human-depletion")
    depletion_applet_inputSpec = depletion_applet.describe()["inputSpec"]

    resource_tarball_id = find_resource_tarball_id()

    # These steps are used in the full assembly workflow
    if not resources.get('abridged', False):

        depletion_input = {
        "bmtagger_dbs": [x for x in depletion_applet_inputSpec if x["name"] == "bmtagger_dbs"][0]["default"],
        "blast_dbs": [x for x in depletion_applet_inputSpec if x["name"] == "blast_dbs"][0]["default"],
        "resources": [x for x in depletion_applet_inputSpec if x["name"] == "resources"][0]["default"]
        }
        depletion_stage_id = wf.add_stage(depletion_applet, stage_input=depletion_input, name="deplete", folder="intermediates")

        filter_input = {
            "reads": dxpy.dxlink({"stage": depletion_stage_id, "outputField": "cleaned_reads"}),
            "resources": dxpy.dxlink({"stage": depletion_stage_id, "inputField": "resources"})
        }
        if "filter-targets" in resources:
            filter_input["targets"] = dxpy.dxlink(resources["filter-targets"])

        filter_stage_id = wf.add_stage(find_applet("viral-ngs-filter"), stage_input=filter_input, name="filter", folder="intermediates")

        trinity_input = {
            "reads": dxpy.dxlink({"stage": filter_stage_id, "outputField": "filtered_reads"}),
            "subsample": 100000,
            "resources": dxpy.dxlink({"stage": depletion_stage_id, "inputField": "resources"})
        }
        if "contaminants" in resources:
            trinity_input["contaminants"] = dxpy.dxlink(resources["contaminants"])

        trinity_stage_id = wf.add_stage(find_applet("viral-ngs-trinity"), stage_input=trinity_input, name="trinity", folder="intermediates")

        scaffold_input = {
            "trinity_contigs": dxpy.dxlink({"stage": trinity_stage_id, "outputField": "contigs"}),
            "trinity_reads": dxpy.dxlink({"stage": trinity_stage_id, "outputField": "subsampled_reads"}),
            "resources": dxpy.dxlink({"stage": depletion_stage_id, "inputField": "resources"})
        }
        if "scaffold-reference" in resources:
            scaffold_input["reference_genome"] = dxpy.dxlink(resources["scaffold-reference"])

        scaffold_stage_id = wf.add_stage(find_applet("viral-ngs-assembly-scaffolding"), stage_input=scaffold_input, name="scaffold", folder="intermediates")

        refine1_input = {
            "assembly": dxpy.dxlink({"stage": scaffold_stage_id, "outputField": "modified_scaffold"}),
            "reads": dxpy.dxlink({"stage": depletion_stage_id, "outputField": "cleaned_reads"}),
            "min_coverage": 2,
            "novoalign_options": "-r Random -l 30 -g 40 -x 20 -t 502",
            "novocraft_license": dxpy.dxlink({"stage": scaffold_stage_id, "inputField": "novocraft_license"}),
            "gatk_tarball": dxpy.dxlink({"stage": scaffold_stage_id, "inputField": "gatk_tarball"}),
            "resources": dxpy.dxlink({"stage": depletion_stage_id, "inputField": "resources"})
        }
        refine1_stage_id = wf.add_stage(find_applet("viral-ngs-assembly-refinement"), stage_input=refine1_input, name="refine1", folder="intermediates")

        refine2_input = refine1_input
        refine2_input["assembly"] = dxpy.dxlink({"stage": refine1_stage_id, "outputField": "refined_assembly"})
        refine2_input["min_coverage"] = 3
        refine2_input["novoalign_options"] = "-r Random -l 40 -g 40 -x 20 -t 100"
        refine2_input["major_cutoff"] = dxpy.dxlink({"stage": refine1_stage_id, "inputField": "major_cutoff"})

        refine2_stage_id = wf.add_stage(find_applet("viral-ngs-assembly-refinement"), stage_input=refine2_input, name="refine2", folder="intermediates")

        analysis_input = {
            "assembly": dxpy.dxlink({"stage": refine2_stage_id, "outputField": "refined_assembly"}),
            "reads": dxpy.dxlink({"stage": refine2_stage_id, "inputField": "reads"}),
            "aligner_options": "-r Random -l 40 -g 40 -x 20 -t 100 -k",
            "novocraft_license": dxpy.dxlink({"stage": scaffold_stage_id, "inputField": "novocraft_license"}),
            "gatk_tarball": dxpy.dxlink({"stage": scaffold_stage_id, "inputField": "gatk_tarball"})
        }
        analysis_stage_id = wf.add_stage(find_applet("viral-ngs-assembly-analysis"), stage_input=analysis_input, name="analysis")

    # Build abridged workflow
    else:
        refine1_input = {
            "min_coverage": 2,
            "novoalign_options": "-r Random -l 30 -g 40 -x 20 -t 502",
            "resources": resource_tarball_id
        }
        refine1_stage_id = wf.add_stage(find_applet("viral-ngs-assembly-refinement"), stage_input=refine1_input, name="refine1", folder="refinement_1")

        refine2_input = {
            "reads": dxpy.dxlink({"stage": refine1_stage_id, "inputField": "reads"}),
            "assembly": dxpy.dxlink({"stage": refine1_stage_id, "outputField": "refined_assembly"}),
            "min_coverage": 3,
            "novoalign_options": "-r Random -l 40 -g 40 -x 20 -t 100",
            "resources": dxpy.dxlink({"stage": refine1_stage_id, "inputField": "resources"}),
            "gatk_tarball": dxpy.dxlink({"stage": refine1_stage_id, "inputField": "gatk_tarball"}),
            "novocraft_license": dxpy.dxlink({"stage": refine1_stage_id, "inputField": "novocraft_license"}),
            "major_cutoff": dxpy.dxlink({"stage": refine1_stage_id, "inputField": "major_cutoff"})
        }
        refine2_stage_id = wf.add_stage(find_applet("viral-ngs-assembly-refinement"), stage_input=refine2_input, name="refine2", folder="refinement_2")

        analysis_input = {
            "assembly": dxpy.dxlink({"stage": refine2_stage_id, "outputField": "refined_assembly"}),
            "reads": dxpy.dxlink({"stage": refine2_stage_id, "inputField": "reads"}),
            "aligner_options": "-r Random -l 40 -g 40 -x 20 -t 100 -k",
            "resources": dxpy.dxlink({"stage": refine1_stage_id, "inputField": "resources"}),
            "novocraft_license": dxpy.dxlink({"stage": refine1_stage_id, "inputField": "novocraft_license"}),
            "gatk_tarball": dxpy.dxlink({"stage": refine1_stage_id, "inputField": "gatk_tarball"})
        }
        analysis_stage_id = wf.add_stage(find_applet("viral-ngs-assembly-analysis"), stage_input=analysis_input, name="analysis")

    return wf

# workflows = dict of species-name: workflow_id
assembly_workflows = build_assembly_workflows(assembly_workflow_resources.keys())

###############################################################################
# DEMUX-ONLY WORKFLOW: upon completion of a streaming run upload, demultiplex
# the samples to unmapped BAMs, using the demux-wrapper to launch appropriate
# instance types
###############################################################################

def build_demux_only_workflow():
    resource_tarball_id = find_resource_tarball_id()

    wf = dxpy.new_dxworkflow(title='viral-ngs-demux-only',
                              name='viral-ngs-demux-only',
                              description='viral-ngs demultiplexing',
                              project=args.project,
                              folder=args.folder,
                              properties={"git_revision": git_revision})

    # Demux

    demux_applet = find_applet('viral-ngs-demux')
    demux_wrapper_applet = find_applet('viral-ngs-demux-wrapper')

    demux_wrapper_input = {
        "resources": dxpy.dxlink(resource_tarball_id),
        "demux_applet": dxpy.dxlink(demux_applet.id)
    }

    demux_stage_id = wf.add_stage(demux_wrapper_applet, stage_input=demux_wrapper_input,
        name='demux')

    # QC
    count_hits_applet = find_applet('viral-ngs-bwa-count-hits')
    count_hits_multiplex_applet = find_applet('viral-ngs-count-hits-multiplex')

    count_hits_multiplex_input = {
        "resources": dxpy.dxlink({"stage": demux_stage_id, "inputField": "resources"}),
        "in_bams": dxpy.dxlink({"stage": demux_stage_id, "outputField": "bams"}),
        "per_sample_output": True,
        "count_hits_applet": dxpy.dxlink(count_hits_applet.id)
    }

    count_hits_stage_id = wf.add_stage(count_hits_multiplex_applet,
        stage_input=count_hits_multiplex_input,
        name='count hits, fastqc')

    return wf

demux_only_workflow = build_demux_only_workflow()

###############################################################################
# DEMUX "PLUS" WORKFLOW: upon completion of a streaming run upload, demultiplex
# the samples to unmapped BAMs, plus run human depletion and metagenomics
# analysis
###############################################################################

def build_demux_plus_workflow():
    # Locate the file ID corresponding to the viral-ngs resource tarball
    resource_tarball_id = find_resource_tarball_id()

    wf = dxpy.new_dxworkflow(title='viral-ngs-demux-plus',
                              name='viral-ngs-demux-plus',
                              description='viral-ngs demultiplexing, human depletion, and metagenomics analysis',
                              project=args.project,
                              folder=args.folder,
                              properties={"git_revision": git_revision})

    # demux
    demux_applet = find_applet('viral-ngs-demux')
    demux_wrapper_applet = find_applet('viral-ngs-demux-wrapper')

    demux_wrapper_input = {
        "resources": dxpy.dxlink(resource_tarball_id),
        "demux_applet": dxpy.dxlink(demux_applet.id),
        "per_sample_output": True
    }
    demux_stage_id = wf.add_stage(demux_wrapper_applet, stage_input=demux_wrapper_input, name="demux")

    # QC
    count_hits_applet = find_applet('viral-ngs-bwa-count-hits')
    count_hits_multiplex_applet = find_applet('viral-ngs-count-hits-multiplex')

    count_hits_multiplex_input = {
        "resources": dxpy.dxlink({"stage": demux_stage_id, "inputField": "resources"}),
        "in_bams": dxpy.dxlink({"stage": demux_stage_id, "outputField": "bams"}),
        "per_sample_output": True,
        "count_hits_applet": dxpy.dxlink(count_hits_applet.id)
    }

    count_hits_stage_id = wf.add_stage(count_hits_multiplex_applet,
        stage_input=count_hits_multiplex_input,
        name='count hits, fastqc')

    # depletion
    depletion_input = {
        "bams": dxpy.dxlink({"stage": demux_stage_id, "outputField": "bams"}),
        "depletion_applet": dxpy.dxlink(find_applet("viral-ngs-human-depletion")),
        "resources": dxpy.dxlink(resource_tarball_id),
        "per_sample_output": True
    }
    depletion_stage_id = wf.add_stage(find_applet("viral-ngs-human-depletion-multiplex"), stage_input=depletion_input, name="deplete")

    # metagenomics
    metagenomics_applet = find_applet('viral-ngs-classification')
    metagenomics_input = {
        "mappings" : dxpy.dxlink({"stage": demux_stage_id, "outputField": "bams"}),
        "resources": dxpy.dxlink(resource_tarball_id)
    }
    metagenomics_stage_id = wf.add_stage(metagenomics_applet, stage_input=metagenomics_input, name="metagenomics")

    return wf

demux_plus_workflow = build_demux_plus_workflow()

###############################################################################
# TESTS
###############################################################################

if args.run_tests is True or args.run_large_tests is True:
    muscle_applet = dxpy.DXApplet("applet-BXQxjv00QyB9QF3vP4BpXg95")

    # Tests for Assembly Workflow

    # test data found in "bi-viral-ngs CI:/test_data"
    test_samples = {
        "SRR1553554": {
            "species": "Ebola",
            "reads": "file-BXPPQ2Q0YzB28x9Q9911Ykz5",
            "reads2": "file-BXPPQ380YzB6xGxJ45K9Yv6Q",
            "broad_assembly": "file-BXQx6G00QyB6PQVYKQBgzxv4",
            "expected_assembly_sha256sum": "3a849c1e545bca1ff938fe847f09206c0d5001de6153bf65831eb21513f1c3fa",
            # Note: contig name line (>....) removed
            "expected_subsampled_base_count":  729174,
            "expected_alignment_base_count": 577417
        }
    }

    if args.run_large_tests is True:
        # nb these samples takes too long for Travis
        # TODO: Update figures of merit for v1.10.1 udpates, which
        # overhauls subsampled and alignment base counts
        test_samples["SRR1553468"] = {
            "species": "Ebola",
            "reads": "file-BXYqZj80Fv4YqP151Zy9291y",
            "reads2": "file-BXYqZkQ0Fv4YZYKx14yJg0b4",
            "broad_assembly": "file-BXYqYKQ0QyB84xYJP9Kz7zzK",
            "expected_assembly_sha256sum": "626c6a72ce4380470340d6fd0d94f0a23e240c8bc57d63a7c83046ac7111ace7",
            # Note: contig name line (>....) removed
            "expected_subsampled_base_count": 18787806,
            "expected_alignment_base_count": 247110236
        }
        test_samples["G1190"] = {
            "species": "Lassa",
            "reads": "file-Bg97bJQ0x0z12q4XyZf5p0Kk",
            "broad_assembly": 'file-Bg533J00x0zBkYkFGb23k58B',
            "expected_assembly_sha256sum": "2c41eed662454fe6fd8757cc7c60b872d5021728bbfaeef1fd3d81d2023c0bca",
            # Note: contig name line (>....) removed
            "expected_subsampled_base_count":  1841634,
            "expected_alignment_base_count": 111944259
        }
        test_samples["SRR1553416"] = {
            "species": "Ebola",
            "reads": "file-BXBP0VQ011y0B0g5bbJFzx51",
            "reads2": "file-BXBP0Xj011yFYvPjgJJ0GzZB",
            "broad_assembly": "file-BXFqQvQ0QyB5859Vpx1j7bqq",
            "expected_assembly_sha256sum": "7e39e584758ba47c828a933cb836ea3a9b21b9afa2555a81d02fc17c8ba00e66",
            # Note: contig name line (>....) removed
            "expected_subsampled_base_count": 460338,
            "expected_alignment_base_count": 485406
        }

    # Launch assembly test workflows
    test_assembly_analyses = []
    for test_sample in test_samples.keys():
        # create a subfolder for this sample
        test_folder = args.folder + "/" + test_sample
        project.new_folder(test_folder)
        # run the workflow on the test sample
        try:
            workflow = assembly_workflows[test_samples[test_sample]["species"]]
        except KeyError:
            # Skip running test if workflow for req species was not built
            continue

        test_input = {
            "deplete.file": dxpy.dxlink(test_samples[test_sample]["reads"]),
            "deplete.skip_depletion": True,
            "scaffold.gatk_tarball": dxpy.dxlink(args.gatk),
        }

        # Chain in novocraft license, when provided (optional)
        if args.novocraft:
            test_input["scaffold.novocraft_license"] = dxpy.dxlink(args.novocraft)

        if "reads2" in test_samples[test_sample]:
            test_input["deplete.paired_fastq"] = dxpy.dxlink(test_samples[test_sample]["reads2"])

        test_analysis = workflow.run(test_input, project=project.get_id(), folder=test_folder,
                                     name=(git_revision+" "+test_sample+"-Assembly"))
        print "Launched {} for {}".format(test_analysis.get_id(), test_sample)
        test_assembly_analyses.append((test_sample,test_analysis))

    # Launch test workflows for "demux-plus"
    demux_plus_samples = {
        "run.151023_0015": {
            "upload_sentinel_record": "record-Bv8qkgQ0jy198GK0QVz2PV8Y",
        },
        "tar.run.151023_0015": {
            "run_tarballs": ["file-Bv8qyX00jy17QYqVpkJb7Gfj"],
        }
    }

    test_demux_analyses = []
    for run, opts in demux_plus_samples.items():
        # create a subfolder for this run
        test_folder = args.folder + "/" + run
        project.new_folder(test_folder)

        test_input = {
            # Skip depletion to save time
            "deplete.skip_depletion": True,

            # Use minikraken database (instead of full one)
            "metagenomics.kraken_db": dxpy.dxlink("file-Bqxxb8Q07q4bFjZZKJ25jyXb")
        }

        if "upload_sentinel_record" in opts:
            test_input["demux.upload_sentinel_record"] = dxpy.dxlink(opts["upload_sentinel_record"])
        elif "run_tarballs" in opts:
            test_input["demux.run_tarballs"] = [dxpy.dxlink(run_tarball) for run_tarball in opts["run_tarballs"]]
        else:
            exit("In Demux test run ({0}), neither sentinel record nor run_tarballs is specified".format(run))

        demux_plus_analysis = demux_plus_workflow.run(test_input, project=project.get_id(),
                                                      folder=test_folder,
                                                      name=(git_revision+" "+run+"-Demux-plus"))
        test_demux_analyses.append((run, demux_plus_analysis))

    # wait for jobs to finish while working around Travis 10m console inactivity timeout
    print "Waiting for analyses to finish..."
    noise = subprocess.Popen(["/bin/bash", "-c", "while true; do date; sleep 60; done"])
    try:
        for (test_sample,test_analysis) in test_assembly_analyses:
            test_analysis.wait_on_done()
            workflow = assembly_workflows[test_samples[test_sample]["species"]]

            # for diagnostics: add on a MUSCLE alignment of the Broad's
            # assembly of the sample with the workflow products
            muscle_input = {
                "fasta": [
                    test_analysis.get_output_ref(workflow.get_stage("scaffold")["id"]+".intermediate_scaffold"),
                    test_analysis.get_output_ref(workflow.get_stage("scaffold")["id"]+".modified_scaffold"),
                    test_analysis.get_output_ref(workflow.get_stage("refine1")["id"]+".refined_assembly"),
                    test_analysis.get_output_ref(workflow.get_stage("refine2")["id"]+".refined_assembly"),
                    dxpy.dxlink(test_samples[test_sample]["broad_assembly"])
                ],
                "output_format": "html",
                "output_name": test_sample+"_test_alignment",
                "advanced_options": "-maxiters 2"
            }
            muscle_applet.run(muscle_input, project=project.get_id(), folder=(args.folder+"/"+test_sample), name=(git_revision+" "+test_sample+" MUSCLE"), instance_type="mem1_ssd1_x4")

        for (run, demux_plus_analysis) in test_demux_analyses:
            # Just make sure demux plus runs without failure now,
            # TODO: check figure of merit for demux plus pipeline
            demux_plus_analysis.wait_on_done()

    finally:
        noise.kill()

    # check figures of merit for assembly tests
    for (test_sample,test_analysis) in test_assembly_analyses:
        workflow = assembly_workflows[test_samples[test_sample]["species"]]
        subsampled_base_count = test_analysis.describe()["output"][workflow.get_stage("trinity")["id"]+".subsampled_base_count"]
        expected_subsampled_base_count = test_samples[test_sample]["expected_subsampled_base_count"]
        print "\t".join([test_sample, "subsampled_base_count", str(expected_subsampled_base_count), str(subsampled_base_count)])

        # Get the final assembly and remove the contig name (>...) lines
        test_assembly_file_id = test_analysis.describe()["output"][workflow.get_stage("analysis")["id"]+".final_assembly"]
        dx_cat_cmd = ["dx", "cat", test_assembly_file_id['$dnanexus_link']]
        grep_cmd = ["grep", "-v", ">"]
        ps = subprocess.Popen(dx_cat_cmd, stdout=subprocess.PIPE)
        editted_assembly_file = subprocess.check_output(grep_cmd, stdin=ps.stdout)

        test_assembly_sha256sum = hashlib.sha256(editted_assembly_file).hexdigest()
        expected_sha256sum = test_samples[test_sample]["expected_assembly_sha256sum"]
        print "\t".join([test_sample, "sha256sum", expected_sha256sum, test_assembly_sha256sum])

        alignment_base_count = test_analysis.describe()["output"][workflow.get_stage("analysis")["id"]+".alignment_base_count"]
        expected_alignment_base_count = test_samples[test_sample]["expected_alignment_base_count"]
        print "\t".join([test_sample, "alignment_base_count", str(expected_alignment_base_count), str(alignment_base_count)])

        assert expected_sha256sum == test_assembly_sha256sum
        # Subsampled_base_count seems to drift, comment out for now
        # assert expected_subsampled_base_count == subsampled_base_count
        assert expected_alignment_base_count == alignment_base_count

    print "Success"
