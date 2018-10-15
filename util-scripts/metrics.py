#!/usr/bin/env python

import argparse
import csv
import sys
import concurrent.futures
import datetime
import re
import multiprocessing

import dxpy

top_level_execution_attributes_to_include = ['id', 'executableName', 'folder', 'name', 'state', 'launchedBy', 'parentAnalysis']

parser = argparse.ArgumentParser(
            description="""This returns output metrics from analyses or jobs run on DNAnexus.
                        If a project ID is given, information is returned for all jobs and analyses within the project."""
            )
parser.add_argument('csvfile', type=argparse.FileType('w'), help='Path of the metrics file to write.')
parser.add_argument('ids', type=str, nargs='+', help='IDs for which information should be returned. These can be project-<ID>, analysis-<ID>, or job-<ID>.')
parser.add_argument('--state', dest='states', nargs='+', choices=["done", "failed", "running", "terminated", "runnable"], default=None, help="Execution states to include when returning information for all jobs or analyses in a project. Note: 'runnable' means the item is waiting to be executed.")
parser.add_argument('--executableName', dest='executable_names', nargs='+', default=None, help="DNAnexus executable names to include. If omitted, all are included.")
parser.add_argument('--created_after', dest='created_after', type=str, default=None, help="The time (or offset) to use as a lower bound for inclusion.")
parser.add_argument('--created_before', dest='created_before', type=str, default=None, help="The time (or offset) to use as an upper bound for inclusion.")
parser.add_argument('--noDescendants', dest='no_descendants', action='store_true', help="Include top-level executions only. This is helpful when specifying a project-ID and child jobs are not desired in the output.")

def available_cpu_count():
    """
    Return the number of available virtual or physical CPUs on this system.
    The number of available CPUs can be smaller than the total number of CPUs
    when the cpuset(7) mechanism is in use, as is the case on some cluster
    systems.

    Adapted from http://stackoverflow.com/a/1006301/715090
    """
    try:
        with open('/proc/self/status') as f:
            status = f.read()
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$', status)
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return min(res, multiprocessing.cpu_count())
    except IOError:
        pass

    return multiprocessing.cpu_count()

if __name__ == "__main__":
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    analysis_ids = filter(lambda s: s.startswith("analysis-"), args.ids)
    project_ids  = filter(lambda s: s.startswith("project-"),  args.ids)
    job_ids      = filter(lambda s: s.startswith("job-"),      args.ids)

    executions = []

    project_job_ids = []
    if project_ids:
        for project_id in project_ids:
            if args.states:
                for state in args.states:
                    project_job_ids.extend([e["id"] for e in dxpy.find_executions(project=project_id, state=state, created_before=args.created_before, created_after=args.created_after)])
            else:
                project_job_ids.extend([e["id"] for e in dxpy.find_executions(project=project_id, created_before=args.created_before, created_after=args.created_after)])

    execution_ids_to_describe = list(set(analysis_ids+project_job_ids+job_ids))

    print("Reading {} total executions...".format(len(execution_ids_to_describe)))

    with concurrent.futures.ThreadPoolExecutor(max_workers=available_cpu_count()) as executor:
        for execution in executor.map(dxpy.describe, execution_ids_to_describe, chunksize=50):
            executions.append(execution)

    all_metrics = []
    keys_seen = set()
    for execution in executions:
        metrics = {}

        if args.no_descendants:
            if "parentAnalysis" in execution and execution["parentAnalysis"] is not None:
                continue
        if args.executable_names:
            if "executableName" in execution and execution["executableName"] not in args.executable_names:
                continue

        metrics=dict([(x,execution[x]) for x in top_level_execution_attributes_to_include if x in execution])
        metrics["created"] = datetime.datetime.utcfromtimestamp(float(execution["created"])/1000).isoformat()
        keys_seen.update(metrics.keys())

        for execution_key in ["output"]:
            if execution_key in execution and execution[execution_key] is not None:
                for key, value in execution[execution_key].items():
                    if type(value) in (int, float):
                        field_name=key.split(".")[-1]
                        metrics[field_name] = value
                        keys_seen.add(field_name)

        all_metrics.append(metrics)

    with args.csvfile as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=sorted(keys_seen))
        writer.writeheader()
        writer.writerows(all_metrics)

    print("Metrics written for {} execution objects.".format(len(all_metrics)))
