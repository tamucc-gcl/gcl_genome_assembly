#!/usr/bin/env python3
"""Own a bounded exclusive-node pool; retire only this pool's idle queued jobs."""
import argparse
import json
import os
from pathlib import Path
import shlex
import signal
import subprocess
import sys
import time
import uuid

from contamination_workers import Results, atomic_json


def command(argv):
    return subprocess.check_output(argv, text=True).strip()


def active_jobs(jobs):
    if not jobs:
        return {}
    # Query the user's jobs: squeue -j can fail when the last requested job just exits.
    text = command(['squeue', '-h', '-u', str(os.getuid()), '-o', '%i|%T|%.128j'])
    found = {}
    for line in text.splitlines():
        identifier, status, name = (field.strip() for field in line.strip().split('|', 2))
        if identifier in jobs and name == jobs[identifier]['name']:
            found[identifier] = status
    return found


def cancel(jobs, states, pending_only=False):
    for identifier, status in states.items():
        if identifier in jobs and (not pending_only or status == 'PENDING'):
            outcome = subprocess.run(['scancel', identifier], capture_output=True, text=True)
            if outcome.returncode and identifier in active_jobs(jobs):
                raise RuntimeError(f'Could not stop owned job {identifier}: {outcome.stderr}')


def run(args):
    plan = json.loads(args.plan.read_text())
    store = Results(plan)
    remaining = list(plan['tasks'])
    args.out.mkdir(parents=True, exist_ok=True)
    registry = args.out / 'pool_jobs.json'
    # A retried coordinator must not leave an earlier copy of its workers running.
    if registry.exists():
        old = json.loads(registry.read_text())
        states = active_jobs(old)
        cancel(old, states)
        if states:
            raise RuntimeError('Previous pool jobs are stopping; resume once they have exited')
    jobs = {}
    interrupted = False

    def terminate(signum, frame):
        nonlocal interrupted
        interrupted = True

    handlers = {s: signal.signal(s, terminate) for s in (signal.SIGTERM, signal.SIGINT)}
    deadline = time.monotonic() + args.timeout
    try:
        remaining = [task for task in remaining if not store.completed(task)]
        if remaining:
            token = uuid.uuid4().hex[:12]
            for worker in range(1, args.nodes + 1):
                name = f'contam_pool_{token}_{worker}'
                summary = (args.out / f'worker_{worker}.json').resolve()
                exitcode = (args.out / f'worker_{worker}.exitcode').resolve()
                script = (args.out / f'worker_{worker}.sh').resolve()
                invocation = [sys.executable, str(args.helper.resolve()), 'worker', '--plan', str(args.plan.resolve()),
                              '--worker', str(worker), '--slots', str(args.slots), '--scratch-root', args.scratch,
                              '--out', str(summary)]
                script.write_text('#!/bin/bash\nset +e\n' + shlex.join(invocation) + '\ncode=$?\n' +
                    f"printf '%s\\n' \"$code\" > {shlex.quote(str(exitcode))}.tmp\n" +
                    f'mv {shlex.quote(str(exitcode))}.tmp {shlex.quote(str(exitcode))}\nexit "$code"\n')
                argv = ['sbatch', '--parsable', '--nodes=1', f'--cpus-per-task={args.slots * int(plan["threads"])}',
                        '--exclusive', '--mem=0', f'--partition={args.partition}', f'--time={args.walltime}',
                        f'--job-name={name}', f'--output={args.out.resolve()}/worker_{worker}.out',
                        f'--error={args.out.resolve()}/worker_{worker}.err', '--export=ALL']
                if args.nodelist:
                    argv.append('--nodelist=' + args.nodelist)
                if interrupted:
                    raise InterruptedError('Pool interrupted during submission')
                identifier = command(argv + [str(script)]).split(';')[0]
                if not identifier.isdigit():
                    raise ValueError(f'Unexpected sbatch response: {identifier}')
                jobs[identifier] = dict(name=name, worker=worker, summary=str(summary), exitcode=str(exitcode))
                atomic_json(registry, jobs)
                print(f'Submitted worker {worker}: Slurm job {identifier}', flush=True)
        while jobs:
            if interrupted:
                raise InterruptedError('Pool interrupted; verified chunk results are retained')
            remaining = [task for task in remaining if not store.completed(task)]
            states = active_jobs(jobs)
            print(f'Remaining chunks: {len(remaining)}; active requests: {len(states)}', flush=True)
            if not remaining:
                cancel(jobs, states, pending_only=True)
            if not states:
                # Exiting jobs may have failed, but the durable checksum-verified results are authoritative.
                break
            if time.monotonic() >= deadline:
                raise TimeoutError('Exclusive search pool timed out; verified chunk results are retained')
            time.sleep(args.poll_seconds)
        remaining = [task for task in remaining if not store.completed(task)]
        summaries = []
        for job in jobs.values():
            path = Path(job['summary'])
            if path.is_file():
                summaries.append(json.loads(path.read_text()))
        atomic_json(args.out / 'pool_summary.json', dict(chunks=len(plan['tasks']), remaining=len(remaining),
                    workers=summaries, jobs=jobs))
        if remaining:
            raise RuntimeError(f'{len(remaining)} chunks remain after the pool exited; inspect worker logs and resume')
    except BaseException:
        cancel(jobs, active_jobs(jobs))
        raise
    finally:
        for signum, handler in handlers.items():
            signal.signal(signum, handler)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--helper', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--nodes', type=int, required=True)
    parser.add_argument('--slots', type=int, required=True)
    parser.add_argument('--partition', default='normal')
    parser.add_argument('--nodelist', default='')
    parser.add_argument('--walltime', required=True)
    parser.add_argument('--scratch', default='/scratch')
    parser.add_argument('--timeout', type=float, default=4 * 86400)
    parser.add_argument('--poll-seconds', type=float, default=30)
    args = parser.parse_args()
    if min(args.nodes, args.slots, args.timeout, args.poll_seconds) <= 0:
        parser.error('Pool sizes and time limits must be positive')
    run(args)


if __name__ == '__main__':
    main()
