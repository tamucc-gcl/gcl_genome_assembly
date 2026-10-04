# Reference-scoring resource comparison

Evidence read directly from the user-supplied archives:

- `chimera-review-20261003-061009-1505905.tar.gz`: pipeline trace and Nextflow log.
- `chimera-review-20261003-080345-1505984.tar.gz`: pipeline trace.

Run 1505905 submitted eight HARMONIZE_SCORE tasks within 0.1 seconds on
2026-10-02 at 15:28:52. The Nextflow log records all eight running concurrently.
All completed successfully. Run 1505984 reused these same task hashes from cache.

| Reference | Slurm job | Reported peak RSS | Runtime |
|---|---:|---:|---:|
| CBau hap2 | 1505907 | 36.2 GB | 37m 29s |
| CMat hap2 | 1505908 | 35.5 GB | 35m 38s |
| CBau hap1 | 1505909 | 36.5 GB | 36m 36s |
| CLim hap2 | 1505910 | 34.7 GB | 36m 8s |
| CTlk hap2 | 1505911 | 37.1 GB | 37m 7s |
| CMat hap1 | 1505912 | 35.8 GB | 34m 31s |
| CLim hap1 | 1505913 | 32.6 GB | 32m 3s |
| CTlk hap1 | 1505914 | 32 GB | 27m 56s |

Values preserve the units displayed by Nextflow. These are process peaks, not
simultaneous totals or complete node-memory measurements.

Commit a56430e changed only a docstring in harmonize_names.py with respect to
the scoring implementation. HARMONIZE_SCORE stages that script as an input and
performs minimap2 alignments and Python scoring in the same task, so this edit
invalidated the expensive task cache despite no change in executable logic.

The user's current isolated task 1506190 reported 34665468K (about 33.1 GiB)
while still running. This is within the previous successful range. There is no
demonstrated increase in alignment memory requirements or configured concurrency.

The cluster reports select/cons_tres with CR_CORE, and failed job allocations
omit memory despite explicit memory requests. Node contention is a plausible
explanation for signal-9 failures, but neither the archival traces nor current
node snapshots establish the historical kill cause. Successful-run node placement
and historical node OOM logs are the remaining evidence for that diagnosis.

Do not block chimera assessment on resolving this scheduler question. Keep the
current isolated diagnostic run progressing. Reference-scoring alignment/scoring
separation is a later cache-efficiency improvement, not a biological prerequisite
or a demonstrated correction for these kills.

Successful-run placement confirmed by the user's Slurm accounting:

- Jobs 1505907–1505913 all ran on crest-c011 from 15:28:53, overlapping for
  more than 32 minutes. Their individual Slurm MaxRSS values sum to about
  235.2 GiB; this is not a simultaneous node-memory measurement.
- Job 1505914 ran on crest-c014 during the same period.
- The failures were on crest-c019 and crest-c020. There is no demonstrated
  increase in configured scoring concurrency. Seven concurrent scoring jobs on
  one node previously succeeded, so concurrency alone cannot explain the change.
- Historical free/available memory, other workloads, node capacity differences,
  and kill logs remain unmeasured. Do not assert a specific node-level cause.

Placement query used:

```bash
sacct -j 1505907,1505908,1505909,1505910,1505911,1505912,1505913,1505914 \
  --format=JobID,State,NodeList,Start,End,MaxRSS -P
```
