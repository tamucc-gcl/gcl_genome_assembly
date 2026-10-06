# Chimera integration rerun: all three scenarios pass

Reviewed the returned eligible, unresolved and manual archives. All have exit status 0, PASS.json and three completed Nextflow processes with no failures or parameter-undefined warnings. Published decisions, coordinate lifts and verification JSON are regular archived files rather than broken work-directory links.

| Scenario | Decision/selection | Output |
|---|---|---|
| eligible | Synthetic complete measurements produce one automatic UNJOIN_UNSUPPORTED action | 2 -> 3 records, 14 bp |
| unresolved | Missing evidence produces UNRESOLVED and no actions | 2 -> 2 records, identical input/output FASTA checksum |
| manual | Explicit reviewed action is applied despite an unresolved automatic assessment | 2 -> 3 records, 14 bp |

Both cutting routes produce identical output FASTA checksums (`d06eb771903f71166f4c23446a0eff2ba57badf6566186264fd63844434e5f4c`), identical piece coordinate lifts and exact-parent-reconstruction verification. The left piece is [0,6), including both gap Ns; the right piece is [6,10). The unrelated record remains [0,4). Applied-action provenance distinguishes automatic selection from human review. All VERIFY_RESULT processes confirm neutral split names, unplaced class, forward orientation and matching name-map/FASTA identities.

This completes the synthetic Nextflow integration check of adjudication and application. It does not validate chromosome-transition discovery, real-data eligibility, Hi-C controls, native graph assessment, post-cut reassignment, finishing or the full CHIMERA workflow. No production CTlk assembly was changed.

Next work is to complete the real evidence measurement/controls producer and post-cut reassignment, then exercise the complete workflow on isolated current CTlk inputs and peer assemblies. H01 remains a reviewed biological test case; the synthetic eligibility fixture does not make H01 auto-eligible. No further rerun of these unchanged synthetic controls is needed.
