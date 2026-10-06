# Integration batch 1511293

The unresolved scenario passed all three processes (adjudication, application, output verification) with exit status zero. The eligible scenario completed adjudication, then failed application. The manual scenario failed application; its simultaneous adjudication task was aborted.

Both failures used `--min-piece-bp 1000000` instead of the fixture's intended 4 bases. Nextflow logs explicitly warn that `chimera_min_piece_bp` and `publish_dir_mode` were undefined in the included modules. Fixture defaults were declared after module inclusion. This is a test harness parameter-binding error, not a failure of the sequence-length guard.

The harness now initializes parameters before module inclusion, and the launcher explicitly supplies the fixture minimum and copy publishing mode. Archive creation dereferences published links so evidence remains readable away from Crest. Production minimum-piece settings remain unchanged.

The unresolved archive includes PASS.json, but decision and verification artifacts were published as links into the excluded work directory. Their contents could not be independently reviewed from this archive. The two successful-cut scenarios still require a fresh Nextflow run. The local command-line integration test passes all three scenarios after the harness repair; it does not replace that cluster rerun.

Resubmit the same array launcher after syncing the changes. Return all three fresh review archives. Do not use this failed batch as evidence that automatic or manual cutting passed Nextflow integration.
