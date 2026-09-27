# Queue logs

Stdout from the sequential MD/analysis queues, kept as the run record: start and finish timestamps
per system, the final `production.log` line, and the MM-GBSA `DELTA TOTAL` as each system completed.

| File | Driver | What it records |
|---|---|---|
| `mao_pose_queue.txt` | [`../../scripts/run_mao_pose_queue.sh`](../../scripts/run_mao_pose_queue.sh) | the four MAO pose-scan MD runs, 26–27 Sep |
| `mao_pose_analysis.txt` | [`../../scripts/run_mao_pose_analysis.sh`](../../scripts/run_mao_pose_analysis.sh) | MM-GBSA and cpptraj for three of them |
| `mao_replicates.txt` | [`../../scripts/run_mao_replicate_queue.sh`](../../scripts/run_mao_replicate_queue.sh) | the two velocity replicates and their analysis, with the replicate spreads printed at the end |

`mao_pose_queue.txt` covers the relaunch only. The first attempt was killed by a machine restart
0.77 ns into `system_MAOA_p3`; that system was re-run from the start, which is why the queue's own
log shows `system_MAOA_p2` being skipped. See the 2026-09-27 entries in
[`../../LOGBOOK.md`](../../LOGBOOK.md).
