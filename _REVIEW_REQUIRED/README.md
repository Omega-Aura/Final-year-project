# Review Required — Human Decision Needed

Items that could not be classified confidently from their content. **Nothing here is deleted or
quarantined**; they are parked pending a decision.

## `plan_dual-target-flavonoid-cadd-pipeline-agai_6b40bf21.json` and `..._2.json`

**What they are.** Two versions of the original pipeline *plan*, machine-generated on
**17 Aug 2026** at the very start of the project:

| File | Created | `task_summary` |
|---|---|---|
| `..._2.json` | 12:07:29 | "Dual-target flavonoid CADD pipeline against TTBK1 and MAO-B" |
| `..._6b40bf21.json` | 12:25:03 | same, **"(revised)"** |

Both carry the same schema: `version`, `created_at`, `task_summary`, `agents`, `phases`,
`desired_outputs`, `feasibility`. The `phases` array holds the intended Phase 0–8 breakdown with
per-step titles and descriptions — e.g. *"Phase 0 — Target validation & plant selection"* →
*"Review target disease biology"* → *"Compile literature evidence for TTBK1 in pathological tau
…"*.

Combined size: **36 KB.**

**Why they are not classified.** They are not scientific outputs — no data, no result, nothing
downstream depends on them, and no script references them. But they are also not junk: they record
**what the project originally set out to do**, which is the only artefact that does. Comparing the
planned phases against what actually happened is how you can see where the project changed course,
and the 18-minute gap between the two versions captures an early revision of scope.

**The judgement call.** Are they worth keeping as provenance, or are they scaffolding that the
[`LOGBOOK.md`](../LOGBOOK.md) has fully superseded? The logbook records what *was* done with far more
fidelity; these record what was *intended* before any of it ran.

**Options:**

| Option | Result |
|---|---|
| Keep as provenance | Move to `09_manuscript/` or a `docs/` folder — useful if you ever want to write about how the project's scope evolved |
| Quarantine | Move to `_ARCHIVE_TO_DELETE/` for deletion with the rest |
| Delete now | 36 KB; the logbook covers everything they would be consulted for |

They cost almost nothing to keep, so the only real argument for removing them is tidiness.

---

## Also worth a decision (not moved)

### `07_mmgbsa/` was an empty directory

The workflow numbering implies MM-GBSA has its own step directory, but it was empty: `MMPBSA.py`
must run beside its matching topologies, so all current results live in `06_md/system_*/`.

It now holds [`md_mmgbsa_summary.csv`](../07_mmgbsa/md_mmgbsa_summary.csv) (all ten systems
collected from the primary `.dat` files), a `README.md`, and the prior-phase MM-GBSA results. If you
would rather the step not exist at all, the summary and README can move into `06_md/` and this
directory can go — but the numbered step is referenced in the workflow index and the manuscript
Methods mapping.

### `06_md/params/cand_003.acpype/` — superseded parameters

Free-conformer GAFF2 parameters, superseded by `params/docked/cand_003.acpype/` (built from the
actual docked pose, which is what every MD run used). Kept deliberately as provenance for the
charge-mismatch incident recorded in the logbook. Remove only if you are confident that history no
longer needs to be demonstrable.
