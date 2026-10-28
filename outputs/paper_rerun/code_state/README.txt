base commit: e90f3f24772080ff8e7f8daae41f588fccf9edaf (origin/main e90f3f2)
frozen: Sun Oct 25 10:59:41 -03 2026
 M extra/ambient_supervised_hopf.py
 M extra/campaign_unit.py
 M extra/period_map.py
 M extra/poly_cs.py
 M extra/s4t2_reconstruction.py
 M extra/test_all.py
 M src/geometry/frames.py
?? experiments/
?? extra/inv_atlas.py
?? paper/
Sun Oct 25 12:02:30 -03 2026: campaign_unit.py refinement guard restricted to accepted rungs; acs_s4_s7/s8 and cspoly_s4_s7/s8 rerun
Sun Oct 25 17:06:54 -03 2026: S3xS1 homotopy classes. Added extra/hopf_invariant.py (+ test), experiments/hopf_class.py, Table 4 class split
  in experiments/make_figures.py / paper_numbers.py, run_all.py jobs acs_s3xs1_s9/s10, hopf_s3xs1_sob_b10, hopf_class
  (synced into code/ for the record; the running scheduler keeps its original job list). Run by hand from code/:
  acs_s3xs1_s9, acs_s3xs1_s10 (16:43), hopf_s3xs1_sob_b10 (16:57), hopf_class (repo, 17:00).
  The repo's extra/ambient_supervised_hopf.py also records base['hopf_invariant'] (n = 2; diagnostic only) -- not
  in code/, whose copy produced the S3xS1 runs.
Sun Oct 25 18:51:12 -03 2026: final state after the main campaign. make_figures.py: Fig. 3 legend below the axes (synced to code/). Diff/stat regenerated (tracked changes + new files experiments/, extra/hopf_invariant.py, extra/test_hopf_invariant.py; excludes extra/inv_atlas.py and paper/). Tests: 63/63 pass. Paper update list: outputs/paper_rerun/PAPER_UPDATES.md
Mon Oct 26 13:25:22 -03 2026: campaign complete (scheduler finished 1578.9 min, failed: []; four hopf_s5xs1_sob_* runs done 13:20). make_figures.py: tab_compute prints <1 for sub-minute jobs (synced to code/). Tables/figures/numbers regenerated (203 macros, 0 failures). Total 76 jobs, 6817 job-minutes. PAPER_UPDATES.md final.
