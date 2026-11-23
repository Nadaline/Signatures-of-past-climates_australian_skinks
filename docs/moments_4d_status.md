# Joint four-population moments: technical checks only

I inspected the folded joint SFS for K1–K4. Its axes contain 81 × 21 × 101 × 81 bins, corresponding to haploid projections of 80, 20, 100 and 80. The audited spectrum has mass 6,894, with 127,358 positive cells and no negative cells. See the [SFS QC](../results/summary_tables/moments_4d/stage9_4d_sfs_qc.csv) and its [audit script](../analysis/08_demography/stage9_audit_4d_sfs.py).

M0 and M1 each completed a one-iteration technical run, labelled `smoke_only`. M0 reached a parameter bound, and M1 recorded eight invalid evaluations. The subsequent pilot stopped at an input assertion checking SFS mass and the number of observed cells, before returning a pilot fit. See the [technical run summary](../results/summary_tables/moments_4d/stage10_smoke_5f4cad1bcc354adb_runs.csv).

I report no joint 4D AIC ranking, dates or evidence for shared contraction. Figure 7 and its supplement use the archived pairwise results. Before attempting the joint pilot again, I would check the runner's exact SFS and observation mask against the audited counts without relaxing the validation gate.
