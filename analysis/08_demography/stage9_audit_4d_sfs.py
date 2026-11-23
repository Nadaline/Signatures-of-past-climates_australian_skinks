#### Packages ####
import os
import csv
import itertools
import numpy as np
import moments

#### Directories ####
project = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sfs_file = os.path.join(project, "results/demography/sfs/stage9_4d/K1_K2_K3_K4_projected_folded_4d.fs")
out_dir = os.path.join(project, "results/demography/moments/stage9_4d/tables")
os.makedirs(out_dir, exist_ok=True)

#### Load genomic frequency spectrum ####
data = moments.Spectrum.from_file(sfs_file)
expected_shape = (81, 21, 101, 81)
if tuple(data.shape) != expected_shape:
    raise ValueError("Unexpected Stage 9 SFS dimensions: " + str(data.shape))
if not bool(data.folded):
    raise ValueError("The Stage 9 empirical SFS must be folded.")

#### Build a topology-free four-way split reference ####
# This spectrum is only a technical reference for the folded-likelihood audit.
# It represents an instantaneous four-way split with no migration and no
# elapsed integration time. It is not a fitted biological demographic model.
sample_sizes = list(data.sample_sizes)
ancestral = moments.Spectrum(moments.LinearSystem_1D.steady_state_1D(sum(sample_sizes)))
expected = moments.Manips.split_1D_to_2D(ancestral, sample_sizes[0], sum(sample_sizes[1:]))
expected = moments.Manips.split_2D_to_3D_2(expected, sample_sizes[1], sum(sample_sizes[2:]))
expected = moments.Manips.split_3D_to_4D_3(expected, sample_sizes[2], sample_sizes[3])
expected = expected.fold()
expected.pop_ids = data.pop_ids
expected.mask = data.mask.copy()

#### Check folded likelihood compatibility ####
dimensions_match = tuple(data.shape) == tuple(expected.shape)
pop_ids_match = list(data.pop_ids) == list(expected.pop_ids)
model_folded = bool(expected.folded)
passed = bool(data.folded) and model_folded and dimensions_match and pop_ids_match

with open(os.path.join(out_dir, "stage9_4d_folding_audit.csv"), "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "dataset", "model_reference", "data_folded", "model_before_fold_status",
        "model_used_for_likelihood_folded", "dimensions_match", "pop_ids_match", "PASS_FAIL", "notes"
    ])
    writer.writerow([
        "K1_K2_K3_K4 projected 4D SFS", "four_way_split_no_migration_technical_reference",
        bool(data.folded), "unfolded after sequential split manipulation", model_folded,
        dimensions_match, pop_ids_match, "PASS" if passed else "FAIL",
        "Reference created only to validate folded dimensions and axes before any optimization."
    ])

with open(os.path.join(out_dir, "stage9_4d_sfs_qc.csv"), "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["metric", "value"])
    writer.writerow(["dimensions", "x".join(str(value) for value in data.shape)])
    writer.writerow(["folded", bool(data.folded)])
    writer.writerow(["sample_sizes", ";".join(str(value) for value in data.sample_sizes)])
    writer.writerow(["pop_ids", ";".join(data.pop_ids)])
    writer.writerow(["sfs_mass", float(data.sum())])
    writer.writerow(["positive_cells", int(np.sum(np.asarray(data) > 0))])
    writer.writerow(["negative_cells", int(np.sum(np.asarray(data) < 0))])
    writer.writerow(["moments_version", moments.__version__])

#### Save diagnostic two-dimensional marginals ####
# These are visual summaries of the same four-dimensional SFS, not replacement
# pairwise inputs. Saving them separately avoids loading the dense 4D object in
# the figure script.
array = np.ma.filled(np.asarray(data), 0)
with open(os.path.join(out_dir, "stage9_4d_pairwise_marginals.csv"), "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["pair", "population_1", "population_2", "allele_count_1", "allele_count_2", "sfs_mass"])
    for first, second in itertools.combinations(range(4), 2):
        drop_axes = tuple(axis for axis in range(4) if axis not in (first, second))
        marginal = array.sum(axis=drop_axes)
        for allele_count_1, allele_count_2 in np.ndindex(marginal.shape):
            writer.writerow([
                data.pop_ids[first] + "_" + data.pop_ids[second], data.pop_ids[first], data.pop_ids[second],
                allele_count_1, allele_count_2, float(marginal[allele_count_1, allele_count_2])
            ])

if not passed:
    raise RuntimeError("The Stage 9 folded-spectrum audit failed.")

print("Stage 9 four-dimensional folded-spectrum audit: PASS")
