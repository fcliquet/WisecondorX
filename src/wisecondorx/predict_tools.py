# WisecondorX

import logging
import os

import numpy as np
from scipy.stats import norm
from sklearn.decomposition import PCA

from wisecondorx.overall_tools import exec_R, get_z_score

"""
Returns gender based on Gaussian mixture
model trained during newref phase.
"""


def predict_gender(sample, trained_cutoff):
    Y_fraction = float(np.sum(sample["24"])) / float(
        np.sum([np.sum(sample[x]) for x in sample.keys()])
    )
    if Y_fraction > trained_cutoff:
        return "M"
    else:
        return "F"


"""
Normalize sample for read depth and apply mask.
"""


def coverage_normalize_and_mask(sample, ref_file, ap):
    by_chr = []

    chrs = range(1, len(ref_file["bins_per_chr{}".format(ap)]) + 1)

    for chr in chrs:
        this_chr = np.zeros(ref_file["bins_per_chr{}".format(ap)][chr - 1], dtype=float)
        min_len = min(
            ref_file["bins_per_chr{}".format(ap)][chr - 1], len(sample[str(chr)])
        )
        this_chr[:min_len] = sample[str(chr)][:min_len]
        by_chr.append(this_chr)
    all_data = np.concatenate(by_chr, axis=0)
    all_data = all_data / np.sum(all_data)
    masked_data = all_data[ref_file["mask{}".format(ap)]]

    return masked_data


"""
Project test sample to PCA space.
"""


def project_pc(sample_data, ref_file, ap):
    pca = PCA(n_components=ref_file["pca_components{}".format(ap)].shape[0])
    pca.components_ = ref_file["pca_components{}".format(ap)]
    pca.mean_ = ref_file["pca_mean{}".format(ap)]

    transform = pca.transform(np.array([sample_data]))

    reconstructed = np.dot(transform, pca.components_) + pca.mean_
    reconstructed = reconstructed[0]
    return sample_data / reconstructed


"""
Defines cutoff that will add bins to a blacklist
depending on the within reference distances.
"""


def get_optimal_cutoff(ref_file, repeats):
    distances = ref_file["distances"]
    cutoff = float("inf")
    for i in range(0, repeats):
        mask = distances < cutoff
        average = np.average(distances[mask])
        stddev = np.std(distances[mask])
        cutoff = average + 3 * stddev
    return cutoff


"""
Within sample normalization. Cycles through a number
of repeats where z-scores define whether a bin is seen
as 'normal' in a sample (in most cases this means
'non-aberrant') -- if it is, it can be used as a
reference for the other bins.
"""


def normalize_repeat(test_data, ref_file, optimal_cutoff, ct, cp, ap):
    results_z = None
    results_r = None
    ref_sizes = None
    test_copy = np.copy(test_data)
    for i in range(3):
        results_z, results_r, ref_sizes = _normalize_once(
            test_data, test_copy, ref_file, optimal_cutoff, ct, cp, ap
        )

        test_copy[ct:][np.abs(results_z) >= norm.ppf(0.99)] = -1
    m_lr = np.nanmedian(np.log2(results_r))
    m_z = np.nanmedian(results_z)

    return results_z, results_r, ref_sizes, m_lr, m_z


def _normalize_once(test_data, test_copy, ref_file, optimal_cutoff, ct, cp, ap):
    masked_bins_per_chr = ref_file["masked_bins_per_chr{}".format(ap)]
    masked_bins_per_chr_cum = ref_file["masked_bins_per_chr_cum{}".format(ap)]
    results_z = np.zeros(masked_bins_per_chr_cum[-1])[ct:]
    results_r = np.zeros(masked_bins_per_chr_cum[-1])[ct:]
    ref_sizes = np.zeros(masked_bins_per_chr_cum[-1])[ct:]
    indexes = ref_file["indexes{}".format(ap)]
    distances = ref_file["distances{}".format(ap)]

    i = ct
    i2 = 0
    for chr in list(range(len(masked_bins_per_chr)))[cp:]:
        start = masked_bins_per_chr_cum[chr] - masked_bins_per_chr[chr]
        end = masked_bins_per_chr_cum[chr]
        chr_data = np.concatenate(
            (
                test_copy[: masked_bins_per_chr_cum[chr] - masked_bins_per_chr[chr]],
                test_copy[masked_bins_per_chr_cum[chr] :],
            )
        )

        for index in indexes[start:end]:
            ref_data = chr_data[index[distances[i] < optimal_cutoff]]
            ref_data = ref_data[ref_data >= 0]
            ref_stdev = np.std(ref_data)
            results_z[i2] = (test_data[i] - np.mean(ref_data)) / ref_stdev
            results_r[i2] = test_data[i] / np.median(ref_data)
            ref_sizes[i2] = ref_data.shape[0]
            i += 1
            i2 += 1

    return results_z, results_r, ref_sizes


"""
The means of sets of within-sample reference
distances can serve as inverse weights for
CBS, Z-scoring and plotting.
"""


def get_weights(ref_file, ap):
    inverse_weights = [np.mean(np.sqrt(x)) for x in ref_file["distances{}".format(ap)]]
    weights = np.array([1 / x for x in inverse_weights])
    return weights


"""
Unmasks results array.
"""


def inflate_results(results, rem_input):
    temp = [0 for x in rem_input["mask"]]
    j = 0
    for i, val in enumerate(rem_input["mask"]):
        if val:
            temp[i] = results[j]
            j += 1
    return temp


"""
Log2-transforms results_r. If resulting elements are infinite,
all corresponding possible positions (at results_r, results_z
and results_w are set to 0 (blacklist)).
"""


def log_trans(results, log_r_median, fix_zero_bins=False):
    for chr in range(len(results["results_r"])):
        results["results_r"][chr] = np.log2(results["results_r"][chr])

    results["results_r"] = [x.tolist() for x in results["results_r"]]

    for c in range(len(results["results_r"])):
        for i, rR in enumerate(results["results_r"][c]):
            if not np.isfinite(rR):
                if fix_zero_bins and results["results_w"][c][i] != 0:
                    results["results_r"][c][i] = -20.0
                    results["results_z"][c][i] = -100.0
                else:
                    results["results_r"][c][i] = 0
                    results["results_z"][c][i] = 0
                    results["results_w"][c][i] = 0
            if results["results_r"][c][i] != 0:
                results["results_r"][c][i] = results["results_r"][c][i] - log_r_median


"""
Applies additional blacklist to results_r, results_z
and results_w if requested.
"""


def apply_blacklist(rem_input, results):
    blacklist = _import_bed(rem_input)

    for chr in blacklist.keys():
        for s_e in blacklist[chr]:
            for pos in range(s_e[0], s_e[1]):
                if len(results["results_r"]) < 24 and chr == 23:
                    continue
                if pos >= len(results["results_r"][chr]) or pos < 0:
                    continue
                results["results_r"][chr][pos] = 0
                results["results_z"][chr][pos] = 0
                results["results_w"][chr][pos] = 0


def _import_bed(rem_input):
    bed = {}
    for line in open(rem_input["args"].blacklist):
        chr_name, s, e = line.strip().split("\t")
        if chr_name[:3] == "chr":
            chr_name = chr_name[3:]
        if chr_name == "X":
            chr_name = "23"
        if chr_name == "Y":
            chr_name = "24"
        chr = int(chr_name) - 1
        if chr not in bed.keys():
            bed[chr] = []
        bed[chr].append(
            [int(int(s) / rem_input["binsize"]), int(int(e) / rem_input["binsize"]) + 1]
        )
    return bed


"""
Executed CBS on results_r using results_w as weights.
Calculates segmental zz-scores.
"""


def exec_cbs(rem_input, results):
    json_cbs_dir = os.path.abspath(rem_input["args"].outid + "_CBS_tmp")

    json_dict = {
        "R_script": str("{}/include/CBS.R".format(rem_input["wd"])),
        "ref_gender": str(rem_input["ref_gender"]),
        "alpha": str(rem_input["args"].alpha),
        "binsize": str(rem_input["binsize"]),
        "seed": str(rem_input["args"].seed),
        "results_r": results["results_r"],
        "results_w": results["results_w"],
        "infile": str("{}_01.json".format(json_cbs_dir)),
        "outfile": str("{}_02.json".format(json_cbs_dir)),
    }

    results_c = _get_processed_cbs(exec_R(json_dict))
    segment_z = get_z_score(results_c, results)
    results_c = [
        results_c[i][:3] + [segment_z[i]] + [results_c[i][3]]
        for i in range(len(results_c))
    ]
    return results_c


def _get_processed_cbs(cbs_data):
    results_c = []
    for i, segment in enumerate(cbs_data):
        chr = int(segment["chr"]) - 1
        s = int(segment["s"])
        e = int(segment["e"])
        r = segment["r"]
        results_c.append([chr, s, e, r])

    return results_c


def _get_aberration_cutoff(beta, ploidy):
    loss_cutoff = np.log2((ploidy - (beta / 2)) / ploidy)
    gain_cutoff = np.log2((ploidy + (beta / 2)) / ploidy)
    return loss_cutoff, gain_cutoff


def _get_aberrant_segments(rem_input, results):
    aberrant = []
    for segment in results["results_c"]:
        chr_name = str(segment[0] + 1)
        ploidy = 2
        if (chr_name == "23" or chr_name == "24") and rem_input["ref_gender"] == "M":
            ploidy = 1
        if rem_input["args"].beta is not None:
            loss_c, gain_c = _get_aberration_cutoff(rem_input["args"].beta, ploidy)
            if float(segment[4]) > gain_c or float(segment[4]) < loss_c:
                aberrant.append(segment)
        elif not isinstance(segment[3], str):
            if abs(float(segment[3])) > rem_input["args"].zscore:
                aberrant.append(segment)
    return aberrant


def _run_cbs_on_segment(rem_input, residuals, weights):
    json_cbs_dir = os.path.abspath(rem_input["args"].outid + "_reseg_tmp")
    # CBS.R expects results_r to have 23 (female) or 24 (male) chromosome lists.
    # We put our segment data in "chr 1" and pad the rest with single-zero lists.
    # Using [0] instead of [] avoids R's 1:0 == c(1,0) bug; the zeros become NA
    # in CBS.R (line 41) and all-NA chromosomes are skipped (lines 56-62).
    n_chrs = 23  # use "F" gender for simplicity
    padded_r = [list(residuals)] + [[0] for _ in range(n_chrs - 1)]
    padded_w = [[float(w) for w in weights]] + [[0] for _ in range(n_chrs - 1)]
    json_dict = {
        "R_script": str("{}/include/CBS.R".format(rem_input["wd"])),
        "ref_gender": "F",
        "alpha": str(rem_input["args"].alpha),
        "binsize": str(rem_input["binsize"]),
        "seed": str(rem_input["args"].seed),
        "results_r": padded_r,
        "results_w": padded_w,
        "infile": str("{}_01.json".format(json_cbs_dir)),
        "outfile": str("{}_02.json".format(json_cbs_dir)),
    }
    cbs_out = exec_R(json_dict)
    sub_segs = []
    if cbs_out:
        for seg in cbs_out:
            # Only take segments from chr 1 (our actual data)
            if int(seg["chr"]) == 1:
                sub_segs.append((int(seg["s"]), int(seg["e"]), seg["r"]))
    return sub_segs


def resegment_aberrations(rem_input, results, min_bins=3):
    aberrant_segments = _get_aberrant_segments(rem_input, results)
    nested = []

    for seg in aberrant_segments:
        chrom, start, end, seg_z, seg_ratio = seg
        seg_len = end - start + 1

        if seg_len < min_bins * 2:
            continue

        bin_ratios = results["results_r"][chrom][start : end + 1]
        bin_weights = results["results_w"][chrom][start : end + 1]

        valid_count = sum(1 for r in bin_ratios if r != 0)
        if valid_count < min_bins * 2:
            continue

        residuals = [r - seg_ratio if r != 0 else 0 for r in bin_ratios]

        sub_segments = _run_cbs_on_segment(rem_input, residuals, bin_weights)

        for sub_s, sub_e, sub_ratio in sub_segments:
            sub_len = sub_e - sub_s + 1
            if sub_len < min_bins:
                continue
            if sub_len > seg_len * 0.9:
                continue

            abs_ratio = seg_ratio + sub_ratio

            if rem_input["args"].beta is not None:
                if abs(sub_ratio) < rem_input["args"].beta / 4:
                    continue
            else:
                if abs(sub_ratio) < 0.1:
                    continue

            abs_start = start + sub_s
            abs_end = start + sub_e

            bin_z = results["results_z"][chrom][abs_start : abs_end + 1]
            valid_z = [z for z in bin_z if z != 0]
            if valid_z:
                sub_z = float(np.mean(valid_z))
            else:
                sub_z = "nan"

            nested.append([chrom, abs_start, abs_end, sub_z, abs_ratio, "nested"])

    return nested


def resegment_all_segments(rem_input, results, min_bins=3):
    """Re-run CBS on ALL segments to detect embedded events.

    Unlike resegment_aberrations() which only processes aberrant segments,
    this function processes every segment using a two-pass approach:

    1. CBS pass: Re-runs CBS on segment residuals (same as resegment_aberrations)
       to find moderate-size embedded events.
    2. Scan pass: For large neutral segments where CBS has limited power, uses a
       z-score-based scanning approach to detect small focal CNVs (e.g. the
       KANSL1 duplication hidden in an 8 Mb neutral segment).

    Sub-segments from neutral parents are filtered using the same aberration
    thresholds (z-score or beta) as primary segments.
    """
    all_segments = results["results_c"]
    nested = []
    reseg_alpha = getattr(rem_input["args"], "resegment_alpha", None)

    # Determine which segments are aberrant for adaptive filtering
    aberrant_set = set()
    for seg in _get_aberrant_segments(rem_input, results):
        aberrant_set.add((seg[0], seg[1], seg[2]))

    for seg in all_segments:
        chrom, start, end, seg_z, seg_ratio = seg[:5]
        seg_len = end - start + 1
        is_aberrant = (chrom, start, end) in aberrant_set

        if seg_len < min_bins * 2:
            continue

        bin_ratios = results["results_r"][chrom][start : end + 1]
        bin_weights = results["results_w"][chrom][start : end + 1]

        valid_count = sum(1 for r in bin_ratios if r != 0)
        if valid_count < min_bins * 2:
            continue

        residuals = [r - seg_ratio if r != 0 else 0 for r in bin_ratios]

        # Pass 1: CBS on residuals (good for moderate-size events)
        sub_segments = _run_cbs_on_segment(rem_input, residuals, bin_weights, alpha_override=reseg_alpha)

        for sub_s, sub_e, sub_ratio in sub_segments:
            sub_len = sub_e - sub_s + 1
            if sub_len < min_bins:
                continue
            if sub_len > seg_len * 0.9:
                continue

            abs_ratio = seg_ratio + sub_ratio

            if is_aberrant:
                if rem_input["args"].beta is not None:
                    if abs(sub_ratio) < rem_input["args"].beta / 4:
                        continue
                else:
                    if abs(sub_ratio) < 0.1:
                        continue
            else:
                if not _passes_aberration_filter(rem_input, chrom, abs_ratio):
                    continue

            abs_start = start + sub_s
            abs_end = start + sub_e

            bin_z = results["results_z"][chrom][abs_start : abs_end + 1]
            valid_z = [z for z in bin_z if z != 0]
            if valid_z:
                sub_z = float(np.mean(valid_z))
            else:
                sub_z = "nan"

            if not is_aberrant and rem_input["args"].beta is None:
                if isinstance(sub_z, str) or abs(sub_z) < rem_input["args"].zscore:
                    continue

            nested.append([chrom, abs_start, abs_end, sub_z, abs_ratio, "nested"])

        # Pass 2: Z-score scanning for small focal events in large neutral segments
        # CBS struggles to detect <50 aberrant bins in >1000-bin segments;
        # scanning for runs of extreme z-scores catches these.
        if not is_aberrant and seg_len >= min_bins * 10:
            scan_results = _scan_for_focal_events(
                results, chrom, start, end, seg_ratio, min_bins, rem_input
            )
            nested.extend(scan_results)

    return nested


def _passes_aberration_filter(rem_input, chrom, abs_ratio):
    """Check if a ratio passes the aberration threshold for a given chromosome."""
    ploidy = 2
    chr_name = str(chrom + 1)
    if (chr_name == "23" or chr_name == "24") and rem_input["ref_gender"] == "M":
        ploidy = 1
    if rem_input["args"].beta is not None:
        loss_c, gain_c = _get_aberration_cutoff(rem_input["args"].beta, ploidy)
        return abs_ratio > gain_c or abs_ratio < loss_c
    return True  # z-score mode: ratio filter not applicable


def _scan_for_focal_events(results, chrom, start, end, seg_ratio, min_bins, rem_input):
    """Scan a segment for focal events using z-score runs.

    Identifies contiguous runs of bins where z-scores consistently exceed
    a per-bin threshold, then evaluates the cluster as a potential CNV.
    Only the core cluster (seed bins bridging small NaN gaps) is used —
    no aggressive expansion that would dilute the signal.
    """
    nested = []
    bin_z = results["results_z"][chrom][start : end + 1]
    bin_r = results["results_r"][chrom][start : end + 1]
    bin_w = results["results_w"][chrom][start : end + 1]
    seg_len = end - start + 1

    # Scan minimum bins: at least 5 to reduce FP from small noise clusters
    scan_min_bins = max(min_bins, 5)

    # Per-bin z-score threshold for seed selection: use the segment threshold
    # (same stringency as primary calling) to limit seeds to strong signals
    bin_z_thresh = rem_input["args"].zscore

    # Find seed bins: valid bins with |z-score| > bin_z_thresh
    seeds = []
    for i in range(seg_len):
        if bin_r[i] != 0 and abs(bin_z[i]) > bin_z_thresh:
            seeds.append(i)

    if not seeds:
        return nested

    # Cluster adjacent seeds, allowing gaps of up to 3 bins (NaN/dropout)
    clusters = []
    current_cluster = [seeds[0]]
    for i in range(1, len(seeds)):
        if seeds[i] - seeds[i - 1] <= 3:
            current_cluster.append(seeds[i])
        else:
            clusters.append(current_cluster)
            current_cluster = [seeds[i]]
    clusters.append(current_cluster)

    for cluster in clusters:
        if len(cluster) < scan_min_bins:
            continue

        # Use the span of the cluster (first seed to last seed)
        left = cluster[0]
        right = cluster[-1]

        sub_len = right - left + 1
        if sub_len < scan_min_bins:
            continue
        if sub_len > seg_len * 0.9:
            continue

        # Compute sub-segment statistics from valid bins in the span
        sub_z_vals = [bin_z[i] for i in range(left, right + 1) if bin_r[i] != 0]
        sub_r_vals = [bin_r[i] for i in range(left, right + 1) if bin_r[i] != 0]
        sub_w_vals = [bin_w[i] for i in range(left, right + 1) if bin_r[i] != 0]

        if len(sub_z_vals) < scan_min_bins:
            continue

        # Require consistent sign (all seeds should have same direction)
        signs = [1 if z > 0 else -1 for z in sub_z_vals]
        dominant_sign = 1 if sum(signs) > 0 else -1
        concordant = sum(1 for s in signs if s == dominant_sign)
        if concordant < len(signs) * 0.7:  # at least 70% same direction
            continue

        mean_z = float(np.mean(sub_z_vals))
        if sub_w_vals and sum(sub_w_vals) > 0:
            mean_r = float(np.average(sub_r_vals, weights=sub_w_vals))
        else:
            mean_r = float(np.mean(sub_r_vals))

        abs_start = start + left
        abs_end = start + right

        # Apply same aberration thresholds as primary segments
        if rem_input["args"].beta is not None:
            if not _passes_aberration_filter(rem_input, chrom, mean_r):
                continue
        else:
            if abs(mean_z) < rem_input["args"].zscore:
                continue

        nested.append([chrom, abs_start, abs_end, mean_z, mean_r, "nested"])

    return nested
