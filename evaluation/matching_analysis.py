import os,sys
import argparse

import random

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd
from matplotlib import pyplot
from collections import Counter
from typing import Iterator

from modules import indexing_Maier, help_functions


def get_match_coverage(seq_len, mers, matches, order, span, partial_matches=dict()):
    """
        Span is k if kmer, length between first and last position if spaced kmer,
        and span between the start of the first strobe and the last nucleotide of the last strobe.
    """
    if not matches:
        return 0

    if order == 1:
        match_list = sorted([ (p, p+span) for (p, h) in mers.items() if h in matches ])
    elif order == 2:
        match_list = [ (p1,p2+span) for ((p1,p2), h) in mers.items() if h in matches ]
        match_list += [ (p1,p1 + span) for ((p1,_), h) in mers.items() if h in partial_matches[0]]
    elif order == 3:
        match_list = [ (p1, p3+span) for ((p1,p2, p3), h) in mers.items() if h in matches ]
        match_list += [(p1, p2 + span) for ((p1, p2, _), h) in mers.items() if h in partial_matches[1]]
        match_list += [(p1, p1 + span) for ((p1, _, _), h) in mers.items() if h in partial_matches[0]]

    match_list = sorted(match_list)

    covered_bases = match_list[0][1] - match_list[0][0]
    max_stop = match_list[0][1]
    for start,stop in match_list:
        if start < max_stop and stop > max_stop:
            covered_bases += stop - max_stop
            max_stop = stop
        elif start >= max_stop:
            covered_bases += stop - start
            max_stop = stop
    # print(covered_bases)
    return covered_bases #round(100*covered_bases/seq_len, 1)


def seq_covered_spaced_kmers(mers, matches, seq, positions):
    """
        Function specific to calculate the coverage of spaced k-mers
        since we have all the "sampled" positions in a spaced k-mer
        we can keep an active set of covered positions as we iterate 
        through the string.
    """
    seq_covered = 0
    if not matches:
        return seq_covered

    spaced_kmer_sampled_pos = sorted(positions)
    covered_pos = [0]*len(seq)
    all_match_pos_vector = sorted([p for (p, k) in mers.items() if k in matches ])

    for p in all_match_pos_vector:
        for j in spaced_kmer_sampled_pos:
            covered_pos[p + j] +=1

    # all positions covered at least once
    c = sum([1 for p in covered_pos if p > 0])
    return c


def get_sequence_coverage(mers, matches, order, k_len, partial_matches=dict()):
    covered_bases = 0
    if not matches:
        return covered_bases

    if order == 1:
        all_pos_vector = sorted([p for (p, k) in mers.items() if k in matches ])
    elif order == 2:
        match_list = [ [p1,p2] for ((p1,p2), k) in mers.items() if k in matches ]
        match_list += [ [p1] for ((p1,_), k) in mers.items() if k in partial_matches[0]]
        all_pos_vector = sorted([p for sublist in match_list for p in sublist])

    elif order == 3:
        match_list = [ [p1,p2, p3] for ((p1,p2, p3), k) in mers.items() if k in matches ]
        match_list += [[p1, p2] for ((p1, p2, _), k) in mers.items() if k in partial_matches[1]]
        match_list += [ [p1] for ((p1, _, _), k) in mers.items() if k in partial_matches[0]]
        all_pos_vector = sorted([p for sublist in match_list for p in sublist])


    prev_p = all_pos_vector[0]
    covered_bases += k_len # for first pos
    if len(all_pos_vector) == 1:
        return covered_bases

    for p in all_pos_vector[1:]:
        if p <= prev_p + k_len - 1:
            covered_bases += p-prev_p
        else:
            covered_bases += k_len

        prev_p = p

    return covered_bases



def get_intervals(mers, matches, order, partial_matches=dict()):
    if not matches:
        return [], []

    if order == 1:
        all_pos_vector = sorted([p for (p, k) in mers.items() if k in matches])
    elif order == 2:
        # match_list = [ [p1,p2] for ((p1,p2), k) in mers.items() if k in matches ]
        # all_pos_vector = sorted([p for sublist in match_list for p in sublist])
        match_list = [[p for p in range(p1, p2+1)] for ((p1, p2), k) in mers.items() if k in matches]
        match_list += [[p for ((p, _), k) in mers.items() if k in partial_matches[0]]]
        all_pos_vector = sorted(set([p for sublist in match_list for p in sublist]))

    elif order == 3:
        # match_list = [ [p1,p2, p3] for ((p1,p2, p3), k) in mers.items() if k in matches ]
        # all_pos_vector = sorted([p for sublist in match_list for p in sublist])
        match_list =  [[p for p in range(p1, p3 + 1)] for ((p1, p2, p3), k) in mers.items() if k in matches]
        match_list += [[p for p in range(p1, p2 + 1)] for ((p1, p2, _), k) in mers.items() if k in partial_matches[1]]
        match_list += [[p for ((p, _, _), k) in mers.items() if k in partial_matches[0]]]
        all_pos_vector = sorted(set([p for sublist in match_list for p in sublist]))

    # ivls = []
    # iv_start = all_pos_vector[0]
    # p_prev = all_pos_vector[0]
    # for p in all_pos_vector[1:]:
    #     if p == p_prev + 1:
    #         p_prev += 1
    #     else:
    #         ivls.append((iv_start, p_prev))
    #         p_prev = p
    #         iv_start = p
    # ivls.append((iv_start, p_prev))

    ivls = []
    iv_start = all_pos_vector[0]
    length = 0
    for p1,p2 in zip(all_pos_vector[:-1], all_pos_vector[1:]):
        if p2 == p1 + 1:
            length += 1
        # elif p2 == p1: # for the strobes
        #     pass
        elif p2 > p1 + 1:
            ivls.append((iv_start, iv_start+length))
            length = 0
            iv_start = p2

    if len(all_pos_vector) > 1:
        if p2 <= p1 + 1:
            ivls.append((iv_start, iv_start+length))
    elif len(all_pos_vector) == 1:
        ivls.append((iv_start, iv_start))
    # print(ivls)
    return ivls, all_pos_vector


def statistics(ivls, seq, k):
    if not ivls:
        return 1, [len(seq)], 0
    seq_covered = 0
    nr_islands = 0
    gap_lengths = []

    prev_stop = 0 #ivls[0][1] + k
    for i, (start, stop) in enumerate(ivls):
        if i == 0:
            seq_covered += (stop-start) + k
            if start > 0:
                nr_islands += 1
                gap_lengths.append(start - prev_stop)

        elif start > prev_stop + k:
            seq_covered += (stop-start) + k
            nr_islands += 1
            # Extra check to see if the gap is at least 1 nt. We may end up here because 
            # it may be that start = prev_stop + k + 1, which leads to a 0nt gap
            # (0nt can occur because of homopolymer stretches)
            # if start > prev_stop + k + 1: 
            gap_lengths.append(start - prev_stop - k)
            assert start - prev_stop >= 0, "found: {0},{1},{2}: {3}".format(start, prev_stop, k, ivls)
        else: # overlapping previous hit
            seq_covered += (stop-start) + k - (k - (start - prev_stop) )

        prev_stop = stop

    if ivls[-1][0] + k - 1 < len(seq):
        nr_islands +=1
        gap_lengths.append(len(seq) - prev_stop - k + 1)

    return nr_islands, gap_lengths, seq_covered


def multi_context_lookup(query_mer: tuple, order_to_ref_seeds: dict) -> tuple | None:
    """
    Multi context query seed lookup
    :param query_mer: tuple containing (<partial hash>, <full hash>) of the multi-context seed
    :param order_to_ref_seeds: dictionary from the search level to reference hash values
    :returns: in case of a full match: tuple (None, <full hash>, <number of hits in reference>)
              in case of a match: tuple (<search level>, <partial hash>, <number of hits in reference>)
              in case of no match: None
    """
    query_values = query_mer
    # print(query_mer)
    order = len(query_values)
    for search_level in range(order - 1, -1, -1):
        query_value = query_values[search_level]
        if query_value in order_to_ref_seeds[search_level]:
            return search_level, query_value, order_to_ref_seeds[search_level]
    return None

def get_multi_context_matches(ref_mers: Iterator[tuple], query_mers: Iterator[tuple], order: int):
    order_to_ref_seeds = dict()
    for i in range(order):
        order_to_ref_seeds[i] = Counter([mer[i] for mer in ref_mers])
    lookup_results = [multi_context_lookup(q, order_to_ref_seeds) for q in query_mers]
    full_hits = set()
    partial_hits = dict()
    for i in range(order - 1):
        partial_hits[i] = set()
    for lookup_result in lookup_results:
        if lookup_result:
            if lookup_result[0] == order - 1:
                full_hits.add(lookup_result[1])
            else:
                partial_hits[lookup_result[0]].add(lookup_result[1])
    # partial_hits = {lookup_result[0]: for lookup_result in lookup_results if lookup_result and lookup_result[0]}
    # full_hits = {lookup_result[1] for lookup_result in lookup_results if lookup_result and lookup_result[0] == order}

    partial_matches = dict()
    covered_ref_mers = set()

    for level, hits in partial_hits.items():
        partial_matches[level] = [ref_mer for ref_mer in ref_mers if ref_mer[level] in hits]
        covered_ref_mers |= set(partial_matches[level])
    full_matches = [ref_mer for ref_mer in ref_mers if ref_mer[order - 1] in full_hits]
    covered_ref_mers |= set(full_matches)
    return full_matches, partial_matches, len(covered_ref_mers)

def analyze_strobemers(seq1, seq2, k_size, order, hash_fcn, w, w_low = 0, w_high = 50):
    # minstrobes
    if order == 2:
        assert k_size % 2 == 0, "Not even kmer length, results will be different"
        if hash_fcn == "randstrobes":
            strobemers1 = indexing_Maier.randstrobes(seq1, k_size, w_low, w_high, w, order = 2 )
            strobemers2 = indexing_Maier.randstrobes(seq2, k_size, w_low, w_high, w, order = 2 )
            # print("randstrobes2",  len(strobemers1), len(strobemers2), len(set(strobemers1.values()) & set(strobemers2.values())))
            # print(strobemers1)
            # print(sorted(set(strobemers1.values()))[:20])
            # print(sorted(set(strobemers2.values()))[:20])
        elif hash_fcn == "minstrobes":
            strobemers1 = indexing_Maier.minstrobes(seq1, k_size, w_low, w_high, w, order = 2)
            strobemers2 = indexing_Maier.minstrobes(seq2, k_size, w_low, w_high, w, order = 2)
            # print("minstrobes2",  len(strobemers2))
        elif hash_fcn == "hybridstrobes":
            strobemers1 = indexing_Maier.hybridstrobes(seq1, k_size, w_low, w_high, w, order = 2)
            strobemers2 = indexing_Maier.hybridstrobes(seq2, k_size, w_low, w_high, w, order = 2)
            # print("minstrobes2",  len(strobemers2))
        elif hash_fcn == "multi-context":
            strobemers1 = indexing_Maier.multi_context(seq1, k_size, w_low, w_high, w, order = 2)
            strobemers2 = indexing_Maier.multi_context(seq2, k_size, w_low, w_high, w, order = 2)
    elif order == 3:
        assert k_size % 3 == 0, "Not div by 3 kmer length, results will be different"
        if hash_fcn == "randstrobes":
            strobemers1 = indexing_Maier.randstrobes(seq1, k_size, w_low, w_high, w, order = 3)
            strobemers2 = indexing_Maier.randstrobes(seq2, k_size, w_low, w_high, w, order = 3)
            # print("randstrobes3",  len(strobemers1), len(strobemers2), len(set(strobemers1.values()) & set(strobemers2.values())))
            # print(strobemers1)

        elif hash_fcn == "minstrobes":
            strobemers1 = indexing_Maier.minstrobes(seq1, k_size, w_low, w_high, w, order = 3)
            strobemers2 = indexing_Maier.minstrobes(seq2, k_size, w_low, w_high, w, order = 3)
            # print("minstrobes3",  len(strobemers2))

        elif hash_fcn == "hybridstrobes":
            strobemers1 = indexing_Maier.hybridstrobes(seq1, k_size, w_low, w_high, w, order = 3)
            strobemers2 = indexing_Maier.hybridstrobes(seq2, k_size, w_low, w_high, w, order = 3)
            # print("minstrobes2",  len(strobemers2))

        elif hash_fcn == "multi-context":
            strobemers1 = indexing_Maier.multi_context(seq1, k_size, w_low, w_high, w, order=3)
            strobemers2 = indexing_Maier.multi_context(seq2, k_size, w_low, w_high, w, order=3)

    # elif order == 4:
    #     assert k_size % 4 == 0, "Not div by 4 kmer length, results will be different"
    #     if hash_fcn == "randstrobes":
    #         strobemers1 = indexing_Maier.randstrobes(seq1, k_size, order = 4, w_1 = w_1, w_2 = w_2, w_3 = w_3)
    #         strobemers2 = indexing_Maier.randstrobes(seq2, k_size, order = 4, w_1 = w_1, w_2 = w_2, w_3 = w_3)
    #     # elif hash_fcn == "minstrobes":
    #     #     strobemers1 = indexing_Maier.minstrobes(seq1, k_size, order = 3, w_1 = w_1, w_2 = w_2)
    #     #     strobemers2 = indexing_Maier.minstrobes(seq2, k_size, order = 3, w_1 = w_1, w_2 = w_2 )
    # print(hash_fcn, order, len(strobemers2))
    full_matches = set(strobemers1.values()) & set(strobemers2.values())
    partial_matches = dict()
    for i in range(order - 1):
        partial_matches[i] = set()
    m = len(full_matches)
    mp = len(strobemers1.values())
    if hash_fcn == "multi-context":
        full_matches, partial_matches, m = get_multi_context_matches(strobemers1.values(), strobemers2.values(), order)
    ivls, all_pos_vector = get_intervals(strobemers1, full_matches, order, partial_matches)
    nr_islands, gap_lengths, c = statistics(ivls, seq1, k_size//order)
    seq_cov = get_sequence_coverage(strobemers1, full_matches, order, k_size//order, partial_matches)
    match_coverage = get_match_coverage(len(seq1), strobemers1, full_matches, order, k_size//order, partial_matches)

    # print("2-spaced minstrobes nr_matches:", len(matches2))  
    # print("Number of gaps (gaps):", nr_islands)
    # print("Avg island size (gaps):", sum(gap_lengths)/len(gap_lengths))
    # print("Seq covered:", seq_covered)
    # print("2-spaced minstrobes intervals:", ivls)
    return m, mp, seq_cov, gap_lengths, all_pos_vector, match_coverage


def analyze_kmers(seq1, seq2, k_size, w):
    #kmers
    kmers_pos1 = indexing_Maier.kmers(seq1, k_size, w)
    kmers_pos2 = indexing_Maier.kmers(seq2, k_size, w)
    # print("kmers", 1, len(kmers_pos1))
    # print("kmers:",  len(kmers_pos2))
    # kmers_pos1 = {p : seq1[i:i+k_size] for p, i in enumerate(range(len(seq1) - k_size +1))}
    # kmers_seq1 = set([seq1[i:i+k_size] for i in range(len(seq1) - k_size +1)])
    # kmers_seq2 = set([seq2[i:i+k_size] for i in range(len(seq2) - k_size +1)])
    # matches  = kmers_seq1 & kmers_seq2
    matches = set(kmers_pos1.values()) & set(kmers_pos2.values())
    m = len(matches)
    mp = len(kmers_pos1.values())
    ivls, all_pos_vector = get_intervals(kmers_pos1, matches, 1)
    nr_islands, gap_lengths, c = statistics(ivls, seq1, k_size)
    seq_cov = get_sequence_coverage(kmers_pos1, matches, 1, k_size)
    match_coverage = get_match_coverage(len(seq1), kmers_pos1, matches, 1, k_size)
    # print("kmers nr_matches:", len(matches))    
    # print("Number of gaps (gaps):", nr_islands)
    # print("Avg island size (gaps):", sum(gap_lengths)/len(gap_lengths))
    # print("Seq covered:", seq_covered)
    # print("kmer intervals:", ivls)
    return m, mp, seq_cov, gap_lengths, all_pos_vector, match_coverage


def analyze_spaced_kmers(seq1, seq2, k_size, span_size, w):
    positions = set(random.sample(range(1, span_size - 1 ), k_size-2))
    positions.add(0)
    positions.add(span_size - 1) # asserts first and last position is sampled so that we have a spaced kmer of length span size
    spaced_kmers_seq1 = indexing_Maier.spaced_kmers(seq1, k_size, span_size, positions, w)
    spaced_kmers_seq2 = indexing_Maier.spaced_kmers(seq2, k_size, span_size, positions, w)
    matches  = set(spaced_kmers_seq1.values()) & set(spaced_kmers_seq2.values())
    m = len(matches)
    mp = len(spaced_kmers_seq1.values())
    ivls, all_pos_vector = get_intervals(spaced_kmers_seq1, matches, 1)
    nr_islands, gap_lengths, _ = statistics(ivls, seq1, span_size)
    # we compute coverage for spaced k-mers with specific function
    seq_cov = seq_covered_spaced_kmers(spaced_kmers_seq1, matches, seq1, positions)
    match_coverage = get_match_coverage(len(seq1), spaced_kmers_seq1, matches, 1, span_size)

    # print("kmers nr_matches:", len(matches))    
    # print("Number of gaps (gaps):", nr_islands)
    # print("Avg island size (gaps):", sum(gap_lengths)/len(gap_lengths))
    # print("Seq covered:", seq_covered)
    # print("kmer intervals:", ivls)
    return m, mp, seq_cov, gap_lengths, all_pos_vector, match_coverage


# def plot_island_distribution(results, mut_freq, outfolder):
#     # pd.set_option("display.precision", 8)
#     filename = os.path.join(outfolder, "{0}.pdf".format(mut_freq))
#     plt.yscale('log', nonposy='clip')
#     # bins = [0.1*i for i in range(300)]
#     for label in results:
#         if label == "kmers" or label == "spaced kmers & dense" or label == "spaced kmers & sparse":
#             if label == "kmers":
#                 flat = [g for l in results[label]["gaps"] for g in l]
#                 pyplot.hist(flat, 100, range=[0, 500], alpha=0.5, label=label)
#         else:
#             for t in results[label]:
#                 flat = [g for l in results[label][t]["gaps"] for g in l]
#                 tmp_label = label + "-{0}".format(t)
#                 # if label == "randstrobes" and t == (3,10,25):
#                 #     pyplot.hist(flat, 100, range=[0, 500], alpha=0.5, label=tmp_label)
#                 pyplot.hist(flat, 100, range=[0, 500], alpha=0.5, label=tmp_label)

#     pyplot.legend(loc='upper right')
#     # pyplot.xlabel("Difference to genome (%)")
#     pyplot.xlabel("Island length")
#     pyplot.ylabel("Count")
#     plt.savefig(filename)
#     plt.savefig(filename)
#     plt.close()

def plot_island_distribution2(results, mut_freq, k_size, outfolder):
    # pd.set_option("display.precision", 8)
    filename = os.path.join(outfolder, "{0}.pdf".format(mut_freq))

    # make data correct format
    data = open(os.path.join(outfolder, "{0}.csv".format(mut_freq)), "w")
    data.write("label\tdp\tmut_freq\n")

    for label in results:
        non_strobemer_protocols = {"kmers & {}".format(k_size), "kmers & {}".format(k_size // 2),
                                   "kmers & {}".format(k_size * 2 // 3), "spaced kmers & sparse",
                                   "spaced kmers & dense"}
        if label in non_strobemer_protocols:
            flat = [g for l in results[label]["gaps"] for g in l]
            for dp in flat:
                data.write("{0}\t{1}\t{2}\n".format(label, dp, mut_freq))
        else:
            for t in results[label]:
                flat = [g for l in results[label][t]["gaps"] for g in l]
                tmp_label = label + "-{0}".format(t)
                for dp in flat:
                    data.write("{0}\t{1}\t{2}\n".format(tmp_label, dp, mut_freq))
    data.close()
    # hue_order = ["randstrobes-(3, 10, 25, 50)", "randstrobes-(2, 15, 25, 50)", "kmers", "minstrobes-(3, 10, 25, 50)", "minstrobes-(2, 15, 25, 50)", "spaced kmers & dense", "spaced kmers & sparse"]
    hue_order = ["randstrobes-(3, 10, 25, 50)", "hybridstrobes-(3, 10, 25, 50)", "kmers", "minstrobes-(2, 15, 25, 50)", "spaced kmers & dense"]
    data = pd.read_csv(data.name, sep='\t')
    # plt.yscale('log', nonposy='clip')
    # bins = [0.1*i for i in range(300)]
    sns.displot(data, x="dp", hue="label", hue_order = hue_order,
                 element="step", log_scale=(True, True)) # , kind="kde", log_scale= True, fill=True, multiple="stack"
    pyplot.legend(loc='upper right')
    # pyplot.xlabel("Difference to genome (%)")
    pyplot.xlabel("Gap length")
    pyplot.ylabel("Count")
    plt.savefig(filename)
    plt.savefig(filename)
    plt.close()


def print_matches(all_pos_vector, method):
    s = set(all_pos_vector)
    for i in range(100):
        if i in s:
            print("X",end='')
        else:
            print(" ",end='')
    print(method)
    print()

import numpy as np
def plot_matches(all_data, method, L, k_size,outfolder):

    data = np.random.randn(5, 2)
    print(data)
    binary_matrices = []
    for all_runs in all_data:
        binary_matrix = []
        for run in all_runs:
            binary_vector = []
            s = set(run)
            for i in range(83):
                if i in s:
                    # print("X",end='')
                    binary_vector.append(1)
                else:
                    # print(" ",end='')
                    binary_vector.append(0)
            binary_matrix.append(binary_vector)
        binary_matrices.append( np.array(binary_matrix)  )
    # print(binary_matrix)

    # np_matrix = np.array(binary_matrix)  
    fig, ax = plt.subplots(3,2,sharex=True, sharey=True)
    fig.suptitle('Match distribution')
    plt.yticks([])
    id_labels = ["1", "2", "3", "4", "5"]
    ax[0, 0].set_title('minstrobes (2,9,10,20)')
    mat = ax[0, 0].imshow(binary_matrices[0], cmap='GnBu', interpolation='nearest')
    # ax[0].set_yticks(range(binary_matrices[0].shape[0]), id_labels)

    ax[1, 0].set_title('minstrobes (3,6,10,20)')
    mat = ax[1, 0].imshow(binary_matrices[1], cmap='GnBu', interpolation='nearest')
    # ax[1].set_yticks(range(binary_matrices[1].shape[0]), id_labels)

    # ax[2].set_yticks(range(binary_matrices[2].shape[0]), id_labels)
    ax[2, 0].set_title('randstrobes (2,9,10,20)')
    mat = ax[2,0].imshow(binary_matrices[2], cmap='GnBu', interpolation='nearest')

    # ax[3].set_yticks(range(binary_matrices[3].shape[0]), id_labels)
    ax[0,1].set_title('randstrobes (3,6,10,20)')
    mat = ax[0,1].imshow(binary_matrices[3], cmap='GnBu', interpolation='nearest')

    # ax[2].set_yticks(range(binary_matrices[2].shape[0]), id_labels)
    ax[1, 1].set_title('hybridstrobes (2,9,10,20)')
    mat = ax[1,1].imshow(binary_matrices[4], cmap='GnBu', interpolation='nearest')

    # ax[3].set_yticks(range(binary_matrices[3].shape[0]), id_labels)
    ax[2,1].set_title('hybridstrobes (3,6,10,20)')
    mat = ax[2,1].imshow(binary_matrices[5], cmap='GnBu', interpolation='nearest')


    # plt.xticks(range(id_matrix.shape[1]), concert_dates)
    # plt.xticks(rotation=30)
    plt.xlabel('Position')

    # # this places 0 or 1 centered in the individual squares
    # for x in range(np_matrix.shape[0]):
    #     for y in range(np_matrix.shape[1]):
    #         ax.annotate(str(np_matrix[x, y])[0], xy=(y, x), 
    #                     horizontalalignment='center', verticalalignment='center')

    # ax = sns.heatmap(data, cbar=False, xticklabels = False, yticklabels=False)
    # ax.tick_params(left=False, bottom=False)
    filename = os.path.join(outfolder, "{0}_ex.pdf".format(method))
    plt.tight_layout()
    plt.savefig(filename)


def get_e_size(all_islands, L, nr_exp):
    # print("all_islands",all_islands)
    sum_of_squares = sum([x**2 for x in all_islands])
    return sum_of_squares/(L*nr_exp)

def print_combined_table(combined_results, mut_freqs, metrics_list, protocols):
    num_metrics = None
    num_protocols = None
    num_freqs = len(combined_results)
    protocol_to_results = {}
    protocols = []
    for mut_freq, results in combined_results.items():
        if not num_protocols:
            num_protocols = len(results)
        else:
            assert(num_protocols == len(results))
        protocols = results.keys()
        for protocol, protocol_results in results.items():
            if not num_metrics:
                num_metrics = len(protocol_results)
            else:
                assert(num_metrics == len(protocol_results))
            if protocol not in protocol_to_results:
                protocol_to_results[protocol] = {}
            protocol_to_results[protocol][mut_freq] = protocol_results
    # Generate header
    dataset_name = "SIM-R"
    table_string = ""
    table_string += "\\begin{table}[]\n"
    num_columns = num_freqs * num_metrics + 2
    column_string = "{" + "".join(["l" for _ in range(num_columns)]) + "}"
    table_string += "\\begin{tabular}" + column_string + "\n"
    table_string += " &  & " + "\\multicolumn{{{}}}{{c}}{{{}}}".format(num_freqs * num_metrics, dataset_name) + "\\\\ \\cline{{{}-{}}}\n".format(3, num_columns)
    table_string += " &  "
    for mut_freq in mut_freqs:
        table_string += "& \\multicolumn{{{}}}{{c}}{{{}}} ".format(num_metrics, mut_freq)
    table_string += "\\\\ \\cline{{{}-{}}}\n".format(3, num_columns)
    table_string += " & "
    for _ in mut_freqs:
        table_string += " & "
        table_string += " & ".join(metrics_list)
    table_string += " \\\\ \\hline\n"
    # Generate result
    print(protocols)
    for protocol, protocol_results in protocol_to_results.items():
        table_string += protocol + " & "
        print(protocol_results[0.01])
        mut_freq_results = [" & ".join([str(round(r, 1)) for r in protocol_results[mut_freq]]) for mut_freq in mut_freqs]
        table_string += " & ".join(mut_freq_results)
        table_string += " \\\\\n"

    print(table_string)


def main(args):
    L = 10000
    # L = 10000
    k_size = 30
    nr_exp = 100
    w = 1 # thinning, w = 1  means no thinning. w =1, 10, 20 was used in the evaluations.
    mut_freqs = [0.01, 0.05, 0.1] #[0.1]
    # mut_freqs = [0.1]
    w_2low = 25
    w_3low = 25
    w_2high = 50
    w_3high = 50

    # w_3strobe = 25
    # w_4strobe = 25
    # experiment_type choose between 'only_subs', 'controlled' or 'all'
    experiment_type = "all" #"controlled" # "all" #"only_subs" # "" # for spaced kmers
    # mut_freq = 0.5 #0.01 #, 0.05, 0.1]
    list_for_illustration = [[],[],[],[],[],[],[],[]]

    combined_results = {}
    protocols = []

    for mut_freq in mut_freqs:
        print("MUTATION RATE:", mut_freq)
        results = {"kmers & {}".format(k_size) : {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc":0},
                    "spaced kmers & dense" : {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc":0},
                    "spaced kmers & sparse" : {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc":0},
                    "minstrobes" : { (2,15,w_2low,w_2high): {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc":0 },
                                     (3,10,w_3low,w_3high): {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc":0} },
                    "randstrobes" : { (2,15,w_2low,w_2high): {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc":0},
                                      (3,10,w_3low,w_3high): {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc":0} },
                    "kmers & {}".format(k_size // 2): {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0},
                    "kmers & {}".format(k_size * 2 // 3): {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc": 0},
                    "multi-context" : { (2,15,w_2low,w_2high): {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc":0},
                                        (3,10,w_3low,w_3high): {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc":0} },
                    "hybridstrobes" : { (2,15,w_2low,w_2high): {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc":0 },
                                        (3,10,w_3low,w_3high): {"m": 0, "mp": 0, "sc": 0, "gaps": [], "mc":0 }}
                   }
        for exp_id in range(nr_exp):
            if exp_id % (nr_exp / 50) == 0:
                print("Performing iteration {} out of {}".format(exp_id, nr_exp))

            seq1 = "".join([random.choice("ACGT") for i in range(L)])

            # controlled or random experiment
            if experiment_type == 'only_subs':
                muts = set(random.sample(range(len(seq1)),int(L*mut_freq)))
                seq2 = "".join([seq1[i] if i not in muts else random.choice([help_functions.reverse_complement(seq1[i])]) for i in range(len(seq1))])
            elif experiment_type == 'controlled':
                # muts = set(range(15,L,15)) # every 15th nt for figure 2 only!
                muts = set(range(20,L,20))
                seq2 = "".join([seq1[i] if i not in muts else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
            elif experiment_type == 'all':
                muts = set(random.sample(range(len(seq1)),int(L*mut_freq)))
                seq2 = "".join([seq1[i] if i not in muts else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
            else:
                print("Wrong experiment label specified")
                sys.exit()



            # kmers
            m,mp,sc,gaps,all_pos_vector, match_coverage = analyze_kmers(seq1, seq2, k_size, w)
            results["kmers & {}".format(k_size)]["m"] += m
            results["kmers & {}".format(k_size)]["sc"] += sc
            results["kmers & {}".format(k_size)]["gaps"].append(gaps)
            results["kmers & {}".format(k_size)]["mc"] += match_coverage
            results["kmers & {}".format(k_size)]["mp"] += mp
            # print("kmers", match_coverage)
            # print_matches(all_pos_vector, "kmers")


            m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_kmers(seq1, seq2, k_size // 2, w)
            results["kmers & {}".format(k_size // 2)]["m"] += m
            results["kmers & {}".format(k_size // 2)]["sc"] += sc
            results["kmers & {}".format(k_size // 2)]["gaps"].append(gaps)
            results["kmers & {}".format(k_size // 2)]["mc"] += match_coverage
            results["kmers & {}".format(k_size // 2)]["mp"] += mp
            # print("kmers", match_coverage)
            # print_matches(all_pos_vector, "kmers")

            m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_kmers(seq1, seq2, k_size * 2 // 3, w)
            results["kmers & {}".format(k_size * 2 // 3)]["m"] += m
            results["kmers & {}".format(k_size * 2 // 3)]["sc"] += sc
            results["kmers & {}".format(k_size * 2 // 3)]["gaps"].append(gaps)
            results["kmers & {}".format(k_size * 2 // 3)]["mc"] += match_coverage
            results["kmers & {}".format(k_size * 2 // 3)]["mp"] += mp
            # print("kmers", match_coverage)
            # print_matches(all_pos_vector, "kmers")

            # Spaced kmers dense
            m,mp,sc,gaps,all_pos_vector, match_coverage = analyze_spaced_kmers(seq1, seq2, k_size, k_size+k_size//2, w)
            results["spaced kmers & dense"]["m"] += m
            results["spaced kmers & dense"]["sc"] += sc
            results["spaced kmers & dense"]["gaps"].append(gaps)
            results["spaced kmers & dense"]["mc"] += match_coverage
            results["spaced kmers & dense"]["mp"] += mp
            # print("spaced kmers & dense", match_coverage)

            # print_matches(all_pos_vector, "Spaced kmers")

            # Spaced kmers sparse
            m,mp,sc,gaps,all_pos_vector, match_coverage = analyze_spaced_kmers(seq1, seq2, k_size, 3*k_size, w)
            results["spaced kmers & sparse"]["m"] += m
            results["spaced kmers & sparse"]["sc"] += sc
            results["spaced kmers & sparse"]["gaps"].append(gaps)
            results["spaced kmers & sparse"]["mc"] += match_coverage
            results["spaced kmers & sparse"]["mp"] += mp
            # print("spaced kmers & sparse", match_coverage)
            # print_matches(all_pos_vector, "Spaced kmers")


            m,mp,sc,gaps,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 2, "minstrobes", w , w_low = w_2low, w_high = w_2high)
            results["minstrobes"][(2,15,w_2low,w_2high)]["m"] += m
            results["minstrobes"][(2,15,w_2low,w_2high)]["sc"] += sc
            results["minstrobes"][(2,15,w_2low,w_2high)]["gaps"].append(gaps)
            results["minstrobes"][(2,15,w_2low,w_2high)]["mc"] += match_coverage
            results["minstrobes"][(2,15,w_2low,w_2high)]["mp"] += mp
            # print_matches(all_pos_vector, "minstrobes2")
            # print("minstrobes2", match_coverage)
            list_for_illustration[0].append(all_pos_vector)
            # print(gaps)

            m,mp,sc,gaps,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 3, "minstrobes", w ,  w_low = w_3low, w_high = w_3high )
            results["minstrobes"][(3,10,w_3low,w_3high)]["m"] += m
            results["minstrobes"][(3,10,w_3low,w_3high)]["sc"] += sc
            results["minstrobes"][(3,10,w_3low,w_3high)]["gaps"].append(gaps)
            results["minstrobes"][(3,10,w_3low,w_3high)]["mc"] += match_coverage
            results["minstrobes"][(3,10,w_3low,w_3high)]["mp"] += mp
            # print_matches(all_pos_vector, "minstrobes3") 
            # print("minstrobes3", match_coverage)
            list_for_illustration[1].append(all_pos_vector)
            # print(gaps)

            m,mp,sc,gaps,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 2, "randstrobes", w ,  w_low = w_2low, w_high = w_2high)
            results["randstrobes"][(2,15,w_2low,w_2high)]["m"] += m
            results["randstrobes"][(2,15,w_2low,w_2high)]["sc"] += sc
            results["randstrobes"][(2,15,w_2low,w_2high)]["gaps"].append(gaps)
            results["randstrobes"][(2,15,w_2low,w_2high)]["mc"] += match_coverage
            results["randstrobes"][(2,15,w_2low,w_2high)]["mp"] += mp
            # print_matches(all_pos_vector, "randstrobes2") 
            # print("randstrobes2", match_coverage)
            list_for_illustration[2].append(all_pos_vector)
            # print(gaps)

            # Tried randstrobe n=3 with w1=17 and w2=40 and it further decreases E-size of gaps over results in paper
            # for higher mutation rates 0.05 and 0.1
            m,mp,sc,gaps,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 3, "randstrobes", w ,  w_low = w_3low, w_high = w_3high )
            results["randstrobes"][(3,10,w_3low,w_3high)]["m"] += m
            results["randstrobes"][(3,10,w_3low,w_3high)]["sc"] += sc
            results["randstrobes"][(3,10,w_3low,w_3high)]["gaps"].append(gaps)
            results["randstrobes"][(3,10,w_3low,w_3high)]["mc"] += match_coverage
            results["randstrobes"][(3,10,w_3low,w_3high)]["mp"] += mp
            # print_matches(all_pos_vector, "randstrobes3") 
            # print("randstrobes3", match_coverage)
            list_for_illustration[3].append(all_pos_vector)
            # print(gaps)

            m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size,
                                                                                 2, "multi-context",
                                                                                 w,
                                                                                 w_low=w_2low,
                                                                                 w_high=w_2high)
            results["multi-context"][(2, 15, w_2low, w_2high)]["m"] += m
            results["multi-context"][(2, 15, w_2low, w_2high)]["sc"] += sc
            results["multi-context"][(2, 15, w_2low, w_2high)]["gaps"].append(gaps)
            results["multi-context"][(2, 15, w_2low, w_2high)]["mc"] += match_coverage
            results["multi-context"][(2, 15, w_2low, w_2high)]["mp"] += mp
            # print_matches(all_pos_vector, "randstrobes3")
            # print("randstrobes3", match_coverage)
            list_for_illustration[4].append(all_pos_vector)
            # print(gaps)

            m, mp, sc, gaps, all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size,
                                                                                 3, "multi-context",
                                                                                 w, w_low=w_3low,
                                                                                 w_high=w_3high)
            results["multi-context"][(3, 10, w_3low, w_3high)]["m"] += m
            results["multi-context"][(3, 10, w_3low, w_3high)]["sc"] += sc
            results["multi-context"][(3, 10, w_3low, w_3high)]["gaps"].append(gaps)
            results["multi-context"][(3, 10, w_3low, w_3high)]["mc"] += match_coverage
            results["multi-context"][(3, 10, w_3low, w_3high)]["mp"] += mp
            # print_matches(all_pos_vector, "randstrobes3")
            # print("randstrobes3", match_coverage)
            list_for_illustration[5].append(all_pos_vector)
            # print(gaps)

            m,mp,sc,gaps,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 2, "hybridstrobes", w , w_low = w_2low, w_high = w_2high)
            results["hybridstrobes"][(2,15,w_2low,w_2high)]["m"] += m
            results["hybridstrobes"][(2,15,w_2low,w_2high)]["sc"] += sc
            results["hybridstrobes"][(2,15,w_2low,w_2high)]["gaps"].append(gaps)
            results["hybridstrobes"][(2,15,w_2low,w_2high)]["mc"] += match_coverage
            results["hybridstrobes"][(2,15,w_2low,w_2high)]["mp"] += mp
            list_for_illustration[6].append(all_pos_vector)

            m,mp,sc,gaps,all_pos_vector, match_coverage = analyze_strobemers(seq1, seq2, k_size, 3, "hybridstrobes", w , w_low = w_3low, w_high = w_3high)
            results["hybridstrobes"][(3,10,w_3low,w_3high)]["m"] += m
            results["hybridstrobes"][(3,10,w_3low,w_3high)]["sc"] += sc
            results["hybridstrobes"][(3,10,w_3low,w_3high)]["gaps"].append(gaps)
            results["hybridstrobes"][(3,10,w_3low,w_3high)]["mc"] += match_coverage
            results["hybridstrobes"][(3,10,w_3low,w_3high)]["mp"] += mp
            list_for_illustration[7].append(all_pos_vector)

            # m,sc,gaps,all_pos_vector = analyze_strobemers(seq1, seq2, 28, 4, "randstrobes", w_1 = 7, w_2 = 10, w_3 = 25 )
            # results["randstrobes"][(4,7,w_4strobe)]["m"] += m 
            # results["randstrobes"][(4,7,w_4strobe)]["sc"] += sc 
            # results["randstrobes"][(4,7,w_4strobe)]["gaps"].append(gaps) 
            # # print_matches(all_pos_vector, "randstrobes3") 
            # # list_for_illustration[4].append(all_pos_vector)
            # # print(gaps)


            # print(len(list_for_illustration))

        # plot_matches(list_for_illustration, "m", L, k_size, args.outfolder)

        plot_island_distribution2(results, mut_freq, k_size, args.outfolder)

        combined_results[mut_freq] = {}
        for protocol in results:
            non_strobemer_protocols = {"kmers & {}".format(k_size), "kmers & {}".format(k_size // 2), "kmers & {}".format(k_size * 2 // 3), "spaced kmers & sparse", "spaced kmers & dense"}
            if protocol in non_strobemer_protocols:
                flat = [g for l in results[protocol]["gaps"] for g in l]
                if flat:
                    # avg_island_len = sum(flat)/len(flat)
                    # print(protocol)
                    e_size = get_e_size(flat, L, nr_exp)
                # else:
                #     avg_island_len = 0
                res = [round(100*results[protocol]["m"]/results[protocol]["mp"], 1), 100*results[protocol]["sc"]/(L*nr_exp), 100*results[protocol]["mc"]/(L*nr_exp), e_size]
                combined_results[mut_freq][protocol] = res
                print(protocol, " & ".join([ str(round(r, 1)) for r in res]) )
            else:
                for params in results[protocol]:
                    flat = [g for l in results[protocol][params]["gaps"] for g in l]
                    if flat:
                        # avg_island_len = sum(flat)/len(flat)
                        # print(protocol, params)
                        e_size = get_e_size(flat, L, nr_exp)
                    # else:
                        # avg_island_len = 0
                    res = [round(100*results[protocol][params]["m"]/results[protocol][params]["mp"], 1), 100*results[protocol][params]["sc"]/(L*nr_exp), 100*results[protocol][params]["mc"]/(L*nr_exp), e_size]
                    combined_results[mut_freq][str(protocol) + " & " + str(params)] = res
                    print(protocol, params, " & ".join([ str(round(r, 1)) for r in res]) )
        protocols = results.keys()

    metrics_list = ["m", "sc", "mc", "E"]
    print_combined_table(combined_results, mut_freqs, metrics_list, protocols)

    # print(results)


    # # random mutation mositions
    # for mut_freq in [0.01, 0.05, 0.1]:
    #     for exp_id in range(10):
    #         seq1 = "".join([random.choice("ACGT") for i in range(L)])
    #         muts = set(random.sample(range(len(seq1)),int(L*mut_freq)))
    #         # muts = set(range(20,1000,20)) #set([20,40,60,80])
    #         # print(muts)
    #         seq2 = "".join([seq1[i] if i not in muts else random.choice(['', help_functions.reverse_complement(seq1[i]), seq1[i] + random.choice("ACGT")]) for i in range(len(seq1))])
    #         print()
    #         print("MUT FREQ:", mut_freq)
    #         positions_matching_kmers(seq1, seq2, k)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('--fasta', type=str,  default=False, help='Path to consensus fastq file(s)')
    # parser.add_argument('--k', type=int, default=13, help='Kmer size')
    # parser.add_argument('--w', type=int, default=20, help='Window size')
    parser.add_argument('--outfolder', type=str,  default="results", help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    # parser.add_argument('--pickled_subreads', type=str, help='Path to an already parsed subreads file in pickle format')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    # if len(sys.argv)==1:
    #     parser.print_help()
    #     sys.exit()

    main(args)