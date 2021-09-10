#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""Main classification workflow.

Notes
-----
Only in this script can functions directly interface with the user by screen
output (via `click`) and file input/output, except for raising errors.
"""

from os import makedirs
from os.path import join, basename, isfile, isdir
from collections import deque, defaultdict
from functools import partial, lru_cache
from typing import Tuple
from pathos.multiprocessing import ProcessingPool as Pool
from multiprocessing import cpu_count
import click
from .util import (
    update_dict, allkeys, sum_dict, scale_factor, scale_dict, round_dict)
from .file import (
    openzip, readzip, path2stem, stem2rank, read_ids, id2file_from_dir,
    id2file_from_map, read_map_uniq, read_map_1st, write_readmap)
from .align import plain_mapper
from .classify import (
    assign_none, assign_free, assign_rank, counter, counter_size,
    counter_strat, counter_size_strat)
from .tree import (
    read_names, read_nodes, read_lineage, read_newick, read_columns,
    fill_root)
from .ordinal import (
    ordinal_mapper, read_gene_coords, whether_prefix, calc_gene_lens)
from .table import prep_table, write_table


def classify(mapper:  object,
             files:     list or dict,
             samples:   list = None,
             fmt:        str = None,
             demux:     bool = None,
             trimsub:    str = None,
             tree:      dict = None,
             rankdic:   dict = None,
             namedic:   dict = None,
             root:       str = None,
             ranks:      str = None,
             rank2dir:  dict = None,
             outzip:     str = None,
             uniq:      bool = False,
             major:      int = None,
             above:     bool = False,
             subok:     bool = False,
             sizes:     dict = None,
             unasgd:    bool = False,
             stratmap:  dict = None,
             chunk:      int = None,
             cache:      int = 1024,
             zippers:   dict = None) -> dict:
    """Core of the classification workflow.

    Parameters
    ----------
    mapper : object
        Mapping module (Plain or Ordinal).
    files : list or dict
        Paths to input alignment files, if multiplexed, or dictionary of file
        paths to sample IDs, if per-sample.
    samples : list of str, optional
        Sample ID list to include.
    fmt : str, optional
        Format of input alignment file. Options:
        - 'b6o': BLAST tabular format.
        - 'sam': SAM format.
        - 'map': Simple map of query <tab> subject.
        If None, program will automatically infer from file content.
    demux : bool, optional
        Whether perform demultiplexing.
    trimsub : str, optional
        Trim subject IDs at the last given delimiter.
    tree : dict, optional
        Taxonomic tree.
    rankdic : dict, optional
        Rank dictionary.
    namedic : dict, optional
        Taxon name dictionary.
    root : str, optional
        Root identifier.
    ranks: list of str, optional
        List of ranks at each of which sequences are to be classified. Can also
        be "none" to omit classification (simply report subject IDs) or "free"
        to perform free-rank classification (LCA of subjects regardless of rank
        will be reported).
    rank2dir : dict, otional
        Write classification map per rank to directory.
    outzip : str, optional
        Output read map compression method (gz, bz2, xz or None).
    uniq : bool, optional
        Assignment must be unique. Otherwise, report all possible assignments
        and normalize counts (for none- and fixed-rank assignments).
    major : int, optional
        In given-rank classification, perform majority-rule assignment based on
        this percentage threshold. Range: [51, 99].
    above : bool, optional
        Allow assigning to a classification unit higher than given rank.
    subok : bool, optional
        In free-rank classification, allow assigning sequences to their direct
        subjects instead of higher classification units, if applicable.
    sizes : dict, optional
        Subject size dictionary.
    unasgd : bool, optional
        Report unassigned sequences.
    stratmap : dict, optional
        Map of sample ID to stratification file.
    chunk : int, optional
        Number of lines per chunk to read from alignment file.
    cache : int, optional
        LRU cache size for classification results at each rank.
    zippers : dict, optional
        External compression programs.

    Returns
    -------
    dict of dict
        Per-rank profiles generated from classification.

    Notes
    -----
    Subject(s) of each query are sorted and converted into a tuple, which is
    hashable, a property necessary for subsequent assignment result caching.
    """
    data = {x: {} for x in ranks}

    # assigners for each rank
    assigners = {}

    # assignment parameters
    kwargs = {'assigners': assigners, 'cache': cache, 'tree': tree, 'rankdic':
              rankdic, 'namedic': namedic, 'root':  root, 'uniq': uniq,
              'major': major and major / 100, 'above': above, 'subok': subok,
              'sizes': sizes, 'unasgd': unasgd, 'rank2dir': rank2dir,
              'outzip': outzip if outzip != 'none' else None}

    # current sample Id
    csample = False

    # parse input alignment file(s) and generate profile(s)
    def classify_file_mp(fp):
        click.echo(f'Parsing alignment file {basename(fp)} ', nl=False)
        # read alignment file into query-to-subject(s) map
        with readzip(fp, zippers) as fh:
            # query and progress counters
            nqry, nstep = 0, -1	
            # parse alignment file by chunk
            for qryque, subque in mapper(fh, fmt=fmt, n=chunk):
                nqry += len(qryque)
                # (optional) strip indices and generate tuples
                subque = deque(map(tuple, map(sorted, strip_suffix(
                    subque, trimsub) if trimsub else subque)))

                # (optional) demultiplex and generate per-sample maps
                rmaps = demultiplex(qryque, subque, samples) if demux else {
                    files[fp] if files else None: (qryque, subque)}

                # assign reads at each rank
                for sample, rmap in rmaps.items():

                    # (optional) read strata of current sample into cache
                    if stratmap and sample != csample:
                        kwargs['strata'] = read_strata(
                            stratmap[sample], zippers)
                        csample = sample

                    # call assignment workflow for each rank
                    for rank in ranks:
                        assign_readmap(*rmap, data, rank, sample, **kwargs)

                # show progress
                istep = nqry // 1000000 - nstep
                if istep:
                    click.echo('.' * istep, nl=False)
                    nstep += istep

        click.echo(' Done.')
        click.echo(f'  Number of sequences classified: {nqry}.')
        return data
    
    #create n processes where n is max 12 cores 
    num_processes = cpu_count()
    if num_processes >= 12:
         classify_pool = Pool(12)
    else:
         classify_pool = Pool(num_processes)
    classify_results = classify_pool.map(classify_file_mp, sorted(files))
    classify_pool.close()
    classify_pool.join()
    
    #aggregate the results from all the processes into the "data" nested dict object
    for results in classify_results:
        for rank,file in results.items():
            data[rank].update(file)

    click.echo('Classification completed.')
    return data
