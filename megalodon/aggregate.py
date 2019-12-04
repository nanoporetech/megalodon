#!/usr/bin/env python3
import os
import sys
import queue
from time import sleep
import multiprocessing as mp

from tqdm import tqdm

from megalodon import logging, mods, variants, megalodon_helper as mh

_DO_PROFILE_AGG_MOD = False
_DO_PROFILE_GET_MODS = False
_DO_PROFILE_AGG_FILLER = False
_DO_PROF = _DO_PROFILE_AGG_MOD or _DO_PROFILE_AGG_FILLER or _DO_PROFILE_GET_MODS
_N_MOD_PROF = 200000


############################################
##### Aggregate Variants and Mod Stats #####
############################################

def _agg_vars_worker(
        locs_q, var_stats_q, var_prog_q, vars_db_fn, write_vcf_lp,
        het_factors, call_mode, valid_read_ids):
    agg_vars = variants.AggVars(vars_db_fn, write_vcf_lp)

    while True:
        try:
            var_loc = locs_q.get(block=False)
        except queue.Empty:
            sleep(0.001)
            continue
        if var_loc is None:
            break

        try:
            var_var = agg_vars.compute_var_stats(
                var_loc, het_factors, call_mode, valid_read_ids)
            var_stats_q.put(var_var)
        except mh.MegaError:
            # something not right with the stats at this loc
            pass
        var_prog_q.put(1)

    return

def _get_var_stats_queue(
        var_stats_q, var_conn, out_dir, ref_names_and_lens, out_suffix,
        write_vcf_lp):
    agg_var_fn = mh.get_megalodon_fn(out_dir, mh.VAR_NAME)
    if out_suffix is not None:
        base_fn, fn_ext = os.path.splitext(agg_var_fn)
        agg_var_fn = base_fn + '.' + out_suffix + fn_ext
    agg_var_fp = variants.VcfWriter(
        agg_var_fn, 'w', ref_names_and_lens=ref_names_and_lens,
        write_vcf_lp=write_vcf_lp)

    while True:
        try:
            var_var = var_stats_q.get(block=False)
            if var_var is None: continue
            agg_var_fp.write_variant(var_var)
        except queue.Empty:
            if var_conn.poll():
                break
            sleep(0.001)
            continue

    while not var_stats_q.empty():
        var_var = var_stats_q.get(block=False)
        agg_var_fp.write_variant(var_var)

    agg_var_fp.close()

    return

def _agg_mods_worker(
        pos_q, mod_stats_q, mod_prog_q, mods_db_fn, mod_agg_info,
        valid_read_ids, write_mod_lp):
    # functions for profiling purposes
    def get_pos_data():
        return pos_q.get(block=False)
    def put_mod_site(mod_site):
        mod_stats_q.put(mod_site)
        return
    def do_sleep():
        sleep(0.0001)
        return

    agg_mods = mods.AggMods(mods_db_fn, mod_agg_info, write_mod_lp)

    while True:
        try:
            pos_data = get_pos_data()
        except queue.Empty:
            do_sleep()
            continue
        if pos_data is None:
            break

        try:
            mod_site = agg_mods.compute_mod_stats(
                pos_data, valid_read_ids=valid_read_ids)
            put_mod_site(mod_site)
        except mh.MegaError:
            # no valid reads cover location
            pass
        mod_prog_q.put(1)

    return

if _DO_PROFILE_AGG_MOD:
    _agg_mods_wrapper = _agg_mods_worker
    def _agg_mods_worker(*args):
        import cProfile
        cProfile.runctx('_agg_mods_wrapper(*args)', globals(), locals(),
                        filename='aggregate_mods.prof')
        return

def _get_mod_stats_queue(
        mod_stats_q, mod_conn, out_dir, mod_names, ref_names_and_lens,
        out_suffix, write_mod_lp, mod_output_fmts):
    def get_mod_site():
        # function for profiling purposes
        return mod_stats_q.get(block=False)
    def do_sleep():
        # function for profiling purposes
        sleep(0.001)
        return

    agg_mod_bn = mh.get_megalodon_fn(out_dir, mh.MOD_NAME)
    if out_suffix is not None:
        agg_mod_bn += '.' + out_suffix
    agg_mod_fps = []
    if mh.MOD_BEDMETHYL_NAME in mod_output_fmts:
        agg_mod_fps.append(mods.ModBedMethylWriter(
            agg_mod_bn, mod_names, 'w'))
    if mh.MOD_VCF_NAME in mod_output_fmts:
        agg_mod_fps.append(mods.ModVcfWriter(
            agg_mod_bn, mod_names, 'w', ref_names_and_lens=ref_names_and_lens,
            write_mod_lp=write_mod_lp))
    if mh.MOD_WIG_NAME in mod_output_fmts:
        agg_mod_fps.append(mods.ModWigWriter(
            agg_mod_bn, mod_names, 'w'))

    while True:
        try:
            mod_site = get_mod_site()
            for agg_mod_fp in agg_mod_fps:
                agg_mod_fp.write_mod_site(mod_site)
        except queue.Empty:
            if mod_conn.poll():
                break
            do_sleep()
            continue

    while not mod_stats_q.empty():
        mod_site = mod_stats_q.get(block=False)
        for agg_mod_fp in agg_mod_fps:
            agg_mod_fp.write_mod_site(mod_site)
    for agg_mod_fp in agg_mod_fps:
        agg_mod_fp.close()

    return

if _DO_PROFILE_GET_MODS:
    _get_mod_stats_queue_wrapper = _get_mod_stats_queue
    def _get_mod_stats_queue(*args):
        import cProfile
        cProfile.runctx('_get_mod_stats_queue_wrapper(*args)',
                        globals(), locals(),
                        filename='get_mods_queue.prof')
        return

def _agg_prog_worker(
        var_prog_q, mod_prog_q, num_vars, num_mods, prog_conn,
        suppress_progress):
    var_bar, mod_bar = None, None
    if num_vars > 0:
        if num_mods > 0 and not suppress_progress:
            mod_bar = tqdm(desc='Mods', unit=' sites', total=num_mods,
                           position=1, smoothing=0, dynamic_ncols=True)
            var_bar = tqdm(desc='Variants', unit=' sites', total=num_vars,
                           position=0, smoothing=0, dynamic_ncols=True)
        elif not suppress_progress:
            var_bar = tqdm(desc='Variants', unit=' sites', total=num_vars,
                           position=0, smoothing=0, dynamic_ncols=True)
    elif num_mods > 0 and  not suppress_progress:
        mod_bar = tqdm(desc='Mods', unit=' sites', total=num_mods,
                       position=0, smoothing=0, dynamic_ncols=True)

    while True:
        try:
            var_prog_q.get(block=False)
            if not suppress_progress:
                if var_bar is not None: var_bar.update(1)
                if mod_bar is not None: mod_bar.update(0)
        except queue.Empty:
            try:
                mod_prog_q.get(block=False)
                if not suppress_progress:
                    if var_bar is not None: var_bar.update(0)
                    if mod_bar is not None: mod_bar.update(1)
            except queue.Empty:
                sleep(0.001)
                if prog_conn.poll():
                    break
                continue

    while not var_prog_q.empty():
        var_prog_q.get(block=False)
        if not suppress_progress: var_bar.update(1)
    while not mod_prog_q.empty():
        mod_prog_q.get(block=False)
        if not suppress_progress: mod_bar.update(1)
    if var_bar is not None:
        var_bar.close()
    if mod_bar is not None:
        mod_bar.close()
    if num_mods > 0 and num_vars > 0 and not suppress_progress:
        sys.stderr.write('\n\n')

    return

def _fill_locs_queue(locs_q, db_fn, agg_class, num_ps, limit=None):
    agg_db = agg_class(db_fn, load_in_mem_indices=False)
    for i, loc in enumerate(agg_db.iter_uniq()):
        locs_q.put(loc)
        if limit is not None and i >= limit: break
    for _ in range(num_ps):
        locs_q.put(None)

    return

if _DO_PROFILE_AGG_FILLER:
    _fill_locs_queue_wrapper = _fill_locs_queue
    def _fill_locs_queue(*args):
        import cProfile
        cProfile.runctx('_fill_locs_queue_wrapper(*args)', globals(), locals(),
                        filename='agg_fill_locs.prof')
        return

def aggregate_stats(
        outputs, out_dir, num_ps, write_vcf_lp, het_factors, call_mode,
        mod_names, mod_agg_info, write_mod_lp, mod_output_fmts,
        suppress_progress, valid_read_ids=None, out_suffix=None):
    if mh.VAR_NAME in outputs and mh.MOD_NAME in outputs:
        num_ps = max(num_ps // 2, 1)

    logger = logging.get_logger('agg')
    num_vars, num_mods, var_prog_q, mod_prog_q = (
        0, 0, queue.Queue(), queue.Queue())
    if mh.VAR_NAME in outputs:
        vars_db_fn = mh.get_megalodon_fn(out_dir, mh.PR_VAR_NAME)
        agg_vars = variants.AggVars(
            vars_db_fn, load_in_mem_indices=False)
        num_vars = agg_vars.num_uniq()
        ref_names_and_lens = agg_vars.vars_db.get_all_chrm_and_lens()
        agg_vars.close()
        logger.info('Spawning variant aggregation processes.')
        # create process to collect var stats from workers
        var_stats_q, var_stats_p, main_var_stats_conn = mh.create_getter_q(
            _get_var_stats_queue, (
                out_dir, ref_names_and_lens, out_suffix, write_vcf_lp))
        # create process to fill variant locs queue
        var_filler_q = mp.Queue(maxsize=mh._MAX_QUEUE_SIZE)
        var_filler_p = mp.Process(
            target=_fill_locs_queue,
            args=(var_filler_q, vars_db_fn, variants.AggVars, num_ps),
            daemon=True)
        var_filler_p.start()
        # create worker processes to aggregate variants
        var_prog_q = mp.Queue(maxsize=mh._MAX_QUEUE_SIZE)
        agg_vars_ps = []
        for _ in range(num_ps):
            p = mp.Process(
                target=_agg_vars_worker,
                args=(var_filler_q, var_stats_q, var_prog_q, vars_db_fn,
                      write_vcf_lp, het_factors, call_mode, valid_read_ids),
                daemon=True)
            p.start()
            agg_vars_ps.append(p)

    if mh.MOD_NAME in outputs:
        mods_db_fn = mh.get_megalodon_fn(out_dir, mh.PR_MOD_NAME)
        agg_mods = mods.AggMods(mods_db_fn, load_in_mem_indices=False)
        num_mods = agg_mods.num_uniq()
        ref_names_and_lens = agg_mods.mods_db.get_all_chrm_and_lens()
        agg_mods.close()
        logger.info('Spawning modified base aggregation processes.')
        # create process to collect mods stats from workers
        mod_stats_q, mod_stats_p, main_mod_stats_conn = mh.create_getter_q(
            _get_mod_stats_queue, (
                out_dir, mod_names, ref_names_and_lens, out_suffix,
                write_mod_lp, mod_output_fmts))
        # create process to fill mod locs queue
        mod_filler_q = mp.Queue(maxsize=mh._MAX_QUEUE_SIZE)
        mod_fill_limit = _N_MOD_PROF if _DO_PROF else None
        mod_filler_p = mp.Process(
            target=_fill_locs_queue,
            args=(mod_filler_q, mods_db_fn, mods.AggMods, num_ps,
                  mod_fill_limit), daemon=True)
        mod_filler_p.start()
        # create worker processes to aggregate mods
        mod_prog_q = mp.Queue(maxsize=mh._MAX_QUEUE_SIZE)
        agg_mods_ps = []
        for _ in range(num_ps):
            p = mp.Process(
                target=_agg_mods_worker,
                args=(mod_filler_q, mod_stats_q, mod_prog_q, mods_db_fn,
                      mod_agg_info, valid_read_ids, write_mod_lp),
                daemon=True)
            p.start()
            agg_mods_ps.append(p)

    # create progress process
    logger.info((
        'Aggregating {} variants and {} modified base sites over reads.\n\t\t' +
        'NOTE: If this step is very slow, ensure the output directory is ' +
        'located on a fast read disk (e.g. local SSD). Aggregation can be ' +
        'restarted using the megalodon/scripts/run_aggregation.py ' +
        'script.').format(num_vars, num_mods))
    main_prog_conn, prog_conn = mp.Pipe()
    prog_p = mp.Process(
        target=_agg_prog_worker,
        args=(var_prog_q, mod_prog_q, num_vars, num_mods, prog_conn,
              suppress_progress),
        daemon=True)
    prog_p.start()

    # join filler processes first
    if mh.VAR_NAME in outputs:
        var_filler_p.join()
        for agg_vars_p in agg_vars_ps:
            agg_vars_p.join()
        # send to conn
        if var_stats_p.is_alive():
            main_var_stats_conn.send(True)
        var_stats_p.join()
    if mh.MOD_NAME in outputs:
        for agg_mods_p in agg_mods_ps:
            agg_mods_p.join()
        if mod_stats_p.is_alive():
            main_mod_stats_conn.send(True)
        mod_stats_p.join()
    if prog_p.is_alive():
        main_prog_conn.send(True)
        prog_p.join()

    return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
