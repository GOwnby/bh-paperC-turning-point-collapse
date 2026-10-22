#!/bin/bash
# Regenerates every number of the Lagrangian (Misner-Sharp) verification study.  Two independent job streams run in parallel (2 cores).
# usage:  bash run_all.sh            (from any directory)   -- target wall clock < 15 min on 2 idle cores; logs in logs/, tables in results/
# Requirements: python3 with numpy, scipy (no other packages).  The report's code is used only through the unmodified copies in ref/.
set -u
cd "$(dirname "$0")/src"
mkdir -p ../logs ../results
export PYTHONWARNINGS=ignore   # numpy overflow warnings of the deliberately destroyed runs (pressure-floor test) would otherwise put machine paths into the logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
T0=$(date +%s)
PY="nice -n 10 python3 -u"
: > ../logs/timing.log
# step <logfile-name> <script> [args]: runs one script, logs its wall time
step() { local name=$1; shift; local t=$(date +%s); $PY "$@" > ../logs/$name.log 2>&1; local rc=$?; echo "$name: $(( $(date +%s) - t )) s (exit $rc)" >> ../logs/timing.log; }
# the verified copies of the report's routines must be unmodified
( cd ../ref && sha256sum -c SHA256SUMS.copy > ../logs/ref_checksums.log 2>&1 ) || echo "WARNING: ref/ checksum mismatch, see logs/ref_checksums.log"
rm -f ../results/eos_table_*.npz ../results/*.json ../results/RESULTS.md
step eos_tables eosmodel.py          # builds the EOS tables from the report's routines
step structure_tov structure.py      # TOV solutions (M, R, compactness)

streamA() {
  step verify_ms_equations     verify_ms_equations.py        # (1) Einstein-tensor identities, negative control
  step control_star_truncation control_star_truncation.py    # the report's rho_c = 1e11 control was truncated
  step test_a_static           test_a_static.py              # (a) static equilibrium, N = 50..800
  step test_b_linear           test_b_linear.py              # (b) linear modes vs Galerkin and independent Sturm-Liouville; growth rate
  step repro_spline            repro_spline.py               # (iii) the spline re-creation: reproduction and diagnosis
  step test_e_conservation     test_e_conservation.py        # (e) mass / energy conservation
  step test_f_forced           test_f_forced.py              # (f) distance of the report's data from GR equilibrium
  step test_i_dust             test_i_dust.py                # (i) exact Oppenheimer-Snyder / LTB dust solutions, order of the scheme, Gamma iteration
  step test_j_gamma_size       test_j_gamma_size.py          # (j) size of the time variation of Gamma^2
}
streamB() {
  step test_c_converge_pure    test_c_converge.py pure       # (c) convergence, N = 100..1600, dt/2, AV
  step test_c_converge_neu     test_c_converge.py neu
  step test_g_long             test_g_long.py                # (c') inviscid to T = 150 L/c; exact-equilibrium control
  step test_h_collapse         test_h_collapse.py            # (c'') viscous continuation: converged collapse metrics, time-step checks
  step test_d_edge             test_d_edge.py                # (d) vacuum-edge sensitivity
}
streamA & PA=$!
streamB & PB=$!
wait $PA $PB
step make_tables make_tables.py
echo "total wall time: $(( $(date +%s) - T0 )) s" | tee ../logs/wall_time.log
echo "--- per-step wall times"; cat ../logs/timing.log
echo "--- PASS/FAIL lines"
grep -h "PASS\|FAIL" ../logs/*.log | sed 's/^ *//' | sort | uniq -c | sort -rn | head -60
