================================================================================
  README — Code Reproduction Guide for Reviewers
  VCM-NIC QZS Complex Operator Paper
================================================================================

This document describes how to reproduce all computational results in the
manuscript "Frequency-Domain Complex Operator-Based Two-Stage Quasi-Zero-
Stiffness Vibration Isolation System: VCM-NIC Electromechanical Coupled
Modeling, Dynamics Shaping and Stability Analysis".

================================================================================
  0. SYSTEM REQUIREMENTS
================================================================================

  - MATLAB R2020a or later (tested on R2022b, R2023a, R2024a)
  - No additional toolboxes required (core MATLAB only)
  - Operating system: Windows / Linux / macOS
  - Disk space: ~500 MB for output figures and data

================================================================================
  1. DIRECTORY STRUCTURE
================================================================================

  E:\项目二\code\
  ├── init_path.m                # Path setup (run once per session)
  ├── README_reproduction.txt    # This file
  │
  ├── lib/                       # Core computational library
  │   ├── nondim_temp2.m         # 15-DOF full HBM-AFT model (residual + Jacobian)
  │   ├── nondim_temp2_op.m      # 10-DOF operator-reduced model
  │   ├── newton.m               # Newton-Raphson solver
  │   ├── branch_follow2.m       # Arc-length continuation (15-DOF)
  │   ├── branch_follow2N.m      # Arc-length continuation (N-DOF generic)
  │   ├── FRF.m                  # FRF sweep via HBM + continuation
  │   ├── L1.m                   # Fixed-frequency force sweep
  │   ├── compute_TF_fast.m      # Force transmissibility computation
  │   ├── compute_floquet_fast.m # Floquet multiplier computation (RK4)
  │   ├── Get_Floquet_Stability.m# Batch Floquet stability evaluation
  │   ├── classify_bifurcation.m # Bifurcation type classifier
  │   ├── compute_nic_power.m    # NIC active power computation
  │   ├── arc_length_frf.m       # Reusable arc-length FRF sweep wrapper
  │   ├── branch_aux2.m          # Augmented residual for 15-DOF continuation
  │   ├── branch_aux2N.m         # Augmented residual for N-DOF continuation
  │   ├── nondim_temp2_1.m       # Variant with pure-resistor branch
  │   └── nondim_temp2_resistor_equiv.m # Equivalent damping derivation
  │
  ├── runners/                   # Top-level run scripts (START HERE)
  │   ├── verify_peak_reduction.m    # Quick: verify 62.6% peak reduction
  │   ├── run_full_optimization.m    # Full 4-phase optimization pipeline
  │   ├── run_phase_d.m              # Regenerate all journal figures
  │   └── run_phase_e.m              # Run all verification & validation tests
  │
  ├── figures/                   # Figure generation scripts
  │   ├── Generate_All_Journal_Figures.m     # Core journal EPS figures
  │   ├── Generate_Chapter4_Figures.m        # Chapter 4 operator figures (PDF)
  │   ├── Generate_All_Chapter5_6_Figures.m  # Chapters 5-6 FRF/Stability (PDF)
  │   └── Plot_Chapter4_1_ComplexOperator.m  # Operator K(Omega) detailed plot
  │
  ├── optimization/              # Optimization pipeline
  │   ├── unified_optimization.m # Main 4-phase two-level optimization
  │   ├── optimization.m         # Standalone H-infinity optimization (legacy)
  │   ├── run_optimization.m     # Entry point wrapper
  │   ├── figures/               # Output figures (optimization)
  │   └── results/               # Saved .mat results
  │
  ├── validation/                # Model validation & verification
  │   ├── Verify_Harmonic_Convergence.m      # 3-harmonic vs 5-harmonic
  │   ├── Verify_Energy_Dissipation_SCI_v2.m # Energy conservation check
  │   ├── Verify_Cubic_Truncation.m          # Cubic AFT truncation check
  │   ├── Validate_HBM_vs_ODE45.m            # HBM vs ODE45 time integration
  │   ├── Validate_01_Jacobian_Check.m       # Analytical vs FD Jacobian
  │   ├── Validate_02_Q_Elimination_Check.m  # q-coefficient elimination check
  │   ├── Validate_03_AFT_Cubic_Projection.m # AFT cubic projection accuracy
  │   ├── Validate_04_HBM_vs_ODE_FullModel.m # Full 15D HBM vs ODE
  │   ├── Validate_05_BranchFollow_Check.m   # Continuation branch validation
  │   ├── Validate_06_Full15_vs_Operator10_FRF.m # 15D vs 10D model equivalence
  │   ├── validate_resistor_equivalent_damping.m  # Resistor equivalent damping
  │   ├── verify_peak1_operator_tuning.m     # Peak 1 tuning study
  │   ├── verify_peak2_by_changing_Csh.m     # Peak 2 tuning via Csh
  │   ├── Compare_Circuit.m                  # Circuit comparison
  │   ├── compare_full_vs_op_FRF.m           # Full vs operator FRF
  │   ├── duibi.m / duibielectromech.m       # Comparison scripts
  │   ├── k_2_5dianlu.m                      # 2.5-circuit study
  │   ├── ode445_1.m                         # HBM vs ODE45 runner
  │   ├── run_op_FRF_only.m                  # Operator-only FRF
  │   └── run_resistor_only_verify.m         # Resistor-only verification
  │
  ├── stability/                 # Stability & bifurcation analysis
  │   ├── Run_Bifurcation_Classification.m   # Bifurcation type (Fold/Flip/NS)
  │   ├── Run_L1_Floquet.m                   # Force sweep + Floquet stability
  │   ├── Run_L1_Stability.m                 # Combined stability analysis
  │   ├── Run_Stability_Boundary_Map.m       # (Fw, Omega) stability boundary
  │   ├── Run_Circuit_Param_Scan.m           # Parameter sensitivity study
  │   └── figures/                           # Output stability figures
  │
  ├── operator/                  # Complex operator K(Omega) analysis
  │   ├── suanzi_op.m            # Operator-level LHS optimization
  │   ├── Parameters_suanzi.m    # Operator parameter analysis (Kr, Ki, Ceq)
  │   ├── suanzi_3ge.m           # Three-operator comparison
  │   ├── kap_e.m                # Parametric scan of kap_e
  │   ├── kap_c.m                # Parametric scan of kap_c
  │   └── sigema_of_dianlu_of_nolinear.m # Sigma parametric scan
  │
  ├── frf_sweep/                 # FRF parameter sweep experiments
  │   ├── BGModel1.m             # Wang(2017) BG model FRF reproduction
  │   ├── Run_Step1b_Linear.m    # Linear model FRF (no circuit)
  │   ├── dianlu_of_nolinear.m   # Circuit parameter FRF sweep
  │   └── nolinear.m             # Nonlinear FRF sweep (K2=0 case)
  │
  ├── boa/                       # Basin of Attraction analysis
  │   └── Run_BoA_Space.m        # BoA in (x2(0), v2(0)) plane
  │
  ├── data/                      # Saved .mat data files
  │   ├── result_resistor_multi.mat
  │   └── result_resistor_validation.mat
  │
  ├── logs/                      # Run log files
  └── output/                    # Generated figures output
      └── journal_figures/       # Final journal-ready figures


================================================================================
  2. QUICK START — REPRODUCE KEY RESULTS
================================================================================

  Step 0: Open MATLAB and navigate to the code directory
  -----------------------------------------------------------
    cd('E:\项目二\code');

  Step 1: Path setup
  -----------------------------------------------------------
    >> init_path()
    This adds all subdirectories to the MATLAB path.

  Step 2: Quick verification of the 62.6% peak reduction ( ~5 min )
  -----------------------------------------------------------
    >> run('runners/verify_peak_reduction.m')
    Expected output: "PEAK REDUCTION: 62.6x%"

  Step 3: Run the full optimization pipeline ( ~1-2 hours )
  -----------------------------------------------------------
    >> run('runners/run_full_optimization.m')
    This reproduces Figures 5-3 through 5-4 and all optimization results.

  Step 4: Regenerate all journal figures ( ~30-60 min )
  -----------------------------------------------------------
    >> run('runners/run_phase_d.m')
    Output: ~18 figures (EPS + PDF) in output/journal_figures/

  Step 5: Run all verification & validation tests ( ~1-2 hours )
  -----------------------------------------------------------
    >> run('runners/run_phase_e.m')
    Expected: All 4 verification steps pass (harmonic convergence,
              bifurcation classification, energy dissipation,
              data consistency).

================================================================================
  3. DETAILED PIPELINE DESCRIPTION
================================================================================

  3.1  Two-Level Optimization Pipeline (unified_optimization.m)
  ------------------------------------------------------------------
    The core computational contribution. Implements a 4-phase workflow:

    Phase 1 — Operator-Level LHS Rapid Screening (seconds)
      - 5000 Latin Hypercube samples over [sigma, kap_e, kap_c] ∈ [0.02, 3]³
      - Evaluates in the K(Omega) operator domain (purely algebraic)
      - Scoring function with 5 weighted objectives:
          * Low-frequency stiffness KL   (minimize)
          * Low-frequency damping CL     (maximize, target > 0.08)
          * Low-frequency inertance ML   (maximize, target > 0.30)
          * High-frequency stiffness KH  (maximize, target > 0.03)
          * High-frequency damping CH    (minimize)
          * Sharpness penalty (Ceq peak ratio > 12)
      - Selects Top-30 candidates (~20x computation reduction)

    Phase 2 — System-Level FRF Verification (minutes)
      - Full HBM + Newton for Top-30 candidates at 180 frequency points
      - Force transmissibility computed at each point
      - Continuation-based initial guess carrying for robust convergence
      - Spot-check Floquet stability at low/peak/high frequencies
      - Selects Top-8 candidates for refinement

    Phase 3 — fminsearch Local Refinement (minutes)
      - Bounded fminsearch via sigmoid parameter transform
      - Objective: J = TF_peak + stability_penalty
      - Max 80 iterations per candidate
      - Returns J-optimal and TF-optimal solutions

    Phase 4 — Multi-Candidate Fine Verification + Comparison Plot
      - 350 frequency points with full Floquet stability scan
      - Baseline computation (pure mechanical, all circuit params = 0)
      - 6-panel comparison figure (TF, Floquet, K(Omega), Ceq, Keq, Meq)
      - 2-panel candidate comparison (TF + stability)
      - NIC power assessment

    Expected Results:
      Optimal parameters:  sigma = 1.1506
                           kap_e = 1.5222
                           kap_c = 0.5743
      Peak TF reduction:   62.6%  (Baseline 0.743 → EMSD 0.278, at Fw=0.008)
      Stable frequency %:  95.4%
      NIC active power:    well within op-amp linear range (< 100 mW)

  ------------------------------------------------------------------
  3.2  Harmonic Convergence Verification (Verify_Harmonic_Convergence.m)
  ------------------------------------------------------------------
    Verifies that 3 harmonics (0, 1, 3) are sufficient by comparing
    against 5-harmonic HBM at 3 critical operating points:
      (a) Near fold bifurcation (Omega ≈ 0.8)
      (b) At peak TF response
      (c) At high excitation (Fw = 0.05)

    Expected: |A5| / |A1| < 1% at all test points, confirming that the
    5th harmonic amplitude is 2+ orders of magnitude below fundamental.
    This is because the symmetric QZS geometry yields only odd harmonics.

  ------------------------------------------------------------------
  3.3  Bifurcation Type Classification (Run_Bifurcation_Classification.m)
  ------------------------------------------------------------------
    Classifies instability mechanisms by analyzing Floquet multiplier
    crossing modes on the complex plane:
      - Real multiplier exits at (+1, 0) → Fold / Saddle-node
      - Real multiplier exits at (-1, 0) → Period-doubling / Flip
      - Complex conjugate pair exits unit circle → Neimark-Sacker
        (secondary Hopf / quasi-periodic)

    Runs force sweeps at 3 representative frequencies (Omega = 0.5, 1.0, 2.0)
    and generates:
      - Annotated amplitude-frequency figure with bifurcation type labels
      - Argand diagram showing dominant multiplier trajectories

  ------------------------------------------------------------------
  3.4  Energy Dissipation & NIC Power (Verify_Energy_Dissipation_SCI_v2.m)
  ------------------------------------------------------------------
    Verifies circuit energy conservation:
      - Computes cycle-averaged active power injected by the NIC
      - Decomposes total dissipation into passive (resistive) and
        active (NIC-powered) components
      - Confirms NIC output stays within standard op-amp linear range
        (typical: ±10 mA, ±10 V for OPA series)

  ------------------------------------------------------------------
  3.5  Model Validation Suite (validation/)
  ------------------------------------------------------------------
    Validate_01:  Analytical Jacobian vs central finite-difference
                  Criterion: relative error < 1e-6
    Validate_02:  HBM q-coefficients satisfy frequency-domain elimination
    Validate_03:  Cubic AFT projection accuracy (analytical benchmark)
    Validate_04:  15-DOF HBM vs ODE45 time-domain steady-state
    Validate_05:  Arc-length continuation branch point validation
    Validate_06:  15-DOF full model vs 10-DOF operator-reduced model
                  (verifies the operator elimination is exact)


================================================================================
  4. PHYSICAL MODEL & PARAMETER CONVENTIONS
================================================================================

  4.1  System Configuration
  ------------------------------------------------------------------
    - Primary mass m1 (upper stage): QZS mechanism with horizontal springs
    - Secondary mass m2 (lower stage): QZS mechanism + electromagnetic coupling
    - VCM (Voice Coil Motor): shunted to RLC + NIC circuit
    - NIC: synthesizes negative resistance to cancel coil resistance

  4.2  Governing Equations (Dimensionless)
  ------------------------------------------------------------------
    Mechanical (2 DOF):
      x1'' + mu*(x1'' - x2'') + F_QZS1(x1, x1') + F_em = Fw*cos(Omega*tau)
      x2'' - (x1'' - x2'') + F_QZS2(x2, x2') - F_em = 0

    Circuit (1 DOF, eliminated via K(Omega)):
      kap_e * q'' + sigma * q' + kap_c * q = theta * (x1' - x2')

    Complex operator:
      K(Omega) = theta^2 * Omega^2 / (kap_e * Omega^2 - j*sigma*Omega - kap_c)
               = Kr(Omega) + j * Ki(Omega)

    Equivalent mechanical parameters:
      Ceq(Omega) = Ki(Omega) / Omega         (equivalent damping)
      Keq(Omega) = Kr(Omega)                  (equivalent stiffness)
      Meq(Omega) = -Re{K(Omega)} / Omega^2    (equivalent mass/inertance)

  4.3  Parameter Vector sysP = [11 elements]
  ------------------------------------------------------------------
    Index  Symbol      Meaning                   Default / Range
    -----  ----------  ------------------------  ---------------------
      1    be1         Upper linear stiffness    1.0 (fixed)
      2    be2         Lower linear stiffness    alpha2 (derived)
      3    mu          Mass ratio m2/m1          0.2
      4    al1         Upper stiffness offset    alpha1 - be1 = -0.95
      5    ga1         Upper cubic coeff.        1.5
      6    ze1         Lower damping ratio       0.05
      7    lam         Coupling theta^2          0.18 (fixed)
      8    kap_e       Dimensionless inductance  variable (optimized: 1.5222)
      9    kap_c       Dim'less inv. capacitance variable (optimized: 0.5743)
     10    sigma       Dim'less resistance       variable (optimized: 1.1506)
     11    ga2         Lower cubic coeff.        1.5

  4.4  Key Dimensionless Parameters
  ------------------------------------------------------------------
    mu   = m2/m1       = 0.2         Mass ratio
    beta = k2_v/k1_v   = 2.0         Lower/upper linear stiffness ratio
    K1   = k1_h/k1_v   = 1.0         Upper horizontal/vertical spring ratio
    K2   = k2_h/k1_v   = 0.2         Lower horizontal/vertical spring ratio
    U    = L0/delta_s   = 2.0         Geometric scale
    Lg   = L2/L0        = 4/9         Rod length ratio
    v    = alpha0        = 2.5         Pre-compression parameter

    Derived:
    alpha1 = v - 2*K1*(1-Lg)/Lg       Upper linear stiffness coefficient
    alpha2 = beta - 2*K2*(1-Lg)/Lg    Lower linear stiffness coefficient
    gamma1 = K1/(U^2*Lg^3)            Upper cubic coefficient
    gamma2 = K2/(U^2*Lg^3)            Lower cubic coefficient

  4.5  HBM Harmonic Truncation
  ------------------------------------------------------------------
    Each DOF is represented by 5 HBM coefficients: [a0; a1; b1; a3; b3]
    corresponding to harmonics 0 (DC), 1 (fundamental), and 3 (3rd).
    Even harmonics are identically zero due to odd symmetry of the cubic
    Duffing nonlinearity in the QZS restoring force.

    Full model:   3 DOF × 5 coefficients = 15 unknowns + Omega/Fw = 16
    Operator model: 2 DOF × 5 coefficients = 10 unknowns + Omega/Fw = 11


================================================================================
  5. REPRODUCING SPECIFIC FIGURES
================================================================================

  Figure 2-1 (K(Omega) with frequency bands):
    >> run('figures/Plot_Chapter4_1_ComplexOperator.m')

  Figure 2-2 (Equivalent Meq, Ceq, Keq):
    >> run('figures/Generate_Chapter4_Figures.m')

  Figures 2-3, 2-4 (Parameter scans, frequency band mapping):
    >> run('operator/Parameters_suanzi.m')
    >> run('operator/kap_e.m')
    >> run('operator/kap_c.m')

  Figures 5-1, 5-2 (Circuit parameter FRF & sensitivity):
    >> run('figures/Generate_All_Chapter5_6_Figures.m')

  Figures 5-3, 5-4 (TF comparison & optimized K(Omega)):
    >> run('runners/run_full_optimization.m')

  Figure 6-1 (Stability boundary map):
    >> run('stability/Run_Stability_Boundary_Map.m')

  Figure 6-2 (Force sweep + Floquet):
    >> run('stability/Run_L1_Stability.m')

  Figure 6-3 (Pareto tradeoff):
    >> run('runners/run_full_optimization.m')  [Phase 4 output]

  Supplementary: Bifurcation Argand diagram:
    >> run('stability/Run_Bifurcation_Classification.m')


================================================================================
  6. EXPECTED COMPUTATIONAL TIME
================================================================================

  Task                                    Approximate Time
  -------------------------------------   ------------------------
  verify_peak_reduction.m                 2-5 minutes
  Verify_Harmonic_Convergence.m           5-10 minutes
  Run_Bifurcation_Classification.m        10-20 minutes
  Generate_Chapter4_Figures.m             5-10 minutes
  Generate_All_Chapter5_6_Figures.m       15-30 minutes
  Generate_All_Journal_Figures.m          10-20 minutes
  unified_optimization.m (full 4-phase)   45-90 minutes
  run_phase_d.m (all figures)             30-60 minutes
  run_phase_e.m (all tests)               60-120 minutes
  Full pipeline from scratch              ~3-5 hours

  Note: Times are approximate for a modern desktop PC (Intel i7/AMD Ryzen 7,
  16 GB RAM). The primary bottleneck is the HBM Newton solver at fine
  frequency grids (350 points × multiple parameter configurations).


================================================================================
  7. TROUBLESHOOTING
================================================================================

  Q: "Undefined function or variable" errors
  A: Run init_path() first. This must be called at the start of each
     MATLAB session.

  Q: Newton solver fails to converge at some frequency points
  A: This is expected near turning points of the FRF. The arc-length
     continuation (branch_follow2.m) is designed to handle this — use
     the run scripts rather than single-point Newton calls. The expected
     convergence rate is >95% on the Omega ∈ [0.2, 6.0] grid.

  Q: Different peak reduction % from 62.6%
  A: Ensure Fw=0.008 and at least 350 log-spaced frequency points in
     [0.2, 6.0]. Coarser grids miss the true TF peak. Verify with
     verify_peak_reduction.m which uses the exact optimized parameters.

  Q: Out of memory errors
  A: Close other MATLAB figures. The largest memory usage occurs during
     the 5000-sample LHS screening in Phase 1 (~1-2 GB for 350×5000 matrix).

  Q: Figures don't match paper exactly
  A: Ensure the random seed is set (rng(42) in unified_optimization.m).
     The LHS sampling is deterministic with a fixed seed.

================================================================================
  8. CODE CITATION & LICENSE
================================================================================

  If you use this code in your research, please cite the corresponding paper.
  The code is provided for research reproducibility purposes.

  Key references:
    [1] Wang et al. (2017) — BG model baseline
    [2] The current manuscript — Complex operator K(Omega) framework

================================================================================
  Last updated: 2026-05-29
  Contact: [corresponding author email]
================================================================================
