# DOCK6.14 Dynamic Referencing Equivalence Testing

## Overview

This project tests whether extracting the Dynamic Referencing (DRef) algorithm into `dynamic_reference.cpp` and `dynamic_reference.h` preserves the computational behavior of DOCK6.14.

The original (baseline) and modularized implementations are tested using identical inputs, scoring parameters, and random seeds. Molecular structures, scores, partial references, and growth-tree outputs are compared to evaluate equivalence.

The test uses system **1O86** with HMS-based Dynamic Referencing.

Source repository: https://github.com/rizzolab/dock6_dev

Branch: `dref-modularization`

## Workflow

Run the following scripts in order.

### 1. Check configuration

```bash
bash 000.config.sh
```

### 2. Set up the test

```bash
bash 001.setup_test.sh /gpfs/scratch/iamanor/Collaborations/Carissa/16_10_2026_Dref_and_Wref/DRef-and-WRef-One-System-Multiple-Anchor-Experiment-for-DOCK6.14/DRef_6.14_1EVE_20_Anchors/1EVE/anc_1/30_10_2026_Testing/dn_generic.in
```

This prepares identical inputs for both implementations using system 1O86.

### 3. Submit calculations

```bash
bash 002.run_test.sh
```

Both implementations run sequentially on the `gpu-long` partition with a 2-day walltime.

### 4. Compare results

```bash
python3 003.compare_results.py
```

### 5. Analyze equivalence

```bash
python3 004.equivalence_analysis.py
```

Scripts `003` and `004` must be adapted to DOCK6.14 outputs before use.

## Results

Calculation outputs are stored in:

- `baseline/system/`
- `modularized/system/`

Comparison reports, statistics, and plots are saved in:

`comparison_results/`
