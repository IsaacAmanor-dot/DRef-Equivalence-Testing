# DRef Equivalence Testing

This workflow tests whether isolating Dynamic Referencing (DRef) from `DN_Build` into a dedicated `Dynamic_Reference` manager changes the behavior of DOCK6.9 De Novo molecular generation.

The inherited and isolated implementations are run using identical receptor files, starting anchor, fragment libraries, scoring parameters, and De Novo settings. The resulting molecular structures and growth outputs are then compared to determine behavioral equivalence.

The initial test uses the `1EQH` receptor system and anchor 1. Calculations are run on the `rn-long-40core` partition.

## Workflow

Run the scripts in the following order.

### 1. Check configuration

```bash
bash 000.config.sh
```

### 2. Set up the test

```bash
bash 001.setup_test.sh
```

### 3. Run the calculations

```bash
bash 002.run_test.sh
```

### 4. Compare the results

```bash
python3 003.compare_results.py
```

The comparison report is written to:

```text
comparison_results/equivalence_report.txt
```
