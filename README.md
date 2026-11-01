# DOCK6.14 Dynamic Referencing Equivalence Testing

This workflow evaluates whether the isolated Dynamic Referencing implementation reproduces the molecular results of the inherited implementation in DOCK6.14.

## Configuration

Edit `000.config.sh` to configure the number of anchors and SLURM resources.

Specify the receptor systems in `LISTOFSYSTEMS.txt`, one system per line.

## Workflow

Generate identical inputs for inherited and isolated calculations:

```bash
bash 001.setup_test.sh
```

Submit two independent SLURM jobs on `rn-long-40core`:

```bash
bash 002.submit_jobs.sh
```

After both jobs finish, compare their results:

```bash
python3 003.compare_results.py
```

## Results

The calculations are organized under `runs/inherited/` and `runs/isolated/` by system and anchor.

The analysis compares molecular structures, coordinates, scores, growth-layer populations, and growth trees.

Reports are saved under `comparison_results/`.
