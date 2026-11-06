# Running fish-ib2d on the ICF cluster

How this project is set up on the University of Edinburgh's Informatics Compute Facility (ICF), and how to run simulations there. It is written for a person working by hand and for a coding agent. Agents: read [Rules for coding agents](#rules-for-coding-agents) before running anything.

Official documentation: [research cluster](https://computing.help.inf.ed.ac.uk/research-cluster), [cluster tips](https://computing.help.inf.ed.ac.uk/cluster-tips), [Slurm](https://computing.help.inf.ed.ac.uk/slurm).

## 1. The cluster at a glance

| | |
|---|---|
| Head nodes | `icf.inf.ed.ac.uk` (= `hastings`, also reachable as `mlp.inf.ed.ac.uk`) and `icf2.inf.ed.ac.uk` (= `stanger`). Ubuntu 24.04. |
| Access | From the School's OpenVPN, a DICE machine, or through the SSH gateways ([external login](https://computing.help.inf.ed.ac.uk/external-login)). The SSH config alias on the user's Mac is `icf`. |
| Scheduler | Slurm. Slurm account: `research`. |
| Head node rule | **Only edit files and submit jobs on the head node.** Simulations, geometry builds and anything else heavy run on compute nodes, through `sbatch` or `srun`. |
| User | `s2908217`, home `/home/s2908217` |

### Partitions

| Partition | Time limit | Use it for |
|---|---|---|
| `ICF-Free` | 2 days | **Batch jobs: always set this one.** Research accounts only. |
| `Open-Research` | 2 days | Alternative to `ICF-Free`, also research accounts only. |
| `Interactive` | 4 h | Short `srun` tests. Keep them under an hour and exit as soon as you're done. |
| `Teaching` | 2 days | **The default if you don't set a partition.** Don't use it for research jobs. |

Every compute node is a GPU machine, and IB2d uses **1 CPU core and no GPU**. So request `--cpus-per-task=1`, never `--gres`, and cap how many jobs run at once. Taking many CPUs on a GPU node can stop other people's GPU jobs from starting there.

## 2. Storage

| Location | What it is | What goes there |
|---|---|---|
| `~/finray_fsi/fish-ib2d/` | The git repo, on the shared (Lustre) filesystem | Code only |
| `~/finray_fsi/sim_data/fish-ib2d/FinRay/<sweep>/` | Shared filesystem, outside the repo | Cluster sweep results (snapshot, logs, run folders) |
| `/disk/scratch/$USER/...` **on a compute node** | Fast local disk inside that node; cleared from time to time | Temporary output while a job runs. The job scripts handle it. |
| `cases/FinRay/results_approach2/` | In the repo, ignored by git | Runs started by hand with `run_case.py` (e.g. tests) |

- **Nothing on the cluster is backed up**, home included. Copy results you can't lose elsewhere (your Mac, AFS, DataStore).
- There's no per-user quota, but the space is shared (about 100 TB free across everyone), so delete what you no longer need.
- The shared filesystem is slow with many small files. Pack files into a `.tgz` before copying, and expect `du` or `find` over results folders to take minutes.
- `/disk/scratch` on the head node is not the scratch disk of the compute nodes. Each node has its own.

## 3. One-time setup (already done)

All of this already exists. Redo it only on a fresh account.

```bash
git clone --recurse-submodules git@github.com:constancefr/fish-ib2d.git ~/finray_fsi/fish-ib2d
#   (already cloned? fetch the solver with: git -C ~/finray_fsi/fish-ib2d submodule update --init)
conda create -n ib2d python=3.11 -y
conda activate ib2d
pip install -r ~/finray_fsi/fish-ib2d/requirements.txt     # includes meshpy
```

- **Conda environment:** `ib2d` lives in `~/.conda/envs/ib2d`, so compute nodes can use it too. In job scripts, activate it with `source /opt/conda/etc/profile.d/conda.sh && conda activate ib2d`.
- **GitHub:** pushes go over SSH, using the key `~/.ssh/id_ed25519` (registered on the GitHub account). The remote is `git@github.com:constancefr/fish-ib2d.git`.
- **VS Code (on the Mac):** connect with Remote-SSH to host `icf`, open the project folder `~/finray_fsi`, and select the `ib2d` interpreter.
- **Not installed:** VisIt (see [§6](#6-movies-and-analysis)).

## 4. Running a FinRay sweep (the normal way)

The scripts are in `cases/FinRay/slurm/`. Each case of a sweep (a function in `cases/FinRay/sweeps.py`) becomes one task of a Slurm job array. Each task runs on a compute node, in a private copy of the code on that node's `/disk/scratch`. When it ends, the task copies its run folder back to `~/finray_fsi/sim_data`.

### 4.1 Set up the sweep

1. Make sure the repo is on the code you want: `cd ~/finray_fsi/fish-ib2d && git status && git log --oneline -1`.
2. Set the base parameters in `cases/FinRay/case_params.json`: the regime (`actuation.FREQUENCY`, `AMPLITUDE`), the run length (`run_length`; `input2d.Tfinal: null` means Tfinal = cycles / frequency), and `input2d` values such as `dt`, which `run_case.py` copies into `input2d`.
3. Pick the sweep, or add one, in `cases/FinRay/sweeps.py`: one function per sweep, named by what it is for (e.g. `stageC0`). Existing: `base` (the base case alone), `smoke` (the base case for 0.05 s), `density`.
4. Check the list of cases without running anything:
   ```bash
   cd ~/finray_fsi/fish-ib2d/cases/FinRay && python run_case.py --sweep <sweep_def> --list
   ```

Uncommitted changes are fine. They're recorded with every run (see 4.4). Still, prefer committing the code you run.

### 4.2 Submit (on the head node)

```bash
~/finray_fsi/fish-ib2d/cases/FinRay/slurm/submit_sweep.sh <sweep_name> <sweep_def>
```

- **`<sweep_def>`** is the function in `sweeps.py` (step 3). It's saved as `sweep_def.txt` in the sweep folder and in every run's `run_info.json` (`inputs.sweep`).
- **`<sweep_name>` must be new.** The script refuses a name that already exists in `~/finray_fsi/sim_data/fish-ib2d/FinRay/`, so earlier results are never overwritten. Use a name that says what the sweep is, e.g. `mult0-256_Nx32_T15`.
- **It records the code at submit time.** It copies the code to `<sweep>/snapshot/` and saves the git state there. Editing the repo afterwards doesn't affect the sweep, even if it's still queued.
- **It submits two jobs.**
  - The **array job**: one task per case, at most 8 running at once.
  - The **index job**: it waits until every task has ended, finished or not, then builds `index.csv` and `index_full.csv`.
- **Optional settings** go before the command, e.g. `MAX_PARALLEL=4 ~/finray_fsi/fish-ib2d/cases/FinRay/slurm/submit_sweep.sh <name> <sweep_def>`:
  - `MAX_PARALLEL`: the number of cases running at once (default 8);
  - `CONDA_ENV`: the conda environment the jobs use (default `ib2d`);
  - `DATA_ROOT`: where sweeps are stored (default `~/finray_fsi/sim_data/fish-ib2d/FinRay`).

**Resources for each task** are set in `sweep_array.sbatch`: `ICF-Free`, 1 CPU, 4 GB of memory, a 2-hour limit. They haven't been measured on this cluster yet. After the first real sweep, check the actual usage with `sacct` (4.3: `Elapsed` = time, `MaxRSS` = peak memory) and set `--mem` and `--time` to about 1.5× that. Smaller requests start sooner.

### 4.3 Monitor and control

```bash
squeue -u $USER                                   # PD = waiting, R = running
tail -f ~/finray_fsi/sim_data/fish-ib2d/FinRay/<sweep>/logs/<jobid>_<task>.out
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,ReqMem,Timelimit,ExitCode
                                                  # per task: state, time and peak memory used (seff isn't installed)
scontrol update JobId=<jobid> ArrayTaskThrottle=4 # change how many run at once
scancel <jobid>                                   # cancel the whole array (or <jobid>_<task> for one task)
```

### 4.4 What you get

```
~/finray_fsi/sim_data/fish-ib2d/FinRay/<sweep>/
├── snapshot/                     code as submitted (+ code_provenance.json / .diff: git state at submit)
├── cases.txt                     array task number -> case name
├── sweep_def.txt                 the sweeps.py function used
├── logs/                         <jobid>_<task>.out per task, index_<jobid>.out
└── results/
    ├── index.csv                 one row per run, chosen columns (INDEX_COLUMNS in run_case.py)
    ├── index_full.csv            one row per run, every field
    ├── <case>_<YYYYmmdd-HHMMSS>/ one folder per run:
    │   ├── viz_IB2d/             .vtk frames (for VisIt / analysis)
    │   ├── run_info.json         inputs, mesh stats, timing, outcome, analysis (postprocess.py), code version, machine
    │   ├── tip_trajectory.csv, fluid_energy.csv   time series from postprocess.py
    │   ├── case_params.json, input2d   the exact inputs used
    │   └── code.diff             uncommitted changes at submit time, if there were any
    └── _incomplete_task<N>/      only if a task was killed before run_case.py could save its output
```

`run_info.json`'s `outcome.status` tells you whether a run finished or diverged (went NaN, with `first_nan_time`). It also has `outcome.interrupted`, which is true when a run was stopped. If a task reaches its time limit, Slurm warns it 5 minutes early: the run is stopped, recorded as interrupted, and copied back.

To rebuild the index by hand, e.g. after moving run folders around:
```bash
cd ~/finray_fsi/fish-ib2d/cases/FinRay && python run_case.py index --all ~/finray_fsi/sim_data/fish-ib2d/FinRay/<sweep>/results
```

### 4.5 Test before a big sweep

Never start a large sweep with untested settings. First:
1. Submit a one-case sweep at the real `Tfinal` (e.g. the `base` sweep, or a one-case sweep added to `sweeps.py`) as e.g. `test_<something>`.
2. Check the log, `outcome.status` in `run_info.json`, and the time and memory used (`sacct`, 4.3). Then set the resources (4.2).
3. Submit the real sweep.

For a very quick check that the code runs at all, use the `smoke` sweep (0.05 s), or an interactive session (§5).

## 5. Interactive runs (testing, debugging)

```bash
srun -p Interactive --cpus-per-task=1 --mem=4G --time=00:30:00 --pty bash
source /opt/conda/etc/profile.d/conda.sh && conda activate ib2d
cd ~/finray_fsi/fish-ib2d/cases/FinRay
python run_case.py --sweep smoke --case 0   # one case; without --case: the whole sweep
exit                                 # give the node back as soon as you're done
```

Runs started this way write to `cases/FinRay/results_approach2/` in the repo, on the shared filesystem. That's fine for tests, but use `submit_sweep.sh` for real sweeps. For the geometry alone, run `MPLBACKEND=Agg python FinRay_Geom.py` (it writes `finray.*`, `mesh_info.json` and `geometry_preview.png`). The fluid-free stiffness calibration (`MPLBACKEND=Agg python verify_structure.py`, seconds to minutes) also belongs here, not on the head node.

## 6. Movies and analysis

- **Movies on the Mac:** VisIt isn't installed on the cluster, so pack each case into one file and copy it to your Mac, then run `make_movie.py` there as the README describes.
  ```bash
  # on the cluster
  cd ~/finray_fsi/sim_data/fish-ib2d/FinRay/<sweep>/results && tar czf <run_folder>.tgz <run_folder>
  # on the Mac
  scp icf:~/finray_fsi/sim_data/fish-ib2d/FinRay/<sweep>/results/<run_folder>.tgz .
  ```
- **Movies on the cluster (optional):** unpack VisIt's Linux build in your home directory (e.g. `~/opt/visit`; no admin rights needed) and add it to `PATH`. Then, on a compute node: `visit -quiet -cli -nowin -s make_movie.py -- <run_folder> out.mp4`. `ffmpeg` is installed on the compute nodes.

## 7. Git conventions in this repo

- **Branches:** work on `main`. `approach1_trimmed` is kept only for comparison. `approach2` and `wing-approach2` were merged into `main` and deleted.
- **Generated files are not tracked:** `finray.*`, `wing.*`, `geometry_preview.png` and `verify_structure.png` are listed in `.gitignore`. Regenerate them with the case's `*_Geom.py`.
- **Runs leave the tracked inputs unchanged.** `case_params.json` and `input2d` are rewritten during each case, then restored afterwards. So `git diff` on them only ever shows your own edits. Don't commit temporary test values.
- **The replication kit** (`replication/`) re-runs past simulations at the commits that produced them. All of those commits are in `main`'s history. Its `replicate_array.sbatch` is **not yet adapted to this cluster**. Before using it, add `#SBATCH --partition=ICF-Free`, the conda activation, and an `OUT=` outside the repo.

## 8. Pitfalls

- **Forgetting `-p ICF-Free`** sends jobs to `Teaching`. Batch jobs should always set the partition.
- **Single quotes keep `$USER` as literal text:** `'/disk/scratch/$USER'` creates a folder literally named `$USER`. Use double quotes.
- **Python imports are slow the first time** on the shared filesystem, e.g. 1–2 minutes for numpy, numba or vtk. That doesn't mean the job is stuck.
- **A leftover `viz_IB2d/` in `cases/FinRay/`**, from a local run that was killed, makes the next `run_case.py` refuse to start. Move that folder into the run's results folder, or delete it.
- **The wing case** has its own `cases/NACA0020_Wing/slurm/submit_sweep.sh` (same as FinRay's, plus `TIME_LIMIT` and `MEM`; data in `~/finray_fsi/sim_data/fish-ib2d/NACA0020_Wing/`). Its runs need about 0.8 GB; Nx 384 runs take ~4 h (set `TIME_LIMIT`). IB2d often *crashes* on divergence instead of writing NaN frames: check `analysis.stable` in `run_info.json`, not the exit code or `outcome.status` alone.

## Rules for coding agents

1. **Never run simulations, geometry builds or tests that use real CPU on the head node.** Run them with `sbatch`, or with `srun -p Interactive` for short checks. Light commands are fine on the head node: `git`, editing, `--list`, `submit_sweep.sh`, reading logs.
2. **Don't submit jobs yourself (the user's preference): hand the user an ordered plan** of what to submit (exact commands, what each job runs, which results to check). Ask before cancelling jobs, and before deleting or moving anything in `~/finray_fsi/sim_data`. Use `sbatch --test-only` to check a job without queuing it.
3. **Always set `--partition=ICF-Free`**, 1 CPU and no GPU for IB2d. Keep at most about 8 tasks running at once unless the user agrees otherwise.
4. **Don't write large output into the repo.** Real sweeps go through `submit_sweep.sh` to `~/finray_fsi/sim_data`.
5. **Check before overwriting:** sweep names must be new, and never delete run folders to make space without asking.
6. **Before committing,** check that `git status` doesn't include generated files or temporary test values in `case_params.json`, `input2d` or `sweeps.py`. Ask before pushing.
7. **Report outcomes from the files,** not from the exit code alone: `run_info.json` (`outcome.status`, `first_nan_time`), the job log, and `sacct`. A run that diverged still exits normally.
