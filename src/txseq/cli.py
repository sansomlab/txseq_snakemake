import subprocess
import shutil

import yaml

from argparse import ArgumentParser as AP
from importlib.resources import files
from pathlib import Path

from rich_argparse import RawDescriptionRichHelpFormatter

from txseq import __version__

RED = "\033[91m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
RESET = "\033[0m"

PLACEHOLDER_ACCOUNT = "slurm_account"

# Single source of truth for available workflows.
# The keys are the .smk filenames (without extension) in txseq/workflow/.
# name -> (one-line summary, [steps])
WORKFLOWS = {
    "ensembl": (
        "Post-process Ensembl annotations into reference files suitable for "
        "RNA-seq analysis (primary contigs only; PAR-masked).",
        [],
    ),
    "hisat_index": (
        "Build a HISAT2 genome index from the prepared reference.",
        [],
    ),
    "salmon_index": (
        "Build a Salmon transcriptome index from the prepared reference.",
        [],
    ),
    "fastqc": (
        "Run FastQC and post-process the output for downstream visualisation.",
        [],
    ),
    "preprocess_quantseq": (
        "Preprocess QuantSeq FASTQ files prior to mapping.",
        [],
    ),
    "hisat": (
        "Align reads and quantify gene expression from FASTQ files using HISAT2.",
        [],
    ),
    "salmon": (
        "Quantify gene expression from FASTQ files using Salmon and tximeta.",
        [],
    ),
    "feature_counts": (
        "Count reads mapped to transcript/gene models using featureCounts "
        "(Subread).",
        [],
    ),
    "bamqc": (
        "Compute BAM QC statistics using Picard tools and custom scripts.",
        [],
    ),
}


def workflows_help():
    """Render a compact workflow catalogue for the --help epilog."""
    width = max(len(name) for name in WORKFLOWS)
    lines = ["[bold]Available workflows:[/]", ""]
    for name, (summary, _steps) in WORKFLOWS.items():
        lines.append(f"  [bold cyan]{name:<{width}}[/]  {summary}")
    lines.append("")
    lines.append("Run a workflow with:  [bold]txseq <workflow> <task>[/]")
    lines.append("  <task> is one of: config | full | <rulename>")
    return "\n".join(lines)


def welcome_message():
    width = 68

    print()
    print("=" * width)
    print(f"Welcome to txseq (v{__version__})!".center(width))
    print("Snakemake workflows for bulk RNA-sequencing analysis.".center(width))
    print("=" * width)
    print()


def get_default_slurm_account():
    """Return the user's default Slurm account, or None if undeterminable."""
    if shutil.which("sacctmgr") is None:
        return None
    try:
        result = subprocess.run(
            ["sacctmgr", "-nP", "show", "user", "$USER", "format=DefaultAccount"],
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (subprocess.SubprocessError, OSError):
        return None
    if result.returncode != 0:
        return None
    account = result.stdout.strip()
    return account or None


def get_valid_slurm_accounts():
    """Return the accounts the user is associated with, or None if undeterminable."""
    if shutil.which("sacctmgr") is None:
        return None
    try:
        result = subprocess.run(
            ["sacctmgr", "-nP", "show", "assoc", "user=$USER", "format=Account"],
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (subprocess.SubprocessError, OSError):
        return None
    if result.returncode != 0:
        return None
    accounts = {line.strip() for line in result.stdout.splitlines() if line.strip()}
    return accounts or None


def _account_error(message):
    print(f"{RED}[ERROR] {message}{RESET}")
    default_account = get_default_slurm_account()
    if default_account:
        print(
            f"{YELLOW}[HINT] Your default Slurm account appears to be "
            f"'{default_account}'. Set 'resources.account' in the workflow "
            f"config accordingly.{RESET}"
        )
    else:
        print(
            f"{YELLOW}[HINT] To find your default Slurm account, run:\n"
            f"           sacctmgr -nP show user $USER format=DefaultAccount{RESET}"
        )
    raise SystemExit(1)


def validate_config(workflow):
    """Pre-flight check of the workflow config before Snakemake is invoked.

    Validates that resources.account is set to a real Slurm account rather than
    the template placeholder, so failures surface here rather than inside a
    verbose Snakemake/Slurm traceback.
    """
    config_path = Path.cwd().joinpath(f"{workflow}.yaml")
    if not config_path.exists():
        print(
            f"{RED}[ERROR] Config file '{config_path.name}' not found in the "
            f"current directory. Run 'txseq {workflow} config' to create it.{RESET}"
        )
        raise SystemExit(1)

    with config_path.open() as fh:
        config = yaml.safe_load(fh) or {}

    resources = config.get("resources") or {}
    account = resources.get("account")

    if not account:
        _account_error(
            f"No Slurm account set under 'resources.account' in "
            f"{config_path.name}."
        )

    if account == PLACEHOLDER_ACCOUNT:
        _account_error(
            f"'resources.account' in {config_path.name} is still the template "
            f"placeholder '{PLACEHOLDER_ACCOUNT}'. Assign a valid Slurm account."
        )

    valid_accounts = get_valid_slurm_accounts()
    if valid_accounts is not None and account not in valid_accounts:
        _account_error(
            f"Slurm account '{account}' in {config_path.name} is not one of "
            f"your associated accounts ({', '.join(sorted(valid_accounts))})."
        )

    print(f"{GREEN}[INFO] Config validated: using Slurm account '{account}'.{RESET}")


def copy_config_file(workflow):
    yaml_dir = files("txseq").joinpath("configfiles")
    # config files are named config_<workflow>.yaml in the package
    src = yaml_dir.joinpath(f"config_{workflow}.yaml")
    dst = Path.cwd().joinpath(f"{workflow}.yaml")
    if not dst.exists():
        dst.write_text(src.read_text())
        print(f"[INFO] Copied {workflow} config file to the current directory.")
    else:
        print(f"{RED}[ERROR] {workflow} config file already exists.{RESET}")
        raise FileExistsError()


def run_snakemake(smkpath, configpath, task, cores, jobs, dry_run, lock):
    cmd = [
        "snakemake",
        "--snakefile",
        str(smkpath),
        "--configfile",
        str(configpath),
        "--cores",
        cores,
        "--jobs",
        jobs,
        "-p",
    ]
    cmd.extend(["--config", f"lock={str(lock).lower()}"])
    # Slurm executor: per-rule resources (mem_mb, runtime, slurm_partition,
    # slurm_account) are read from each rule's `resources:` block, which the
    # workflows populate from the `resources:` section of the config yaml.
    cmd.extend(["--executor", "slurm"])
    if dry_run:
        cmd.append("--dry-run")
    cmd.extend([task])
    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    welcome_message()

    vinfo = f"%(prog)s {__version__}"

    p = AP(
        description="Snakemake workflows for bulk RNA-sequencing analysis.",
        epilog=workflows_help(),
        formatter_class=RawDescriptionRichHelpFormatter,
    )
    p.add_argument(
        "workflow",
        choices=WORKFLOWS.keys(),
        metavar="workflow",
        help="The workflow to run (see the list below).",
    )
    p.add_argument("task", help="The task to run, [config|full|<rulename>].")
    p.add_argument("--dry", action="store_true", help="Perform a dry run only.")
    p.add_argument("--cores", default="all", help="Number of cores to use.")
    p.add_argument("--jobs", default="1", help="Number of parallel jobs to run.")
    p.add_argument("--lock", action="store_true", help="Lock output by chmod.")
    p.add_argument("-v", "--version", action="version", version=vinfo)
    args = p.parse_args()

    smkpath = files("txseq").joinpath("workflow", f"{args.workflow}.smk")

    if args.lock and args.task == "config":
        print(f"{GREEN}[INFO] `--lock` is ignored when running config task.{RESET}")

    if args.task == "config":
        copy_config_file(args.workflow)
    else:
        validate_config(args.workflow)
        configpath = Path.cwd().joinpath(f"{args.workflow}.yaml")
        run_snakemake(
            smkpath,
            configpath,
            args.task,
            str(args.cores),
            str(args.jobs),
            args.dry,
            str(args.lock).lower(),
        )

    print(f"{GREEN}Completed successfully!{RESET}")
