import getpass
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
DEFAULT_PARTITION = "short"

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


def _sacctmgr(*args):
    """Run sacctmgr and return its stdout lines, or None if unavailable/failed.

    Arguments are passed as a list (no shell), so the user name is resolved in
    Python rather than relying on shell expansion of $USER.
    """
    if shutil.which("sacctmgr") is None:
        return None
    try:
        result = subprocess.run(
            ["sacctmgr", "-nP", *args],
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (subprocess.SubprocessError, OSError):
        return None
    if result.returncode != 0:
        return None
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def get_default_slurm_account():
    """Return the user's default Slurm account, or None if undeterminable."""
    lines = _sacctmgr("show", "user", getpass.getuser(), "format=DefaultAccount")
    return lines[0] if lines else None


def get_valid_slurm_accounts():
    """Return the accounts the user is associated with, or None if undeterminable."""
    lines = _sacctmgr("show", "assoc", f"user={getpass.getuser()}", "format=Account")
    return {line.lower() for line in lines} if lines else None


def _account_error(message):
    print(f"{RED}[ERROR] {message}{RESET}")
    default_account = get_default_slurm_account()
    if default_account:
        print(
            f"{YELLOW}[HINT] Your default Slurm account appears to be "
            f"'{default_account}'. Set 'slurm.account' in the workflow "
            f"config accordingly.{RESET}"
        )
    else:
        print(
            f"{YELLOW}[HINT] To list the Slurm accounts available to you, run:\n"
            f"           sacctmgr -nP show assoc user=$USER format=Account{RESET}"
        )
    raise SystemExit(1)


def _slurm_settings(config, config_name):
    """Return (account, partition) from the config.

    The canonical location is the top-level `slurm:` section. Earlier config
    templates placed these under `resources:` as `account`/`partition` or
    `slurm_account`/`slurm_partition`; these are still honoured, with a
    warning, so that existing run directories continue to work.
    """
    slurm = config.get("slurm") or {}
    resources = config.get("resources") or {}

    account = slurm.get("account")
    partition = slurm.get("partition")

    legacy = []
    if account is None:
        for key in ("slurm_account", "account"):
            if resources.get(key) is not None:
                account = resources[key]
                legacy.append(f"resources.{key}")
                break
    if partition is None:
        for key in ("slurm_partition", "partition"):
            if resources.get(key) is not None:
                partition = resources[key]
                legacy.append(f"resources.{key}")
                break

    if legacy:
        print(
            f"{YELLOW}[WARNING] {config_name} uses the deprecated key(s) "
            f"{', '.join(legacy)}. Move the Slurm account and partition to a "
            f"top-level 'slurm:' section (see 'txseq <workflow> config').{RESET}"
        )

    account = str(account).strip() if account is not None else None
    partition = str(partition).strip() if partition is not None else DEFAULT_PARTITION
    return account or None, partition or DEFAULT_PARTITION


def validate_config(workflow):
    """Pre-flight check of the workflow config before Snakemake is invoked.

    Validates that a Slurm account is set to a real account rather than the
    template placeholder, so failures surface here rather than inside a
    verbose Snakemake/Slurm traceback. Returns (account, partition).
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

    account, partition = _slurm_settings(config, config_path.name)

    if not account:
        _account_error(
            f"No Slurm account set under 'slurm.account' in {config_path.name}."
        )

    if account == PLACEHOLDER_ACCOUNT:
        _account_error(
            f"'slurm.account' in {config_path.name} is still the template "
            f"placeholder '{PLACEHOLDER_ACCOUNT}'. Assign a valid Slurm account."
        )

    valid_accounts = get_valid_slurm_accounts()
    if valid_accounts is not None and account.lower() not in valid_accounts:
        _account_error(
            f"Slurm account '{account}' in {config_path.name} is not one of "
            f"your associated accounts ({', '.join(sorted(valid_accounts))})."
        )

    print(
        f"{GREEN}[INFO] Config validated: Slurm account '{account}', "
        f"partition '{partition}'.{RESET}"
    )
    return account, partition


def copy_config_file(workflow):
    yaml_dir = files("txseq").joinpath("configfiles")
    # config files are named config_<workflow>.yaml in the package
    src = yaml_dir.joinpath(f"config_{workflow}.yaml")
    dst = Path.cwd().joinpath(f"{workflow}.yaml")
    if not dst.exists():
        text = src.read_text()
        # Pre-fill the user's default Slurm account where it can be determined.
        default_account = get_default_slurm_account()
        if default_account:
            text = text.replace(
                f"account: {PLACEHOLDER_ACCOUNT}", f"account: {default_account}", 1
            )
        dst.write_text(text)
        print(f"[INFO] Copied {workflow} config file to the current directory.")
        if default_account:
            print(
                f"[INFO] Set 'slurm.account' to your default Slurm account "
                f"'{default_account}'; edit {dst.name} to use a different one."
            )
    else:
        print(f"{RED}[ERROR] {workflow} config file already exists.{RESET}")
        raise FileExistsError()


def run_snakemake(
    smkpath, configpath, task, cores, jobs, dry_run, lock, account, partition
):
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
    # Slurm executor. Per-rule resources (mem_mb, runtime) come from each
    # rule's `resources:` block. The account and partition are cluster-level
    # settings, so they are supplied as default resources: this applies them
    # to every submitted job, including rules that define their own
    # resources, without each rule having to repeat them.
    cmd.extend(
        [
            "--default-resources",
            # Values are quoted because Snakemake evaluates default resources
            # as Python expressions; an unquoted account name could otherwise
            # resolve to a builtin or a number.
            f"slurm_account={account!r}",
            f"slurm_partition={partition!r}",
        ]
    )
    # --default-resources and --config accept several values; a following
    # option is needed to close them before the positional target.
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
        account, partition = validate_config(args.workflow)
        configpath = Path.cwd().joinpath(f"{args.workflow}.yaml")
        run_snakemake(
            smkpath,
            configpath,
            args.task,
            str(args.cores),
            str(args.jobs),
            args.dry,
            str(args.lock).lower(),
            account,
            partition,
        )

    print(f"{GREEN}Completed successfully!{RESET}")
