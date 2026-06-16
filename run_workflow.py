import argparse
import subprocess
import json
import tempfile
import os
import sys

def main():
    parser = argparse.ArgumentParser(description="Wrapper for Snakemake workflow")

    # Workflow config
    parser.add_argument("--directories", nargs="+", required=True,
                        help="List of input directories, MUST be absolute paths")
    parser.add_argument("--dietmicrobe", action="store_true", help="Enable dietmicrobe based analysis")
    parser.add_argument("--dietmicrobehost", action="store_true", help="Enable dietmicrobehost based analysis")
    parser.add_argument("--metabolome", action="store_true", help="Include if comparison to known metabolome is wanted")
    parser.add_argument("--e-weights", action="store_true", help="Enable edge weights: read abundance")
    parser.add_argument("--n-weights", action="store_true", help="Enable node weights: food frequency")
    parser.add_argument("--include-orgs", action="store_true", help="Include organisms")
    parser.add_argument("--abundance-col", type=str, default="Abundance_RPKs",
                        help="Column name for abundance")
    parser.add_argument("--all-food", action="store_true", help="Enable use of all foods from foodb")

    # Snakemake execution options
    parser.add_argument("--cores", type=int, default=1, help="Number of cores to use")
    parser.add_argument("--profile", type=str, help="Snakemake profile to use")
    parser.add_argument("--dry-run", "-n", action="store_true", help="Perform a dry-run")
    parser.add_argument("--force", nargs="*", default=None,
                        help="Force rerun of rules, e.g., '--force all' or '--force rule1 rule2'")
    parser.add_argument("--rerun-incomplete", action="store_true",
                        help="Rerun rules that have incomplete output files")
    parser.add_argument("--unlock", action="store_true", help="Unlock Snakemake working directory")

    args = parser.parse_args()

    # Convert directories list to comma-separated string
    directories_str = ",".join(args.directories)

    # Prepare config dict
    config_args = {
        "directories": directories_str,
        "dietmicrobe": args.dietmicrobe,
        "dietmicrobehost": args.dietmicrobehost,
        "metabolome": args.metabolome,
        "e_weights": args.e_weights,
        "n_weights": args.n_weights,
        "include_orgs": args.include_orgs,
        "abundance_col": args.abundance_col,
        "all_food": args.all_food
    }

    # Write config to a temporary JSON file
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as tmp:
            json.dump(config_args, tmp, indent=2)
            tmp_path = tmp.name

        # Build base Snakemake command
        if args.unlock:
            cmd = f"snakemake --snakefile Snakefile --unlock --configfile {tmp_path}"
        else: 
            cmd = f"snakemake --snakefile Snakefile --cores {args.cores} --configfile {tmp_path}"

        # Optional flags
        if args.profile:
            cmd += f" --profile {args.profile}"
        if args.dry_run:
            cmd += " -n"
        if args.force is not None:
            if "all" in args.force:
                cmd += " --forceall"
            else:
                rules_to_force = " ".join(args.force)
                cmd += f" --forcerun {rules_to_force}"
        if args.rerun_incomplete:
            cmd += " --rerun-incomplete"

        # Show command
        print("🚀 Running Snakemake command:")
        print(cmd)
        print(f"🧾 Using config file: {tmp_path}")

        # Execute
        subprocess.run(cmd, shell=True, check=True)

    except subprocess.CalledProcessError as e:
        print(f"❌ Snakemake failed with return code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)

    finally:
        # Clean up temp config file
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)

if __name__ == "__main__":
    main()
