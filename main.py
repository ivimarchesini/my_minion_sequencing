"""CLI wrapper for my_minion_sequencing.

Supports `--select-run` to launch the Textual TUI and optionally run the
remote scp (executed via `ssh minit@... "scp -r ..."` so prompts are
forwarded to your terminal). Also supports `--list` to print available runs
non-interactively.
"""

from __future__ import annotations

import argparse
import getpass
import sys
import subprocess
import asyncio


def main() -> int:
    parser = argparse.ArgumentParser(description="my_minion_sequencing helper CLI")
    parser.add_argument("--select-run", action="store_true", help="Launch TUI to select a pod5_pass folder and optionally run scp")
    parser.add_argument("--list", action="store_true", help="List pod5_pass folders (non-interactive)")
    args = parser.parse_args()

    if args.select_run:
        try:
            import select_run
            from select_run import RunSelectorApp
        except Exception as e:
            print("Failed to import select_run. Install dependencies: pip install textual paramiko", file=sys.stderr)
            print(e, file=sys.stderr)
            return 2

        app = RunSelectorApp()
        app.run()

        sel = getattr(app, "selected_path", None)
        scp_cmd = getattr(app, "suggested_scp", None)

        if not sel:
            print("No selection made.")
            return 0

        # Ask whether to run scp now
        reply = input("Run scp now (will SSH to minit and execute scp on that host)? [y/N]: ").strip().lower()
        if reply != "y":
            print("Skipping scp. You can run the suggested command shown earlier.")
            return 0

        # Execute scp locally so your terminal will prompt for both hosts.
        # Build args to avoid shell quoting issues; use -3 so scp connects to
        # both remotes from this machine.
        try:
            # Recompute run_name and paths to be safe
            sel_path = sel
            seg, _ = select_run.parse_run_info(sel_path)
            run_name = seg or "run"
            src = sel_path.rstrip("/") + "/*"
            dst = select_run.DEFAULT_REMOTE_TARGET.rstrip("/") + f"/{run_name}/"

            scp_args = ["scp", "-3", "-r", f"{select_run.SSH_USER}@{select_run.SSH_HOST}:{src}", dst]

            print("Running:", " ".join(scp_args))
            proc = subprocess.run(scp_args)
            if proc.returncode != 0:
                return proc.returncode

            # Ask whether to run the processing on the Ramses cluster now
            reply2 = input("Run processing on Ramses now (ssh + srun)? [y/N]: ").strip().lower()
            if reply2 != "y":
                print("Skipping remote run. You can ssh to the cluster manually to start processing.")
                return 0

            # Determine ssh target from the DEFAULT_REMOTE_TARGET in select_run
            # DEFAULT_REMOTE_TARGET looks like: user@host:/path/to/dir
            try:
                target = select_run.DEFAULT_REMOTE_TARGET.split(":", 1)[0]
            except Exception:
                target = "imarches@ramses4.itcc.uni-koeln.de"

            # Build the run name again to be safe
            sel_path = sel
            seg, _ = select_run.parse_run_info(sel_path)
            run_name = seg or "run"

            # Build the remote srun command. $(date +%F) will be evaluated on the remote host.
            remote_cmd = (
                f"srun -A virology -p interactive --gpus=4 --time=05:00:00 --mem=24gb "
                f"-J nanopore.{run_name}.$(date +%F) --pty bash -i -c \"mi_nanopore_basecalling {run_name}\"")

            print("Running remote command on", target)
            print(remote_cmd)
            # Call ssh and let it handle interactive prompts (password or key)
            rc = subprocess.run(["ssh", target, remote_cmd]).returncode
            return rc
        except KeyboardInterrupt:
            print("Interrupted by user.")
            return 130

    if args.list:
        try:
            from select_run import fetch_remote_paths
        except Exception as e:
            print("Failed to import select_run. Install dependencies: pip install paramiko", file=sys.stderr)
            print(e, file=sys.stderr)
            return 2

        pw = getpass.getpass("Password for minit@10.212.1.44 (leave empty to try key auth): ")
        try:
            lines = asyncio.run(asyncio.to_thread(fetch_remote_paths, pw or None))
        except Exception as e:
            print(f"Error fetching remote paths: {e}", file=sys.stderr)
            return 3

        if not lines:
            print("<no results>")
            return 0

        for line in lines:
            print(line)

        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
