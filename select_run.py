#!/usr/bin/env python3
"""
TUI for selecting a remote `pod5_pass` folder via SSH.

Usage:
  - Install dependencies: `pip install textual paramiko`
  - Run: `python select_run.py`

The app will ask for a password (leave blank to try key-based auth).
It connects to `minit@10.212.1.44`, runs the find command, lists results
and lets the user select one with arrow keys and Enter. The selected
path is printed to stdout and the app exits.
"""
from __future__ import annotations

import asyncio
import sys
from typing import List

import os
import re
import datetime
import paramiko
from textual.app import App, ComposeResult
from textual.containers import Container
from textual.widgets import Header, Footer, Input, Button, Static, ListView, ListItem
from textual import events


SSH_HOST = "10.212.1.44"
SSH_USER = "minit"
FIND_CMD = "find /data/run_* -type d -name 'pod5_pass' | sort -u"
# Default remote target used in README examples. Keep overridable by copying the
# suggested command and editing if necessary.
DEFAULT_REMOTE_TARGET = "imarches@ramses4.itcc.uni-koeln.de:/projects/virology/nanopore/input"


def fetch_remote_paths(password: str | None) -> List[str]:
    """Blocking SSH call using paramiko to run FIND_CMD and return lines."""
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    # If password is empty string, use None -> try key-based auth
    pw = password if password else None

    client.connect(hostname=SSH_HOST, username=SSH_USER, password=pw, look_for_keys=True, allow_agent=True, timeout=15)
    stdin, stdout, stderr = client.exec_command(FIND_CMD)
    out = stdout.read()
    err = stderr.read()
    client.close()

    if err:
        # decode error for raising
        raise RuntimeError(err.decode(errors="replace"))

    text = out.decode(errors="replace")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return lines


def extract_run_segment(path: str) -> str:
    """Return the run_XXXX segment from a path, or a fallback label.

    Examples:
      /data/run_20250723/no_sample_id/... -> run_20250723
    """
    try:
        parts = path.strip("/").split("/")
        for part in parts:
            if part.startswith("run_"):
                return part
        # fallback: use the second component if path starts with data/
        if len(parts) >= 2:
            return parts[1]
    except Exception:
        pass
    return os.path.basename(path) or path


def parse_run_info(path: str) -> tuple[str, datetime.date | None]:
    """Return (run_segment, parsed_date) or (label, None).

    parsed_date is a datetime.date when we find a segment like run_YYYYMMDD.
    """
    seg = extract_run_segment(path)
    m = re.match(r"run_(\d{8})$", seg)
    if m:
        ymd = m.group(1)
        try:
            dt = datetime.datetime.strptime(ymd, "%Y%m%d").date()
            return seg, dt
        except Exception:
            return seg, None
    return seg, None


class Message(Static):
    pass


class RunListItem(ListItem):
    """List item that stores both display label and full path."""

    def __init__(self, label: str, path: str) -> None:
        super().__init__(Static(label))
        self.path = path
        self.label = label


class RunSelectorApp(App):
    CSS_PATH = None
    BINDINGS = [("q","quit","Quit")]
    # populated when user selects an entry
    selected_path: str | None = None
    suggested_scp: str | None = None

    def compose(self) -> ComposeResult:
        yield Header(show_clock=False)
        with Container():
            yield Message("Enter password for minit@10.212.1.44 (leave blank to use SSH key):")
            yield Input(password=True, placeholder="Password (typing hidden)", id="pw")
            yield Button(label="Connect", id="connect")
            yield Message("")
            yield ListView(id="results")
        yield Footer()

    async def on_button_pressed(self, event: Button.Pressed) -> None:  # type: ignore[override]
        if event.button.id == "connect":
            inp = self.query_one(Input)
            password = inp.value or ""
            msg = self.query(Message)
            msg[1].update("Connecting... (will try key auth if password blank)")

            try:
                # Run blocking SSH call in thread
                lines = await asyncio.to_thread(fetch_remote_paths, password)
            except Exception as e:
                msg[1].update(f"Error: {e}")
                return

            results = self.query_one("#results", ListView)
            results.clear()
            if not lines:
                results.append(ListItem(Static("<no results>")))
                msg[1].update("No matching folders found.")
                return

            # Build entries with parsed dates, sort by date desc (newest first)
            entries: list[tuple[str, str, datetime.date | None]] = []
            for line in lines:
                seg, dt = parse_run_info(line)
                entries.append((line, seg, dt))

            # Use a very old date for None so they sort to the end
            oldest = datetime.date(1, 1, 1)
            entries.sort(key=lambda e: e[2] or oldest, reverse=True)

            today = datetime.date.today()
            for path, seg, dt in entries:
                # build human readable delta
                age_label = ""
                if dt is None:
                    age_label = "(unknown date)"
                else:
                    days = (today - dt).days
                    if days <= 0:
                        age_label = "today"
                    elif days < 30:
                        age_label = f"{days} days ago"
                    elif days < 365:
                        months = days // 30
                        age_label = f"{months} months ago"
                    else:
                        years = days // 365
                        age_label = f"{years} years ago"

                display = f"{seg} — {age_label}"
                results.append(RunListItem(display, path))

            msg[1].update("Use arrow keys to pick a folder and press Enter to accept.")

    async def on_key(self, event: events.Key) -> None:  # catch Enter from user selection
        if event.key == "enter":
            results = self.query_one("#results", ListView)
            # Try to get the highlighted index; fallback to first child
            idx = getattr(results, "index", None)
            if idx is None:
                try:
                    idx = 0 if len(results.children) > 0 else None
                except Exception:
                    idx = None

            if idx is None:
                return

            try:
                item = results.children[idx]
                # Prefer the stored full path if available
                selection_text = getattr(item, "path", None)
                if selection_text is None:
                    # fallback to the displayed text
                    text_widget = item.query_one(Static)
                    selection = text_widget.renderable if hasattr(text_widget, "renderable") else text_widget
                    selection_text = str(selection)
            except Exception:
                selection_text = ""

            if selection_text:
                # Store selection and suggested command on the app instance
                self.selected_path = selection_text

                seg, _ = parse_run_info(selection_text)
                run_name = seg or "run"
                src = selection_text.rstrip("/") + "/*"
                dst = DEFAULT_REMOTE_TARGET.rstrip("/") + f"/{run_name}/"

                # Suggest using rsync to copy only new/updated files and keep
                # folders in sync without retransmitting existing files.
                # Local rsync (pull from minit via ssh and write to ramses path):
                scp_local = f"rsync -av --update --progress -e ssh {SSH_USER}@{SSH_HOST}:{src} {dst}"
                # Alternate: run rsync on the minit host to push to ramses
                scp_remote = f"ssh {SSH_USER}@{SSH_HOST} \"rsync -av --update --progress {src} {dst}\""
                self.suggested_scp = scp_local

                # Print for users running `select_run.py` directly
                print(selection_text)
                print()
                print("Suggested local scp command (runs on your machine and will prompt for both hosts):")
                print(scp_local)
                print()
                print("Alternate (run scp on the minit host):")
                print(scp_remote)
                print()

                await self.action_quit()


if __name__ == "__main__":
    # Basic check for dependencies
    try:
        import textual  # noqa: F401
        import paramiko  # noqa: F401
    except Exception as e:
        print("Missing dependency. Install with: pip install textual paramiko")
        print(e)
        sys.exit(2)

    RunSelectorApp().run()
