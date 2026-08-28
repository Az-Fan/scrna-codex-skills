# Long-running execution with tmux

Use tmux when an execution is likely to exceed 10 minutes, has uncertain duration, or could be interrupted by a disconnected remote session. Keep dependency checks and the skill dry run in the foreground. Launch only the already-reviewed `--execute` command in tmux.

Use the bundled supervisor instead of assembling shell redirection manually:

```bash
python3 scripts/run_in_tmux.py \
  --session <project>-<skill-number>-<stage> \
  --cwd <project-root> \
  --log <stage-output>/tmux.log \
  --status <stage-output>/tmux_status.json \
  -- \
  python3 <installed-skill>/scripts/run.py --config <config> --execute
```

The command above is a dry run of the supervisor and only prints the resolved plan. Review it, then add `--execute` before the command separator to create the detached session. Session names must contain only letters, digits, dots, underscores, and hyphens, with at most 64 characters. Do not put secrets in session names or command arguments.

The supervisor refuses an active duplicate session. If it reports a duplicate, inspect that session and its files instead of starting another copy:

```bash
tmux ls
tmux capture-pane -pt <session>
tail -f <stage-output>/tmux.log
cat <stage-output>/tmux_status.json
tmux attach -t <session>
```

`tmux.log` is the supervisor stream and is separate from a skill's own `run.log`. `tmux_status.json` moves from `running` to `completed` or `failed` and records the exit code. Treat the skill run manifest and required artifacts as the final completion contract; a disappearing tmux session alone does not prove success.

Do not terminate a running session unless the user explicitly requests cancellation or continuing would be unsafe. On Wisp-managed projects, use the registered Wisp execution context to run these commands; do not create an independent SSH connection.
