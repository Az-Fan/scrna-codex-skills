# Compatibility

The client computer needs Codex and/or WispScience only. Analysis dependencies are resolved on the registered server where pixi runs.

Tested server baseline: Linux x86_64, Python 3.8.10, R 4.4.1, Seurat 5.x, and the project pixi locks. No GPU is required.

This skill uses the `03-integration` server pixi environment because the optional Harmony branch requires `harmony`. The `clustree` package is optional: when it is missing or incompatible with the installed `ggplot2`, the resolution tree is skipped, a `*_clustree_error.txt` note is written, and the rest of the run continues.

Run `python3 scripts/check_dependencies.py --pixi-root ~/projects/scrna_envs` on the server before execution. The checker and executor resolve the same exact pixi `Rscript`; system R fallback is disabled. The checker is read-only. Missing optional packages disable only the corresponding branch.
