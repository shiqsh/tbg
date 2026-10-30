#!/usr/bin/env python3
"""Run TBG targets after invalidating fpm 0.13's cached executables.

fpm can reuse an older executable after its objects were compiled by another
invocation. Preserve the previous executable, then let fpm link the target from
its current objects. Scientific sources, object files and dependencies stay put.
"""

import fnmatch
import os
from pathlib import Path
import re
import sys


def prepare_executables(project, args):
    """Only move fpm-generated executables selected by this invocation."""
    options = args[:args.index("--")] if "--" in args else args
    if not options or any(x in options for x in ("--help", "--version", "--list")):
        return []

    patterns = ["*"] if "--all" in options else options[0].split(",")
    build_name = os.environ.get("FPM_BUILD_DIR") or "build"
    for index, option in enumerate(options):
        if option == "--build-dir" and index + 1 < len(options):
            build_name = options[index + 1]
        elif option.startswith("--build-dir="):
            build_name = option.split("=", 1)[1]
    build_dir = Path(build_name)
    if not build_dir.is_absolute():
        build_dir = project / build_dir
    if not build_dir.is_dir():
        return []

    target_kind = "example" if "--example" in options else "app"
    moved = []
    for config in sorted(build_dir.iterdir()):
        # Limit changes to fpm's compiler/hash output directories.
        if not config.is_dir() or not re.fullmatch(r".+_[0-9A-Fa-f]{16}", config.name):
            continue
        outputs = config / target_kind
        if not outputs.is_dir():
            continue
        for executable in sorted(outputs.iterdir()):
            if executable.name.endswith((".frt-previous", ".log")):
                continue
            if executable.is_symlink() or not executable.is_file():
                continue
            if not os.access(executable, os.X_OK):
                continue
            if not any(fnmatch.fnmatchcase(executable.name, p) for p in patterns):
                continue
            previous = executable.with_name(executable.name + ".frt-previous")
            executable.replace(previous)
            moved.append(previous)
    return moved


def main():
    project = Path(__file__).resolve().parents[1] / "fortran"
    args = sys.argv[1:]
    os.chdir(project)
    try:
        moved = prepare_executables(project, args)
    except OSError as error:
        print(f"[frt] 无法准备重新链接，已停止运行：{error}", file=sys.stderr)
        return 1
    if moved:
        print("[frt] 已保留旧执行文件为 .frt-previous；本次将重新链接所选程序。", flush=True)
    # exec preserves fpm's exit status and the terminal's Ctrl-C handling.
    os.execvp("fpm", ["fpm", "run", "--runner", "/usr/bin/time -l", "--target", *args])


if __name__ == "__main__":
    sys.exit(main())
