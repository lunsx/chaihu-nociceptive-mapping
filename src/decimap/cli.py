"""Command line interface for DECIMAP."""
import argparse
import sys

from . import __version__
from .benchmark import run_benchmark
from .degree_aware import (
    run_calibration,
    run_degree_analysis,
    run_portability,
)
from .simulation import run_simulation
from .validate import run_validation


def _print_frame(frame, output):
    if output:
        frame.to_csv(output, index=False, encoding="utf-8-sig")
        print(f"Wrote {output}")
    print(frame.to_string(index=False))


def _cmd_benchmark(args):
    _print_frame(run_benchmark(), args.output)
    return 0


def _cmd_degree(args):
    _print_frame(run_degree_analysis(), args.output)
    return 0


def _cmd_calibrate(args):
    _print_frame(run_calibration(), args.output)
    return 0


def _cmd_portability(args):
    _print_frame(run_portability(), args.output)
    return 0


def _cmd_simulate(args):
    frame = run_simulation(n_reps=args.reps, n_null=args.null)
    _print_frame(frame, args.output)
    return 0


def _cmd_validate(args):
    ok, messages = run_validation()
    for message in messages:
        print(message)
    print("VALIDATION PASSED" if ok else "VALIDATION FAILED")
    return 0 if ok else 1


def build_parser():
    parser = argparse.ArgumentParser(
        prog="decimap",
        description="De-circular, evidence-graded computational mapping tools.",
    )
    parser.add_argument("--version", action="version", version=f"DECIMAP {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    specs = [
        ("benchmark", "run the symmetric de-circularized benchmark", _cmd_benchmark, False),
        ("degree", "run the degree-aware de-circularized enrichment", _cmd_degree, False),
        ("calibrate", "run negative-control empirical calibration", _cmd_calibrate, False),
        ("portability", "run the conditional cross-formula portability test", _cmd_portability, False),
        ("simulate", "run the synthetic-graph simulation benchmark", _cmd_simulate, True),
        ("validate", "verify bundled data and benchmark values", _cmd_validate, False),
    ]
    for name, help_text, func, is_sim in specs:
        p = sub.add_parser(name, help=help_text)
        if func is not _cmd_validate:
            p.add_argument("--output", help="optional path to write the result CSV")
        if is_sim:
            p.add_argument("--reps", type=int, default=500, help="simulation replicates")
            p.add_argument("--null", type=int, default=1000, help="per-test null draws")
        p.set_defaults(func=func)

    p_check = sub.add_parser("check", help="alias of validate")
    p_check.set_defaults(func=_cmd_validate)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
