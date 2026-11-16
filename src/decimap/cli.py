"""Command line interface for DECIMAP."""
import argparse
import sys

from . import __version__
from .benchmark import run_benchmark
from .validate import run_validation


def _cmd_benchmark(args):
    frame = run_benchmark()
    if args.output:
        frame.to_csv(args.output, index=False, encoding="utf-8-sig")
        print(f"Wrote {args.output}")
    print(frame.to_string(index=False))
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

    p_bench = sub.add_parser("benchmark", help="run the symmetric de-circularized benchmark")
    p_bench.add_argument("--output", help="optional path to write the metrics CSV")
    p_bench.set_defaults(func=_cmd_benchmark)

    p_validate = sub.add_parser("validate", help="verify bundled data and benchmark values")
    p_validate.set_defaults(func=_cmd_validate)

    p_check = sub.add_parser("check", help="alias of validate")
    p_check.set_defaults(func=_cmd_validate)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
