import argparse
import sys

from typing import NoReturn

def main() -> NoReturn:

    args = parse_args(sys.argv[1:])

    print(args.year)

def parse_args(args: list):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--year",
        required=True,
        help="Year to run model for",
    )
    args = parser.parse_args(args)
    return args

if __name__ == "__main__":
    main()