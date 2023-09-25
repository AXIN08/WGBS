import pandas as pd
import os
import argparse
import sys
import logging
import yaml
from pathlib import Path
from typing import Any,Dict


def parser() -> argparse.ArgumentParser:
    main_parser = argparse.ArgumentParser()

    main_parser.add_argument(
        "-f",
        "--samplefile",
        help="Path to the sample with sample name and raw data path",
        type=str,
        metavar="PATH",
        required=True
    )

    main_parser.add_argument(
        "-o",
        "--outdir",
        help="Path to output directory",
        type=str,
        metavar="PATH",
        default="."
    )

    main_parser.add_argument(
        "-C",
        "--config_output",
        help="Path to the output config.yaml file",
        type=str,
        metavar="PATH"
    )
    main_parser.add_argument(
        "-r",
        "--ref_dir",
        help="Path to the reference directory",
        type=str,
        metavar="PATH",
        default="/home/liuxin/data_pasteur/12_epigenome/CHM13_ref"
    )
    return main_parser

def parse_args(args: Any) -> argparse.ArgumentParser:
    return parser().parse_args(args)

def args_to_dict(args: argparse.ArgumentParser) -> Dict[str, Any]:
    result_dict = {
            #"config": os.path.abspath(args.config_output),
            "outdir": os.path.abspath(args.outdir),
            "threads": {
                "soapnuke": 8,
                "bis_align": 6,
                "bis_dedul": 8,
                "bis_call": 8},
            "ref": args.ref_dir,
            "sample":os.path.abspath(args.samplefile)
            }
    logging.debug(result_dict)
    return result_dict

def dict_to_yaml(indict: Dict[str, Any]) -> str:
    return yaml.dump(indict, default_flow_style=False)

def main(args: argparse.ArgumentParser) -> None:
    logging.debug("Building configuration file:")
    config_params = args_to_dict(args)
    output_path = Path(args.config_output) / "config.yaml"

    with output_path.open("w") as config_yaml:
        logging.debug(f"Saving results to {str(output_path)}")
        config_yaml.write(dict_to_yaml(config_params))


if __name__ == "__main__":
    # Parsing command line
    args = parse_args(sys.argv[1:])
    #makedirs("logs/prepare")

    try:
        logging.debug("Preparing configuration")
        main(args)
    except Exception as e:
        logging.exception("%s", e)
        raise
