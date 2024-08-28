#!/usr/bin/env python3

"""
Example unbiased simulation with pysages and lammps.

For a list of possible options for running the script pass `-h` as argument from the
command line, or call `get_args(["-h"])` if the module was loaded interactively.
"""

import argparse
import sys
import numpy as np
from lammps import lammps

import pysages
from pysages.backends import SamplingContext
from pysages.colvars import DihedralAngle
from pysages.methods import Unbiased, HistogramLogger


def generate_context(args="", script="lj.lmp", store_freq=1):
    """
    Returns a lammps simulation defined by the contents of `script` using `args` as
    initialization arguments.
    """
    context = lammps(cmdargs=args.split())
    context.file(script)
    # Allow for the retrieval of the unwrapped positions
    context.command(f"fix unwrap all store/state {store_freq} xu yu zu")
    return context


def get_args(argv):
    """Process the command-line arguments to this script."""

    available_args = [
        ("time-steps", "t", int, 2e6, "Number of simulation steps"),
        ("kokkos", "k", bool, True, "Whether to use Kokkos acceleration"),
    ]
    parser = argparse.ArgumentParser(description="Example script to run pysages with lammps")

    for name, short, T, val, doc in available_args:
        if T is bool:
            action = "store_" + str(val).lower()
            parser.add_argument("--" + name, "-" + short, action=action, help=doc)
        else:
            convert = (lambda x: int(float(x))) if T is int else T
            parser.add_argument("--" + name, "-" + short, type=convert, default=T(val), help=doc)

    return parser.parse_args(argv)


def main(argv):
    """Example simulation with pysages and lammps."""
    args = get_args(argv)

    context_args = {"store_freq": args.time_steps}
    if args.kokkos:
        # Passed to the lammps constructor as `cmdargs` when running the script
        # with the --kokkos (or -k) option
        context_args["args"] = "-k on g 1 -sf kk -pk kokkos newton on neigh half"

    # Setting the collective variable, method, and running the simulation
    cvs = [
        DihedralAngle([4,6,8,14]),DihedralAngle([6,8,14,16])
    ]
    method = Unbiased(cvs)
    callback = HistogramLogger(2000)
    sampling_context = SamplingContext(method, generate_context, callback=callback, context_args=context_args)
    result = pysages.run(sampling_context, args.time_steps)
    bins = 50
    histo, xedges, yedges = np.histogram2d(
        callback.data[:,0],
        callback.data[:,1],
        bins=bins,
        range=[[-np.pi,np.pi],[-np.pi,np.pi]]
        )
    xedges = (xedges[1:] + xedges[:-1]) / 2
    yedges = (yedges[1:] + yedges[:-1]) / 2
    mesh = np.dstack(np.meshgrid(xedges, yedges)).reshape(-1, 2)
    np.savetxt("hist-fromlogger.csv", np.hstack([mesh, histo.reshape(-1, 1)]))

if __name__ == "__main__":
    main(sys.argv[1:])
