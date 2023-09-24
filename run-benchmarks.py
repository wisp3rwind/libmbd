#! /usr/bin/env python3

import os
import subprocess
import time

size_and_reps = {
    "plain": [
        (1, 4000),
        (2, 4000),
        (3, 4000),
        (4, 4000),
        (5, 4000),
        (6, 500),
        # (7, 500),
        # (8, 500),
        # (9, 200),
        # (10, 200),
        # (15, 20),
        # (20, 20),
    ],
    "rsscs": [
        (1, 5000),
        (2, 5000),
        (3, 5000),
        (4, 2000),
        (5, 1000),
        (6, 500),
        # (7, 200),
        # (8, 150),
        # (9, 80),
        # (10, 50),
        # (15, 5),
        # (20, 2),
    ],
}

# uses localtime
timestamp = time.strftime("%Y-%m-%dT%H-%M-%S")

try:
    os.mkdir("benchmark-results")
except FileExistsError:
    pass

commit = subprocess.check_output(["git", "describe", "--long", "--dirty"])
commit = commit.decode("utf-8").strip("\n")
dirname = f"benchmark-results/{timestamp}_{commit}"
os.mkdir(dirname)

for method, nr in size_and_reps.items():
    for n, reps in nr:
        for periodic in [False, True]:
            print(f"=== {n=} {reps=} {periodic=}")

            cmd = [
                "./build/tests/mbd_benchmark",
                "--name", method,
                "-N", str(n),
                "--repetitions", str(reps),
            ]
            cmd.extend(["--fname", f"{dirname}/method={method}_n={n}_repetitions={reps}.csv"])
            if periodic:
                cmd.append("--periodic")

            print(cmd)

            subprocess.call(cmd, env={"OMP_NUM_THREADS": "1"} | os.environ)
