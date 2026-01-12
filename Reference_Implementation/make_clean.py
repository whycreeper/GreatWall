import os
import subprocess
import shlex

BASE_DIR = os.getcwd()

for entry in os.listdir(BASE_DIR):
    subdir = os.path.join(BASE_DIR, entry)

    if not os.path.isdir(subdir):
        continue

    print(f"\n===== Cleaning {entry} =====")

    try:
        subprocess.run(
            f"cd {shlex.quote(entry)} && make clean",
            shell=True,
            check=True
        )

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] make clean failed in {entry}: {e}")
