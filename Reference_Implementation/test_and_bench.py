import os
import subprocess

BASE_DIR = os.getcwd()

for entry in os.listdir(BASE_DIR):
    subdir = os.path.join(BASE_DIR, entry)

    if not os.path.isdir(subdir):
        continue

    print(f"\n===== Processing {entry} =====")

    try:
        print("Running make...")
        subprocess.run(
            f"cd {entry} && make",
            shell=True,
            executable="/bin/bash",
            check=True
        )

        if os.path.exists(os.path.join(subdir, "tests/test_pylon")):
            print("Running test_pylon...")
            subprocess.run(
                f"cd {entry} && ./tests/test_pylon",
                shell=True,
                check=True
            )

        if os.path.exists(os.path.join(subdir, "tests/test_sign")):
            print("Running test_sign...")
            subprocess.run(
                f"cd {entry} && ./tests/test_sign",
                shell=True,
                check=True
            )

        if os.path.exists(os.path.join(subdir, "tests/test_speed")):
            print("Running test_speed...")
            subprocess.run(
                f"cd {entry} && ./tests/test_speed",
                shell=True,
                check=True
            )

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed in {entry}: {e}")
