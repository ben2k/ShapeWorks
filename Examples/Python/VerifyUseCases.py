#!/usr/bin/env python3
import time
import subprocess

passed = []
failed = []

def run_case(use_case):
    command = f"python RunUseCase.py {use_case}"
    ret = subprocess.call(command.split())
    if ret:
        print(f"{use_case} : FAILED")
        failed.append(use_case)
    else:
        print(f"{use_case} : PASSED")
        passed.append(use_case)


start = time.time()
# run_case("ellipsoid --tiny_test")
# run_case("ellipsoid_mesh --tiny_test")
# run_case("ellipsoid_fd --tiny_test")
# run_case("ellipsoid_cut --tiny_test")
# run_case("lumps --tiny_test")
# run_case("left_atrium --tiny_test")
# run_case("femur --tiny_test")
# run_case("femur_cut --tiny_test")
# run_case("deep_ssm --tiny_test")
# run_case("supershapes_1mode_contour --tiny_test")
# run_case("thin_cavity_bean --tiny_test")
# run_case("ellipsoid_multiple_domain --tiny_test")
# run_case("ellipsoid_multiple_domain_mesh --tiny_test")
# run_case("ellipsoid_pca --tiny_test")

run_case("ellipsoid --verify")

end = time.time()

print("\n---------------------------------------------")
print("Testing Results:")
print("---------------------------------------------")
print(f"The following use cases passed ({len(passed)})")
for item in passed:
    print(f"{item} : PASSED")

if len(failed) > 0:
    print("\n---------------------------------------------")
    print(f"The following use cases failed ({len(failed)})")
    for item in failed:
        print(f"{item} : FAILED")

total = len(passed) + len(failed)
pass_percent = len(passed) / total

print("\n---------------------------------------------")
print(f"\n{pass_percent:.2%} tests passed, {len(failed)} failed out of {total}\n")


print(f"Total Test time: {end-start:.0f} seconds\n")
