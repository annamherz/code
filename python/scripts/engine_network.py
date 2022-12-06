import sys

eng = sys.argv[1]

try:
    file_path = sys.argv[2]
    print("using the folder of the file as the location to write the output file...")
    output_file = f"{file_path.rsplit('.',1)[0]}_{eng.lower()}.dat"

except:
    print("assuming the file path is just network_combined.dat as not provided...")
    file_path = "network_combined.dat"
    output_file = f"network_combined_{eng.lower()}.dat"

with open(file_path, "r") as file:
    with open(output_file, "w") as f:
            for line in file:
                    if eng in line:
                            f.write(f"{line}")

