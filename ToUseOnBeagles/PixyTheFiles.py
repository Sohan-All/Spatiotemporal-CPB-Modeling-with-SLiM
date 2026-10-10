import subprocess
import os
import sys

from sample_exclusions import EXCLUDED_SAMPLES, exclusion_id

# Only years whose sample set changed need re-running (2026-10-09: 2023 only, CLAUDE.md 7.9.16).
YEARS = ("2023",)


def write_pixy_popfile(year):
    """popFile<year> minus sample_exclusions.EXCLUDED_SAMPLES -> popFile<year>_pixy. pixy has no
    exclude option: a sample absent from the populations file is simply not used."""
    out = f"popFile{year}_pixy"
    dropped = []
    with open(f"popFile{year}") as fin, open(out, "w") as fout:
        for line in fin:
            sid = line.split()[0] if line.split() else ""
            if sid in EXCLUDED_SAMPLES:
                dropped.append(sid)
                continue
            fout.write(line)
    print(f"{year}: sample_exclusions [{exclusion_id()}] -> dropped {dropped or 'none'}")
    return out


def main():
    for year in YEARS:
        popFile = write_pixy_popfile(year)
        for i in range(1, 18):
            inVCF = f"chr{i}_cpb.vcf.gz"
            outFolder = f"statsChr{i}_{year}"
            prefix = f"chr{i}_{year}"
            print(f"calculating pi for chr {i} year {year}")
            os.makedirs(outFolder, exist_ok=True)
            # Remove the previous run's output first: if pixy fails, AverageData must find nothing
            # rather than silently re-read stale files (2026-10-09: stale 2023 output with H41 at n=7).
            for stat in ("pi", "dxy", "fst"):
                stale = os.path.join(outFolder, f"{prefix}_{stat}.txt")
                if os.path.exists(stale):
                    os.remove(stale)
            
            with open("temp.txt", 'w') as output_file:
                cmd = f"zcat {inVCF} | tail -n 1"
                subprocess.run(cmd, shell=True, stdout=output_file, check=True)

            ending_base = 0
            # read temp.txt's second line and take the second element
            with open('temp.txt', 'r') as f:
                lines = f.readlines()
            
                line = lines[0].strip()
                parts = line.split()
                ending_base = parts[1]
            
            subprocess.run(['pixy', '--stats', 'pi', 'dxy', 'fst', '--vcf', inVCF, '--populations', popFile, "--window_size", ending_base, "--n_cores", "8", "--output_folder", outFolder, "--output_prefix", prefix, "--bypass_invariant_check", "--silent"], check=True)
    
if __name__ == "__main__":
    main()
