# parse_flagstat.py
import pandas as pd
import sys
import glob

def parse_flagstat(flagstat_files, output_file):
    summary_data = []

    for flagstat_file in flagstat_files:
        sample_name = flagstat_file.split('/')[-1].replace(".flagstat.txt", "")

        with open(flagstat_file, 'r') as f:
            lines = f.readlines()

            # Basic checks to ensure we have enough lines
            if len(lines) < 13:
                print(f"Error: {flagstat_file} does not have the expected number of lines.")
                continue

            try:
                total_reads = int(lines[0].split()[0])
                mapped_reads = int(lines[4].split()[0])

                # Check if the line has the expected format before accessing
                if '(' in lines[4] and '%' in lines[4]:
                    percent_mapped = float(lines[4].split('(')[1].split('%')[0])
                else:
                    percent_mapped = float('nan')  # Assign NaN if the format is unexpected

                properly_paired = int(lines[8].split()[0])
                if '(' in lines[8] and '%' in lines[8]:
                    percent_properly_paired = float(lines[8].split('(')[1].split('%')[0])
                else:
                    percent_properly_paired = float('nan')

                singletons = int(lines[10].split()[0])
                if '(' in lines[10] and '%' in lines[10]:
                    percent_singletons = float(lines[10].split('(')[1].split('%')[0])
                else:
                    percent_singletons = float('nan')

                mate_diff_chr = int(lines[11].split()[0])
                mate_diff_chr_mapq5 = int(lines[12].split()[0])

            except (IndexError, ValueError) as e:
                print(f"Error parsing file {flagstat_file}: {e}")
                continue

            summary_data.append({
                "Sample": sample_name,
                "Total_Reads": total_reads,
                "Mapped_Reads": mapped_reads,
                "Percent_Mapped": percent_mapped,
                "Properly_Paired_Reads": properly_paired,
                "Percent_Properly_Paired": percent_properly_paired,
                "Singletons": singletons,
                "Percent_Singletons": percent_singletons,
                "Mate_Diff_Chr": mate_diff_chr,
                "Mate_Diff_Chr_(mapQ>=5)": mate_diff_chr_mapq5
            })

    df = pd.DataFrame(summary_data)
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    # First argument: output file path
    output_file = sys.argv[1]
    # Remaining arguments: input flagstat files
    flagstat_files = sys.argv[2:]
    parse_flagstat(flagstat_files, output_file)
