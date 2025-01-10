import random
import argparse

def random_sample(input_file, output_prefix, sample_sizes, seed):
    """
    Randomly sample IDs from an input file with reproducibility.
    
    :param input_file: Path to the input file containing IDs (one per line)
    :param output_prefix: Prefix for the output files
    :param sample_sizes: List of sizes to sample
    :param seed: Random seed for reproducibility
    """
    # Set random seed for reproducibility
    random.seed(seed)
    
    # Load all IDs into a list
    with open(input_file, 'r') as f:
        ids = [line.strip() for line in f]
    
    # Check if requested sample sizes are valid
    if max(sample_sizes) > len(ids):
        raise ValueError("Sample size exceeds the number of available IDs.")
    
    # Generate sampled subsets and save them to files
    for size in sample_sizes:
        sampled_ids = random.sample(ids, size)
        output_file = f"{output_prefix}_{size}.txt"
        with open(output_file, 'w') as f:
            f.write("\n".join(sampled_ids))
        print(f"Generated file: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Randomly sample IDs from a list with reproducibility.")
    parser.add_argument("-i", "--input", required=True, help="Input file containing IDs (one per line).")
    parser.add_argument("-o", "--output", required=True, help="Prefix for the output files.")
    parser.add_argument("-s", "--sizes", required=True, nargs="+", type=int,
                        help="List of sizes to sample, e.g., 20000 40000 60000.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility (default: 42).")
    
    args = parser.parse_args()
    
    try:
        random_sample(args.input, args.output, args.sizes, args.seed)
    except ValueError as e:
        print(f"Error: {e}")

