import os
import argparse


def rename_with_parent(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.xyz') or file.endswith('.out'):
                parent_dir = os.path.basename(root)
                new_name = f"{parent_dir}__{file}"
                old_path = os.path.join(root, file)
                new_path = os.path.join(root, new_name)

                os.rename(old_path, new_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Rename .xyz and .out files by prepending the parent directory name')
    parser.add_argument('-directory', type=str, required=True, help='Directory to search for files recursively')

    args = parser.parse_args()

    rename_with_parent(args.directory)
