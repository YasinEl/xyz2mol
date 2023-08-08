import os
import argparse

def rename_files(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if "MDtrj" in file and not file.endswith('.xyz'):
                old_path = os.path.join(root, file)
                new_path = old_path + '.xyz'
                os.rename(old_path, new_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Rename files with MDtrj in their name by adding .xyz extension')
    parser.add_argument('-directory', type=str, required=True, help='Directory to search for files recursively')

    args = parser.parse_args()

    rename_files(args.directory)
