#!/usr/bin/env python3

import os

def rename_mdtrj_files(directory):
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if "MDtrj" in filename and not filename.endswith('.xyz'):
                old_path = os.path.join(dirpath, filename)
                new_path = old_path + ".xyz"
                os.rename(old_path, new_path)
                print(f"Renamed {old_path} to {new_path}")

if __name__ == "__main__":
    directory = input("Enter the directory path: ")
    rename_mdtrj_files(directory)
