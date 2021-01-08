import os

from pathlib import Path

directory = Path('/home/ec2-user/shm_webpage')/'backend'/'data_generated'/'timehistoryfiles'/'results'

files = os.listdir(directory)
files = [f for f in files if not f.startswith('.')]

def remove_files(files):
    for file in files:
        path_to_file = directory/file
        print(directory)
        os.remove(path_to_file)

print(directory)

print(files)

remove_files(files)

print(files)