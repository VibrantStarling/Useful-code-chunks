# Useful-code-chunks
Pieces of code that make life easier, but I would otherwise probably forget the day after I have worked it out

## File juggling
> How to mass rename files with specific prefixes in multiple different directories so file names are more informative. (Just incase you forget to add the option in programmes like metabat and have about 500 files that need prefixes)

```
# ADDING PREFIXES METABAT BINS IN DIFFERENT FOLDERS
# needs a list of directory names/prefixes you want to add to files further in. 
# E.g directory for SRA data named DRR215885 with metabat bins inside called bin.1.fa bin.2.fa ... etc. 
# This adds DRR215885 as a prefix to each bin to become DRR215885_bin.1.fa, DRR215885_bin.2.fa ... etc.

for i in $(cat names); do cd ${i}/*metabat*/; for file in $(ls *.fa); do  name="${i}_${file}"; mv -vi ./${file} ./${name};  done; cd ../../; done
```
