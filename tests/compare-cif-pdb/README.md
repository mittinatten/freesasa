# Test runner to compare CIF and PDB formats

Checks that stdout and stderr is the same for the two formats.
Downloads the files to `./data` the first time the script is run.
Written in F#. Run with

```
> dotnet run ../data/diverse_pdbs_2000.txt
```
