# Archiving the software on GitHub and Zenodo

Zenodo can automatically archive GitHub releases and mint a DOI for each release.
See https://help.zenodo.org/docs/github/archive-software/github-upload/

## Metadata files

- `CITATION.cff` enables the GitHub citation box. 
- `.zenodo.json` enables Zenodo software metadata files: https://help.zenodo.org/docs/github/describe-software/ If both files are present, Zenodo uses `.zenodo.json` and ignores `CITATION.cff` when both are present. 

Update the author list and affiliations in both files before creating the first public release.

## Creating a Zenodo archived release

Zenodo guide to archiving a GitHub release: https://help.zenodo.org/docs/github/archive-software/github-upload/ 

High level steps.

1. Connect your GitHub account to Zenodo and enable the repository.
2. Create a GitHub release tag, for example `v2026.01.02`.
3. Wait for Zenodo to process the release and mint the DOI.
4. Update `CITATION.cff` with the release DOI if you want it shown directly in GitHub's citation prompt.
