# Zenodo archival and DOI for DECIMAP

The GitHub release v2.0.0 exists. A Zenodo DOI is minted by linking Zenodo to the
GitHub account; it cannot be created from the repository alone.

## What Zenodo does

Zenodo connects to GitHub through OAuth and watches the selected repository. Each
time a GitHub release is published, Zenodo stores a frozen snapshot of the
repository at that tag and mints a DOI. The repository can stay private, but the
Zenodo record and its snapshot are public. After the DOI exists, paste it into
README.md and the manuscript data availability statement.

## Steps

1. Sign in to https://zenodo.org (use the GitHub login).
2. Open the account menu and choose GitHub, or go to
   https://zenodo.org/account/settings/github/.
3. Authorize Zenodo to access GitHub if prompted.
4. Find `chaihu-nociceptive-mapping` in the repository list and turn on the toggle.
   For a private repository, confirm that the archived snapshot will be public.
5. Open the release entry and choose to archive the existing release v2.0.0
   (new releases are otherwise picked up automatically).
6. Edit the Zenodo record title, authors (Fengkai Liu, Tianhui Gao), and add the
   MIT license; save.
7. Copy the DOI (format 10.5281/zenodo.<number>) and the badge URL.
8. Replace the placeholder in README.md and the manuscript data availability
   section, then rebuild the Word package.

## Notes

- The DOI points to the immutable v2.0.0 snapshot. Later releases receive separate
  DOIs that Zenodo groups under a single concept DOI.
- If the repository must never be public even as a snapshot, Zenodo is not the
  option; use an internal archive and state that code is available on request.
