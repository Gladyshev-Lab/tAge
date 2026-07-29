# Releasing tAge to Zenodo

How to publish a new tAge version to the **paper-cited** Zenodo record after a
GitHub release. This is a **manual, per-release** process — see "Why manual"
below.

## Key facts

- **Paper-cited concept DOI:** `10.5281/zenodo.19039289` — always resolves to
  the *latest* version. The paper cites the 1.0.0 version DOI
  (`10.5281/zenodo.19039290`); each version also gets its own frozen version DOI.
- **Clock models live in a separate record** (`10.5281/zenodo.18763485`).
  Releasing the R package does **not** touch it.
- Use **"New version"**, never "Edit": a published record's files are
  immutable; "Edit" only changes metadata.

## Prerequisites

1. Changes merged to `main`.
2. `DESCRIPTION` `Version:` bumped, `NEWS.md` updated.
3. A **GitHub release + tag** `X.Y.Z` created (web UI or
   `git tag -a X.Y.Z -m "tAge X.Y.Z" && git push --tags`). Note: tags have
   **no `v` prefix**.

## Steps

1. **Build the archive from the tag** (respects `.gitignore`, so the large
   `.h5ad` and other untracked files are excluded):

   ```bash
   git archive --format=zip -o tAge-X.Y.Z.zip X.Y.Z
   ```

   Sanity-check it does not contain junk:

   ```bash
   unzip -l tAge-X.Y.Z.zip | grep -Ei '\.h5ad|\.venv|\.pdf' && echo "REMOVE THESE" || echo "clean"
   ```

2. Open the current record: <https://zenodo.org/records/19039290> → **log in**
   (lab account) → click **"New version"**.

3. In the new draft: **remove** the previous `.zip`, **upload**
   `tAge-X.Y.Z.zip`.

4. Set **Version** = `X.Y.Z`. Update Publication date and, if useful, a one-line
   Description pointing to the GitHub release notes.

   > **DOI dialog:** leave the default **"No, I need one"** — Zenodo mints a
   > fresh version DOI for this release under the same concept (`19039289`) on
   > publish. Do **not** pick "Yes, I already have one" (that is for external
   > pre-existing DOIs; pasting the concept or old version DOI is wrong).
   > Reserving a DOI is optional.

5. Click **Publish**. The concept DOI `10.5281/zenodo.19039289` now points to
   `X.Y.Z`; a new version DOI is minted for this release.

6. Verify: open the concept DOI <https://doi.org/10.5281/zenodo.19039289> — it
   should show `X.Y.Z`.

## Do / Don't

- ✅ Always **New version** (keeps all prior version DOIs valid; concept DOI
  follows latest).
- ❌ **Do not delete/withdraw** old versions — that would break the version DOI
  cited in the paper.
- ❌ **Do not enable GitHub auto-sync** for this repo. The paper-cited record was
  created manually; auto-sync cannot append to concept `19039289` and would
  create a separate, fragmented DOI lineage. GitHub already serves the latest
  code at <https://github.com/Gladyshev-Lab/tAge>.

## Why manual

The paper cites a manually-created Zenodo record. Zenodo's GitHub integration
manages its own concept DOIs and cannot continue a manually-created concept, so
keeping `10.5281/zenodo.19039289` current requires a manual "New version" each
release. It takes a few minutes and keeps the citation stable.

## One-off: flag the buggy 1.0.0 version

(Only needed once, because 1.0.0 had correctness bugs — see `NEWS.md` 1.1.0.)
Open the **1.0.0** record → **Edit** (metadata only; does not change its DOI) →
prepend a warning to the Description pointing users to the latest version /
concept DOI. Save.
