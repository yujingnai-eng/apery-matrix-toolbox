# Publishing checklist

This directory is a **standalone public release candidate**, not the full notes vault.

## 1. Pre-flight

```sh
python3 scripts/check_public_tree.py
python3 run_tool.py --all
```

Both must pass.

Check manually:

- [ ] No third-party PDFs or referee reports.
- [ ] No personal notes, Obsidian config, or archive ZIPs.
- [ ] No absolute home-directory paths.
- [ ] No tokens, keys, or private URLs.
- [ ] `CITATION.cff` author metadata is correct.
- [ ] `LICENSE` and `LICENSE-docs` are the intended licenses.
- [ ] Repository name/description/visibility are correct.

## 2. Create a private GitHub repository first

```sh
cd ~/Desktop/apery-matrix-toolbox
git init
git add .
git commit -m "Apéry matrix toolbox v0.1.0"
gh repo create apery-matrix-toolbox --private --source=. --remote=origin --push
```

Private first is deliberate: it lets you inspect the rendered README and data
before making anything public.

## 3. Make it public

Only after the private repo looks right:

```sh
gh repo edit --visibility public --accept-visibility-change-consequences
```

Or use GitHub Settings → General → Danger Zone → Change visibility.

## 4. Release

```sh
git tag -a v0.1.0 -m "Apéry matrix toolbox v0.1.0"
git push origin v0.1.0
gh release create v0.1.0 --generate-notes
```

## 5. Optional DOI

Connect the GitHub release to Zenodo, then add the DOI to `CITATION.cff` and
the README. Do not mint a DOI before the author list is final.

## 6. What not to publish

Do not add:

- `03_历史归档/` or `02_离线审阅/`;
- `cm_order_audit/` 38 MB bundle;
- third-party PDFs;
- external referee reports;
- personal notes or absolute paths;
- credentials, tokens, or private URLs.

If you later want to publish more data, use a separate data repository or a
release asset, and keep the code repository small.
