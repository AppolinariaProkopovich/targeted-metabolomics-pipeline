# Safe publication to GitHub

Run these commands from the repository root after installing and authenticating the GitHub CLI (`gh auth login`).

```bash
git init -b main

git add \
  .gitignore \
  DATA_POLICY.md \
  PUBLISH.md \
  README.md \
  R/metabolomics_analysis.R \
  config/config.example.R \
  data/README.md \
  results/README.md \
  run.sh \
  scripts/install_dependencies.R

git status --short
git commit -m "Add code-only metabolomics analysis pipeline"
gh repo create targeted-metabolomics-pipeline --public --source=. --remote=origin --push
```

Before committing, confirm that `git status --short` does **not** show `config/config.R`, files under `data/`, files under `results/`, or any CSV/PDF output.

No licence is included by default. Add one only after deciding how others may reuse the code.
