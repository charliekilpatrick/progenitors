# Contributing to Progenitors

Thank you for your interest in contributing. This document covers how to submit issues, get in touch, and contribute code or documentation.

---

## Submitting issues

We use **GitHub Issues** for bug reports, feature requests, and documentation improvements.

### Before you open an issue

1. **Search existing issues** — Your question or bug may already be reported or answered. Search open and closed issues and discussions.
2. **Use the right place** — For general questions or usage help, consider opening a **Discussion** instead of an issue. Reserve issues for actionable bugs, feature proposals, or doc changes.

### When opening an issue, please include

**For bug reports:**

- A short, clear title.
- Steps to reproduce (commands or code).
- What you expected vs what actually happened.
- Your environment: OS, Python version (`python --version`), and how you installed the package (pip, conda, editable).
- Relevant error messages or tracebacks (use code blocks).
- If the bug depends on credentials or private data, describe the setup without sharing secrets.

**For feature requests or documentation:**

- A short summary and the goal (e.g. “Support reading catalog X” or “Clarify env vars in README”).
- Why it would be useful (use case or user story).
- For features: possible API or behavior, and whether you’re willing to implement it.

### Issue templates (optional)

If the repository has issue templates, choose the one that best fits (e.g. Bug report, Feature request) and fill in the sections. Otherwise, use the guidelines above and structure your description clearly.

---

## Contact

| Purpose | Contact |
|--------|---------|
| **Maintainer** | **Charlie Kilpatrick** — ckilpatrick@northwestern.edu |
| **Bug reports & features** | Prefer [GitHub Issues](https://github.com/YOUR_ORG/progenitors/issues) so others can benefit and track progress. |
| **Sensitive or private questions** | Email the maintainer directly. |
| **Security concerns** | Contact the maintainer by email; do not open a public issue. |

Replace `YOUR_ORG/progenitors` with the actual GitHub org/repo URL if different.

---

## Contributing code or documentation

### Getting started

1. **Fork the repository** on GitHub and clone your fork locally.
2. **Add the upstream remote** (optional but recommended):
   ```bash
   git remote add upstream https://github.com/YOUR_ORG/progenitors.git
   ```
3. **Create a branch** for your work (branch from `main` or the default branch):
   ```bash
   git checkout -b your-name/short-description
   ```

### Development setup

From the **repository root**:

```bash
pip install -e ".[dev]"
```

This installs the package in editable mode with dev dependencies (pytest, pytest-cov, Sphinx). Optional extras: `[sed]`, `[sed_hst]`, `[sed_fit]` as needed.

**Verify:** `progenitors --help` and `pytest tests/ -v`

### Before submitting a pull request

1. **Run the test suite**  
   - Full (excluding slow/remote):  
     `pytest -m "not slow and not remote" -v`  
   - With coverage:  
     `pytest -m "not slow and not remote" --cov=progenitors --cov-report=term-missing`
2. **Build and docs (if you changed code or docs)**  
   - `python -m build` (from repo root)  
   - `sphinx-build -b html docs docs/_build` (if you changed Sphinx sources)
3. **Keep the scope focused** — One logical change per PR (one bug fix or one feature). Large changes can be split into multiple PRs.
4. **Update documentation** — If you change behavior or add options, update the README, relevant module docs, or [docs/](docs/README.md) as appropriate.
5. **Sync with upstream** — Rebase or merge `main` into your branch and fix any conflicts:
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

### Submitting the pull request

1. Push your branch to your fork and open a **Pull Request** against `main` (or the repo’s default branch).
2. In the PR description:
   - Reference any related issue (e.g. “Fixes #123”).
   - Briefly describe the change and why it’s needed.
   - Note any breaking changes or migration steps.
3. Respond to review feedback and keep the PR up to date with the base branch if requested.

### Code and style

- **Python:** The project uses Python 3.12+. Follow existing style in the codebase (e.g. numpy-style docstrings where used). Running the test suite and linters (if configured) is expected.
- **Documentation:** Docstrings and README/docs should stay consistent with the rest of the repo (see main [README.md](README.md) and [docs/README.md](docs/README.md)).
- **Tests:** New behavior should be covered by tests under `tests/`. Use pytest markers `@pytest.mark.slow` or `@pytest.mark.remote` for tests that are slow or require network/external services.

---

## Summary

| I want to…           | Do this |
|----------------------|--------|
| Report a bug         | Open a [GitHub Issue](https://github.com/YOUR_ORG/progenitors/issues) with steps to reproduce and environment details. |
| Suggest a feature    | Open an issue or discussion describing the use case and (if applicable) proposed API. |
| Ask a question       | Open a Discussion or email the maintainer for private/sensitive topics. |
| Contribute code/docs | Fork, branch, develop, run tests, open a PR; see above for details. |
| Contact maintainer   | Charlie Kilpatrick — ckilpatrick@northwestern.edu |

Replace `YOUR_ORG/progenitors` in links with the actual repository URL.
