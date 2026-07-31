# tAge — documentation and artwork
#
#   make docs      build the pkgdown site into docs/
#   make serve     serve that build locally, to look at it before pushing
#   make check     validate _pkgdown.yml against the package's topics
#   make logo      regenerate logo-full.png + logo.png, sync them to tage-python
#   make favicons  regenerate pkgdown/favicon/ from man/figures/logo.png
#   make clean     remove docs/
#
# The published site is https://gladyshev-lab.github.io/tAge/, built by
# .github/workflows/pkgdown.yaml on every push to main. Building here is for
# previewing before you push - docs/ is git-ignored.

RSCRIPT ?= Rscript

# Predictions in the vignettes go through reticulate. This venv only needs
# joblib, pandas, scikit-learn and numpy (see inst/python/requirements.txt).
RETICULATE_PYTHON ?= $(CURDIR)/.venv/bin/python

# logo.py needs Pillow and numpy; pypandoc supplies a working pandoc, see below.
PYTHON ?= $(CURDIR)/../.venv/bin/python

# System pandoc can be broken (missing Haskell shared library), so R
# Markdown rendering borrows the binary that ships with pypandoc. Set PANDOC_DIR
# yourself to override, or to an empty value to just use whatever is on PATH.
PANDOC_DIR ?= $(shell $(PYTHON) -c \
	"import pypandoc, os; print(os.path.dirname(pypandoc.get_pandoc_path()))" 2>/dev/null)

R_ENV = $(if $(PANDOC_DIR),PATH="$(PANDOC_DIR):$$PATH" RSTUDIO_PANDOC="$(PANDOC_DIR)") \
	RETICULATE_PYTHON="$(RETICULATE_PYTHON)"

# Preview server. Bound to loopback by default: this is a throwaway server for
# one person, and Chronos is not the machine the browser runs on. Reach it with
# an SSH tunnel (printed by `make serve`), or set PREVIEW_BIND=0.0.0.0 if your
# browser can already reach this host directly.
PREVIEW_PORT ?= 9655
PREVIEW_BIND ?= 0.0.0.0

.PHONY: docs serve check logo favicons clean

docs:
	@echo "==> pkgdown: tAge"
	@rm -r $(CURDIR)/docs 2>/dev/null || true
	@$(R_ENV) $(RSCRIPT) -e 'pkgdown::build_site(preview = FALSE, install = FALSE)' >/dev/null
	@echo "    $(CURDIR)/docs/index.html"

serve:
	@test -f $(CURDIR)/docs/index.html || \
		{ echo "no docs/ yet — run 'make docs' first"; exit 1; }
	@echo "==> http://0.0.0.0:$(PREVIEW_PORT)/  (Ctrl-C to stop)"
	@echo "    tunnel from your machine:"
	@echo "      ssh -N -L $(PREVIEW_PORT):0.0.0.0:$(PREVIEW_PORT) $$USER@$$(hostname)"
	@cd $(CURDIR)/docs && $(PYTHON) -m http.server $(PREVIEW_PORT) --bind $(PREVIEW_BIND)

check:
	@$(R_ENV) $(RSCRIPT) -e 'pkgdown::check_pkgdown()'

logo:
	@$(PYTHON) logo.py

# Hits realfavicongenerator.net, so it needs network. Only rerun after the mark
# in man/figures/logo.png changes - the 7 files it writes are never hand-edited.
favicons:
	@$(R_ENV) $(RSCRIPT) -e 'pkgdown::build_favicons(overwrite = TRUE)'

clean:
	@rm -r $(CURDIR)/docs 2>/dev/null || true
	@echo "==> cleaned"
