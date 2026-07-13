# Minimal makefile for Sphinx documentation

SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build
PYTHONPATH    ?= ../build-py

html:
	@PYTHONPATH="$(PYTHONPATH)" $(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

open: html
	@open "$(BUILDDIR)/html/index.html"

help:
	@PYTHONPATH="$(PYTHONPATH)" $(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: html open help Makefile

%: Makefile
	@PYTHONPATH="$(PYTHONPATH)" $(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
