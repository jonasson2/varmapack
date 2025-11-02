# Always make PDF; enable SyncTeX for editor↔PDF jumps
$pdf_mode = 1;
$pdflatex = 'pdflatex -interaction=nonstopmode -synctex=1 %O %S';

# Use biber automatically if a .bcf is present, else fall back to bibtex
$bibtex_use = 2;

# macOS PDF viewer (Preview or Skim). Skim auto-reloads nicely:
# $pdf_previewer = 'open %O %S';                # Preview
$pdf_previewer = 'open -a Skim %S';             # Skim (recommended)
$pdf_update_method = 4;                          # Smooth auto-refresh for Skim

# Quiet biber
$biber = 'biber --quiet %O %S';