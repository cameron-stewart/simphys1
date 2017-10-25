(TeX-add-style-hook
 "sim1"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("scrartcl" "a4paper" "11pt" "bibtotoc")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("nag" "l2tabu" "orthodox") ("inputenc" "utf8") ("fontenc" "T1") ("amsmath" "intlimits") ("units" "ugly")))
   (add-to-list 'LaTeX-verbatim-environments-local "python")
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "nag"
    "scrartcl"
    "scrartcl11"
    "inputenc"
    "color"
    "fontenc"
    "lmodern"
    "amsmath"
    "hyperref"
    "grffile"
    "units"
    "url"
    "breakurl"
    "xspace"
    "xcolor"
    "booktabs"
    "listings"
    "bm"
    "courier"
    "graphicx")
   (TeX-add-symbols
    "pythonstyle")
   (LaTeX-add-labels
    "fig:cannonball1"
    "fig:cannonball2"
    "fig:cannonball3")
   (LaTeX-add-xcolor-definecolors
    "deepblue"
    "deepred"
    "deepgreen"))
 :latex)

