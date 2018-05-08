(TeX-add-style-hook
 "ECE742Project"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "amsmath"
    "mathtools"
    "graphicx")
   (LaTeX-add-labels
    "SyEquation"
    "SzEquation"
    "DzEzRelation"
    "fig:RelErrSigmaxm3"
    "fig:RelErrm"
    "fig:RelErrSigmax"))
 :latex)

