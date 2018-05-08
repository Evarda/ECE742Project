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
    "fig:GaussianPulse"
    "fig:GaussianPropegation"
    "fig:GaussianBoundary"
    "fig:GaussianReflection"
    "fig:RelErrSigmaxm3"
    "fig:RelErrm"))
 :latex)

