(TeX-add-style-hook
 "ECE742Project"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "amsmath"
    "mathtools")
   (LaTeX-add-labels
    "DzEzRelation"))
 :latex)

