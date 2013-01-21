(set-logic QF_NRA)
(set-info :source | KeYmaera example: ternary4, node 54
Andre Platzer, Jan-David Quesel, and Philipp Rümmer. Real world verification. In Renate A. Schmidt, editor, International Conference on Automated Deduction, CADE'09, Montreal, Canada, Proceedings, volume 5663 of LNCS, pages 485(- 501.) Springer, 2009.
 |)
(set-info :smt-lib-version 2.0)
(declare-fun x () Real)
(declare-fun z () Real)
(declare-fun y () Real)
(assert (not (>= (+ (* (- 2) x) (* x x) (* 2 (* z z)) (* 2 (* z y)) (* 2 (* x x z)) (* (- 2) (* x z y)) (* (- 2) (* z z y)) (* x x x x) (* 2 (* z z y y))) (- 1))))
(check-sat)
