(set-logic QF_NRA)
(set-info :source | KeYmaera example: ternary2, node 44
Andre Platzer, Jan-David Quesel, and Philipp Rümmer. Real world verification. In Renate A. Schmidt, editor, International Conference on Automated Deduction, CADE'09, Montreal, Canada, Proceedings, volume 5663 of LNCS, pages 485(- 501.) Springer, 2009.
 |)
(set-info :smt-lib-version 2.0)
(declare-fun x () Real)
(declare-fun y () Real)
(declare-fun z () Real)
(assert (not (>= (+ (* (- 4) (* x y z)) (* x x x x) (* y y y y) (* z z z z)) (- 1))))
(check-sat)
