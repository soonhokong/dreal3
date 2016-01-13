(set-logic QF_NRA)
(declare-fun        a1 () Real [-1, 1])
(declare-fun forall x  () Real [-10, 10])
(assert
 (forall
  ((x Real)
   )
  (and
   (or (not (< a1 x)) (= x (abs x)))
   (or (not (< x a1)) (= (- 0.0 x) (abs x)))
   )))
(check-sat)
(exit)
