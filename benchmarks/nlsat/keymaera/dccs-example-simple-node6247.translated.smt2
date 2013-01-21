(set-logic QF_NRA)
(declare-fun eps () Real)
(declare-fun B () Real)
(declare-fun v1uscore2dollarsk!3 () Real)
(declare-fun t4uscore0dollarsk!0 () Real)
(declare-fun v2 () Real)
(declare-fun x1uscore2dollarsk!2 () Real)
(declare-fun x2uscore2dollarsk!1 () Real)
(declare-fun x1 () Real)
(declare-fun x2 () Real)
(declare-fun v1 () Real)
(declare-fun A () Real)
(assert (or (= B 0.0)
            (and (or (= v1uscore2dollarsk!3 0.0)
                     (= B 0.0)
                     (and (not (>= B 0.0)) (not (<= v1uscore2dollarsk!3 0.0)))
                     (and (not (<= B 0.0)) (not (>= v1uscore2dollarsk!3 0.0))))
                 (not (= v1uscore2dollarsk!3 0.0)))
            (= B 0.0)
            (= (+ v1uscore2dollarsk!3 (* (- 1.0) B t4uscore0dollarsk!0)) 0.0)
            (and (not (<= (+ v1uscore2dollarsk!3
                             (* (- 1.0) B t4uscore0dollarsk!0))
                          0.0))
                 (not (<= B 0.0)))
            (and (not (>= (+ v1uscore2dollarsk!3
                             (* (- 1.0) B t4uscore0dollarsk!0))
                          0.0))
                 (not (>= B 0.0)))
            (and (<= B 0.0)
                 (or (= B 0.0)
                     (= (+ v1uscore2dollarsk!3 (* (- 1.0) B eps)) 0.0)
                     (and (not (>= (+ v1uscore2dollarsk!3 (* (- 1.0) B eps))
                                   0.0))
                          (not (<= B 0.0)))
                     (and (not (<= (+ v1uscore2dollarsk!3 (* (- 1.0) B eps))
                                   0.0))
                          (not (>= B 0.0))))
                 (not (= (+ v1uscore2dollarsk!3 (* (- 1.0) B eps)) 0.0)))))
(assert (>= t4uscore0dollarsk!0 0.0))
(assert (<= v1uscore2dollarsk!3 v2))
(assert (not (<= x2uscore2dollarsk!1 x1uscore2dollarsk!2)))
(assert (>= v1uscore2dollarsk!3 0.0))
(assert (not (<= x2 x1)))
(assert (>= v1 0.0))
(assert (>= v2 0.0))
(assert (<= v1 v2))
(assert (not (<= B 0.0)))
(assert (not (<= eps 0.0)))
(assert (not (<= A 0.0)))
(assert (not (= v1uscore2dollarsk!3 0.0)))
(assert (not (<= (+ v1uscore2dollarsk!3 (* A eps)) v2)))
(assert (not (>= (+ (* (- 1.0) B t4uscore0dollarsk!0) v1uscore2dollarsk!3) 0.0)))
(assert (or (not (>= t4uscore0dollarsk!0 0.0))
            (and (>= (+ v1uscore2dollarsk!3 (* (- 1.0) B t4uscore0dollarsk!0))
                     0.0)
                 (>= (+ eps (* (- 1.0) t4uscore0dollarsk!0)) 0.0))))
(assert (or (and (<= eps 0.0) (not (= eps 0.0)))
            (>= (+ eps (* (- 1.0) t4uscore0dollarsk!0)) 0.0)))
(assert (or (not (>= t4uscore0dollarsk!0 0.0))
            (and (>= v1uscore2dollarsk!3 0.0) (>= eps 0.0))))
(check-sat)
(exit)
