(set-logic QF_NRA)
(set-info :precision 0.01)
(declare-fun x () Real [2, 10])
(declare-fun y () Real [-10, 10])
(assert
        (and
                (= y (* (cos x) 2))
                (= y (log x))
        )
)
(check-sat)
(exit)
