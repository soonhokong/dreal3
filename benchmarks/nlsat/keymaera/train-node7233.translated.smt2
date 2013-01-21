(set-logic QF_NRA)
(declare-fun MAeuscore2dollarsk!4 () Real)
(declare-fun zuscore2dollarsk!1 () Real)
(declare-fun MAduscore3dollarsk!8 () Real)
(declare-fun MAeuscore3dollarsk!7 () Real)
(declare-fun MAauscore3dollarsk!9 () Real)
(declare-fun b () Real)
(declare-fun vuscore2dollarsk!0 () Real)
(declare-fun MAauscore2dollarsk!6 () Real)
(declare-fun MAduscore2dollarsk!5 () Real)
(declare-fun Pduscore2dollarsk!2 () Real)
(declare-fun z () Real)
(declare-fun MAa () Real)
(declare-fun MAd () Real)
(declare-fun v () Real)
(declare-fun MAe () Real)
(declare-fun ep () Real)
(declare-fun amax () Real)
(declare-fun Pd () Real)
(declare-fun Pa () Real)
(declare-fun Pe () Real)
(declare-fun Pauscore2dollarsk!3 () Real)
(assert (<= zuscore2dollarsk!1 MAeuscore2dollarsk!4))
(assert (>= MAduscore3dollarsk!8 0.0))
(assert (<= MAauscore3dollarsk!9 MAeuscore3dollarsk!7))
(assert (<= (+ (* vuscore2dollarsk!0 vuscore2dollarsk!0)
               (* (- 1.0) MAduscore3dollarsk!8 MAduscore3dollarsk!8))
            (* 2.0 b (+ MAauscore3dollarsk!9 (* (- 1.0) zuscore2dollarsk!1)))))
(assert (<= (* vuscore2dollarsk!0 vuscore2dollarsk!0)
            (* 2.0 b (+ MAeuscore3dollarsk!7 (* (- 1.0) zuscore2dollarsk!1)))))
(assert (>= zuscore2dollarsk!1 MAauscore2dollarsk!6))
(assert (<= (+ (* vuscore2dollarsk!0 vuscore2dollarsk!0)
               (* (- 1.0) MAduscore2dollarsk!5 MAduscore2dollarsk!5))
            (* 2.0 b (+ MAauscore2dollarsk!6 (* (- 1.0) zuscore2dollarsk!1)))))
(assert (<= (* vuscore2dollarsk!0 vuscore2dollarsk!0)
            (* 2.0 b (+ MAeuscore2dollarsk!4 (* (- 1.0) zuscore2dollarsk!1)))))
(assert (>= vuscore2dollarsk!0 0.0))
(assert (>= MAduscore2dollarsk!5 0.0))
(assert (>= Pduscore2dollarsk!2 0.0))
(assert (<= MAauscore2dollarsk!6 MAeuscore2dollarsk!4))
(assert (<= (+ (* v v) (* (- 1.0) MAd MAd)) (* 2.0 b (+ MAa (* (- 1.0) z)))))
(assert (<= (* v v) (* 2.0 b (+ MAe (* (- 1.0) z)))))
(assert (>= v 0.0))
(assert (not (<= ep 0.0)))
(assert (not (<= b 0.0)))
(assert (not (<= amax 0.0)))
(assert (>= MAd 0.0))
(assert (>= Pd 0.0))
(assert (<= v Pd))
(assert (>= z Pa))
(assert (<= z Pe))
(assert (not (>= zuscore2dollarsk!1 Pauscore2dollarsk!3)))
(assert (not (<= vuscore2dollarsk!0 MAduscore2dollarsk!5)))
(check-sat)
(exit)
