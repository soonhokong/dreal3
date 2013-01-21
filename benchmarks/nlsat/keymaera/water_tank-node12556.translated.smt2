(set-logic QF_NRA)
(declare-fun xuscore2dollarsk!3 () Real)
(declare-fun t10uscore0dollarsk!0 () Real)
(declare-fun stuscore2dollarsk!1 () Real)
(declare-fun yuscore2dollarsk!2 () Real)
(assert (or (not (>= t10uscore0dollarsk!0 0.0))
            (<= (+ t10uscore0dollarsk!0 xuscore2dollarsk!3) 2.0)))
(assert (>= t10uscore0dollarsk!0 0.0))
(assert (not (<= xuscore2dollarsk!3 2.0)))
(assert (= stuscore2dollarsk!1 1.0))
(assert (>= yuscore2dollarsk!2 1.0))
(assert (<= yuscore2dollarsk!2 12.0))
(assert (<= yuscore2dollarsk!2 (+ 10.0 xuscore2dollarsk!3)))
(assert (not (= stuscore2dollarsk!1 3.0)))
(assert (not (<= (+ t10uscore0dollarsk!0 yuscore2dollarsk!2) 12.0)))
(assert (or (and (>= xuscore2dollarsk!3 2.0) (not (= xuscore2dollarsk!3 2.0)))
            (<= (+ t10uscore0dollarsk!0 xuscore2dollarsk!3) 2.0)))
(assert (or (not (>= t10uscore0dollarsk!0 0.0)) (<= xuscore2dollarsk!3 2.0)))
(check-sat)
(exit)
