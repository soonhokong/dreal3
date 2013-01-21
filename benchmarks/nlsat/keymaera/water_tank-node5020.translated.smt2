(set-logic QF_NRA)
(declare-fun yuscore2dollarsk!2 () Real)
(declare-fun t1uscore0dollarsk!0 () Real)
(declare-fun stuscore2dollarsk!1 () Real)
(assert (or (not (>= t1uscore0dollarsk!0 0.0))
            (<= (+ t1uscore0dollarsk!0 yuscore2dollarsk!2) 10.0)))
(assert (>= t1uscore0dollarsk!0 0.0))
(assert (not (<= 10.0 yuscore2dollarsk!2)))
(assert (= stuscore2dollarsk!1 0.0))
(assert (>= yuscore2dollarsk!2 1.0))
(assert (<= yuscore2dollarsk!2 12.0))
(assert (not (= stuscore2dollarsk!1 1.0)))
(assert (not (= stuscore2dollarsk!1 3.0)))
(assert (not (<= (+ t1uscore0dollarsk!0 yuscore2dollarsk!2) 12.0)))
(assert (or (and (>= yuscore2dollarsk!2 10.0) (not (= yuscore2dollarsk!2 10.0)))
            (<= (+ t1uscore0dollarsk!0 yuscore2dollarsk!2) 10.0)))
(assert (or (not (>= t1uscore0dollarsk!0 0.0)) (<= yuscore2dollarsk!2 10.0)))
(check-sat)
(exit)
