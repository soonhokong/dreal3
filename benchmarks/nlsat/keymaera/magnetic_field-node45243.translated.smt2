(set-logic QF_NRA)
(declare-fun xuscore5dollarsk!1 () Real)
(declare-fun a () Real)
(declare-fun vyuscore5dollarsk!3 () Real)
(declare-fun vxuscore5dollarsk!2 () Real)
(declare-fun buscore2dollarsk!5 () Real)
(declare-fun yuscore5dollarsk!0 () Real)
(declare-fun stateuscore2dollarsk!4 () Real)
(declare-fun vx () Real)
(declare-fun vy () Real)
(declare-fun x () Real)
(declare-fun y () Real)
(declare-fun b () Real)
(declare-fun yuscore2dollarsk!6 () Real)
(declare-fun xuscore2dollarsk!7 () Real)
(assert (= vyuscore5dollarsk!3 (+ (- 2.0) (* a (+ (- 2.0) xuscore5dollarsk!1)))))
(assert (= (+ (* vxuscore5dollarsk!2 vxuscore5dollarsk!2)
              (* vyuscore5dollarsk!3 vyuscore5dollarsk!3))
           8.0))
(assert (= vxuscore5dollarsk!2
           (+ 2.0
              (* (- 1.0) a (+ 2.0 yuscore5dollarsk!0))
              (* 4.0 buscore2dollarsk!5 (+ 1.0 (* (- 1.0) a))))))
(assert (= stateuscore2dollarsk!4 2.0))
(assert (= vyuscore5dollarsk!3 (- 2.0)))
(assert (= vxuscore5dollarsk!2 2.0))
(assert (= (* a (+ xuscore5dollarsk!1 yuscore5dollarsk!0))
           (+ (* 4.0 buscore2dollarsk!5) (* (- 4.0) a buscore2dollarsk!5))))
(assert (= stateuscore2dollarsk!4 1.0))
(assert (= stateuscore2dollarsk!4 0.0))
(assert (= vx 2.0))
(assert (= vy (- 2.0)))
(assert (= x 0.0))
(assert (= y 0.0))
(assert (= b 0.0))
(assert (not (= (* a (+ xuscore2dollarsk!7 yuscore2dollarsk!6))
                (+ (* 4.0 buscore2dollarsk!5) (* (- 4.0) a buscore2dollarsk!5)))))
(assert (not (= vxuscore5dollarsk!2 (- 2.0))))
(check-sat)
(exit)
