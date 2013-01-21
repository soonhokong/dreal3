(set-logic QF_NRA)
(declare-fun c2uscore2dollarsk!13 () Real)
(declare-fun x2uscore7dollarsk!8 () Real)
(declare-fun omuscore3dollarsk!12 () Real)
(declare-fun d1uscore7dollarsk!10 () Real)
(declare-fun c1uscore2dollarsk!14 () Real)
(declare-fun x1uscore7dollarsk!9 () Real)
(declare-fun d2uscore7dollarsk!11 () Real)
(declare-fun y2uscore7dollarsk!4 () Real)
(declare-fun e1uscore7dollarsk!6 () Real)
(declare-fun y1uscore7dollarsk!5 () Real)
(declare-fun e2uscore7dollarsk!7 () Real)
(declare-fun z2uscore6dollarsk!0 () Real)
(declare-fun f1uscore6dollarsk!2 () Real)
(declare-fun z1uscore6dollarsk!1 () Real)
(declare-fun f2uscore7dollarsk!3 () Real)
(declare-fun protectedzone () Real)
(declare-fun y2uscore3dollarsk!17 () Real)
(declare-fun x2uscore3dollarsk!19 () Real)
(declare-fun y1uscore3dollarsk!18 () Real)
(declare-fun x1uscore3dollarsk!20 () Real)
(declare-fun z2uscore3dollarsk!15 () Real)
(declare-fun z1uscore3dollarsk!16 () Real)
(declare-fun y2uscore2dollarsk!26 () Real)
(declare-fun x2uscore2dollarsk!24 () Real)
(declare-fun y1uscore2dollarsk!25 () Real)
(declare-fun x1uscore2dollarsk!23 () Real)
(declare-fun z2uscore2dollarsk!22 () Real)
(declare-fun z1uscore2dollarsk!21 () Real)
(declare-fun y2 () Real)
(declare-fun x2 () Real)
(declare-fun y1 () Real)
(declare-fun x1 () Real)
(declare-fun z2 () Real)
(declare-fun z1 () Real)
(assert (= d1uscore7dollarsk!10
           (* (- 1.0)
              omuscore3dollarsk!12
              (+ x2uscore7dollarsk!8 (* (- 1.0) c2uscore2dollarsk!13)))))
(assert (= d2uscore7dollarsk!11
           (* omuscore3dollarsk!12
              (+ x1uscore7dollarsk!9 (* (- 1.0) c1uscore2dollarsk!14)))))
(assert (= e1uscore7dollarsk!6
           (* (- 1.0)
              omuscore3dollarsk!12
              (+ y2uscore7dollarsk!4 (* (- 1.0) c2uscore2dollarsk!13)))))
(assert (= e2uscore7dollarsk!7
           (* omuscore3dollarsk!12
              (+ y1uscore7dollarsk!5 (* (- 1.0) c1uscore2dollarsk!14)))))
(assert (= f1uscore6dollarsk!2
           (* (- 1.0)
              omuscore3dollarsk!12
              (+ z2uscore6dollarsk!0 (* (- 1.0) c2uscore2dollarsk!13)))))
(assert (= f2uscore7dollarsk!3
           (* omuscore3dollarsk!12
              (+ z1uscore6dollarsk!1 (* (- 1.0) c1uscore2dollarsk!14)))))
(assert (>= (+ (* (+ x1uscore3dollarsk!20 (* (- 1.0) y1uscore3dollarsk!18))
                  (+ x1uscore3dollarsk!20 (* (- 1.0) y1uscore3dollarsk!18)))
               (* (+ x2uscore3dollarsk!19 (* (- 1.0) y2uscore3dollarsk!17))
                  (+ x2uscore3dollarsk!19 (* (- 1.0) y2uscore3dollarsk!17))))
            (* protectedzone protectedzone)))
(assert (>= (+ (* (+ y1uscore3dollarsk!18 (* (- 1.0) z1uscore3dollarsk!16))
                  (+ y1uscore3dollarsk!18 (* (- 1.0) z1uscore3dollarsk!16)))
               (* (+ y2uscore3dollarsk!17 (* (- 1.0) z2uscore3dollarsk!15))
                  (+ y2uscore3dollarsk!17 (* (- 1.0) z2uscore3dollarsk!15))))
            (* protectedzone protectedzone)))
(assert (>= (+ (* (+ x1uscore3dollarsk!20 (* (- 1.0) z1uscore3dollarsk!16))
                  (+ x1uscore3dollarsk!20 (* (- 1.0) z1uscore3dollarsk!16)))
               (* (+ x2uscore3dollarsk!19 (* (- 1.0) z2uscore3dollarsk!15))
                  (+ x2uscore3dollarsk!19 (* (- 1.0) z2uscore3dollarsk!15))))
            (* protectedzone protectedzone)))
(assert (>= (+ (* (+ x1uscore2dollarsk!23 (* (- 1.0) y1uscore2dollarsk!25))
                  (+ x1uscore2dollarsk!23 (* (- 1.0) y1uscore2dollarsk!25)))
               (* (+ x2uscore2dollarsk!24 (* (- 1.0) y2uscore2dollarsk!26))
                  (+ x2uscore2dollarsk!24 (* (- 1.0) y2uscore2dollarsk!26))))
            (* protectedzone protectedzone)))
(assert (>= (+ (* (+ y1uscore2dollarsk!25 (* (- 1.0) z1uscore2dollarsk!21))
                  (+ y1uscore2dollarsk!25 (* (- 1.0) z1uscore2dollarsk!21)))
               (* (+ y2uscore2dollarsk!26 (* (- 1.0) z2uscore2dollarsk!22))
                  (+ y2uscore2dollarsk!26 (* (- 1.0) z2uscore2dollarsk!22))))
            (* protectedzone protectedzone)))
(assert (>= (+ (* (+ x1uscore2dollarsk!23 (* (- 1.0) z1uscore2dollarsk!21))
                  (+ x1uscore2dollarsk!23 (* (- 1.0) z1uscore2dollarsk!21)))
               (* (+ x2uscore2dollarsk!24 (* (- 1.0) z2uscore2dollarsk!22))
                  (+ x2uscore2dollarsk!24 (* (- 1.0) z2uscore2dollarsk!22))))
            (* protectedzone protectedzone)))
(assert (>= (+ (* (+ x1 (* (- 1.0) y1)) (+ x1 (* (- 1.0) y1)))
               (* (+ x2 (* (- 1.0) y2)) (+ x2 (* (- 1.0) y2))))
            (* protectedzone protectedzone)))
(assert (>= (+ (* (+ y1 (* (- 1.0) z1)) (+ y1 (* (- 1.0) z1)))
               (* (+ y2 (* (- 1.0) z2)) (+ y2 (* (- 1.0) z2))))
            (* protectedzone protectedzone)))
(assert (>= (+ (* (+ x1 (* (- 1.0) z1)) (+ x1 (* (- 1.0) z1)))
               (* (+ x2 (* (- 1.0) z2)) (+ x2 (* (- 1.0) z2))))
            (* protectedzone protectedzone)))
(assert (not (>= (+ (* (+ x1uscore7dollarsk!9 (* (- 1.0) z1uscore6dollarsk!1))
                       (+ (* 2.0 d1uscore7dollarsk!10)
                          (* (- 2.0) f1uscore6dollarsk!2)))
                    (* (+ x2uscore7dollarsk!8 (* (- 1.0) z2uscore6dollarsk!0))
                       (+ (* 2.0 d2uscore7dollarsk!11)
                          (* (- 2.0) f2uscore7dollarsk!3))))
                 0.0)))
(check-sat)
(exit)
