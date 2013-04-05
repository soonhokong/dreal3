(set-logic QF_NRA_ODE)
(declare-fun z_0_1_t () Real)
(declare-fun z_0_1_0 () Real)
(declare-fun z_0_2_t () Real)
(declare-fun z_0_2_0 () Real)
(declare-fun z_1_1_t () Real)
(declare-fun z_1_1_0 () Real)
(declare-fun z_1_2_t () Real)
(declare-fun z_1_2_0 () Real)
(declare-fun z_2_1_t () Real)
(declare-fun z_2_1_0 () Real)
(declare-fun z_2_2_t () Real)
(declare-fun z_2_2_0 () Real)
(declare-fun z_3_1_t () Real)
(declare-fun z_3_1_0 () Real)
(declare-fun z_3_2_t () Real)
(declare-fun z_3_2_0 () Real)
(declare-fun y_0_1_t () Real)
(declare-fun y_0_1_0 () Real)
(declare-fun y_0_2_t () Real)
(declare-fun y_0_2_0 () Real)
(declare-fun y_1_1_t () Real)
(declare-fun y_1_1_0 () Real)
(declare-fun y_1_2_t () Real)
(declare-fun y_1_2_0 () Real)
(declare-fun y_2_1_t () Real)
(declare-fun y_2_1_0 () Real)
(declare-fun y_2_2_t () Real)
(declare-fun y_2_2_0 () Real)
(declare-fun y_3_1_t () Real)
(declare-fun y_3_1_0 () Real)
(declare-fun y_3_2_t () Real)
(declare-fun y_3_2_0 () Real)
(declare-fun x_0_1_t () Real)
(declare-fun x_0_1_0 () Real)
(declare-fun x_0_2_t () Real)
(declare-fun x_0_2_0 () Real)
(declare-fun x_1_1_t () Real)
(declare-fun x_1_1_0 () Real)
(declare-fun x_1_2_t () Real)
(declare-fun x_1_2_0 () Real)
(declare-fun x_2_1_t () Real)
(declare-fun x_2_1_0 () Real)
(declare-fun x_2_2_t () Real)
(declare-fun x_2_2_0 () Real)
(declare-fun x_3_1_t () Real)
(declare-fun x_3_1_0 () Real)
(declare-fun x_3_2_t () Real)
(declare-fun x_3_2_0 () Real)
(declare-fun v_0_1_t () Real)
(declare-fun v_0_1_0 () Real)
(declare-fun v_0_2_t () Real)
(declare-fun v_0_2_0 () Real)
(declare-fun v_1_1_t () Real)
(declare-fun v_1_1_0 () Real)
(declare-fun v_1_2_t () Real)
(declare-fun v_1_2_0 () Real)
(declare-fun v_2_1_t () Real)
(declare-fun v_2_1_0 () Real)
(declare-fun v_2_2_t () Real)
(declare-fun v_2_2_0 () Real)
(declare-fun v_3_1_t () Real)
(declare-fun v_3_1_0 () Real)
(declare-fun v_3_2_t () Real)
(declare-fun v_3_2_0 () Real)
(declare-fun time_1 () Real)
(declare-fun time_4 () Real)
(declare-fun time_3 () Real)
(declare-fun time_6 () Real)
(declare-fun time_5 () Real)
(declare-fun time_8 () Real)
(declare-fun time_7 () Real)
(define-ode 1 (= d/dt[x_0_1] (((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_0_1 / (z_0_1 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_0_1 / (z_0_1 + 0.5)))))) - (0.00005 * (1.0 - (z_0_1 / 20.0)))) * x_0_1) + (0.0 * x_0_1))))
(define-ode 1 (= d/dt[y_0_1] ((((0.00005 * (1.0 - (z_0_1 / 20.0))) * x_0_1) + (((0.0242 * (1.0 - (1.0 * (z_0_1 / 20.0)))) - 0.0168) * y_0_1)) + (0.0 * y_0_1))))
(define-ode 1 (= d/dt[z_0_1] (((0.0 - 1.0) * (z_0_1 / 62.5)) + (0.0 * z_0_1))))
(define-ode 1 (= d/dt[v_0_1] ((((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_0_1 / (z_0_1 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_0_1 / (z_0_1 + 0.5)))))) - (0.00005 * (1.0 - (z_0_1 / 20.0)))) * x_0_1) + (0.0 * x_0_1)) + ((((0.00005 * (1.0 - (z_0_1 / 20.0))) * x_0_1) + (((0.0242 * (1.0 - (1.0 * (z_0_1 / 20.0)))) - 0.0168) * y_0_1)) + (0.0 * y_0_1)))))
(define-ode 4 (= d/dt[x_1_2] (((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_1_2 / (z_1_2 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_1_2 / (z_1_2 + 0.5)))))) - (0.00005 * (1.0 - (z_1_2 / 20.0)))) * x_1_2) + (0.0 * x_1_2))))
(define-ode 4 (= d/dt[y_1_2] ((((0.00005 * (1.0 - (z_1_2 / 20.0))) * x_1_2) + (((0.0242 * (1.0 - (1.0 * (z_1_2 / 20.0)))) - 0.0168) * y_1_2)) + (0.0 * y_1_2))))
(define-ode 4 (= d/dt[z_1_2] (((20.0 - z_1_2) / 62.5) + (0.0 * z_1_2))))
(define-ode 4 (= d/dt[v_1_2] ((((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_1_2 / (z_1_2 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_1_2 / (z_1_2 + 0.5)))))) - (0.00005 * (1.0 - (z_1_2 / 20.0)))) * x_1_2) + (0.0 * x_1_2)) + ((((0.00005 * (1.0 - (z_1_2 / 20.0))) * x_1_2) + (((0.0242 * (1.0 - (1.0 * (z_1_2 / 20.0)))) - 0.0168) * y_1_2)) + (0.0 * y_1_2)))))
(define-ode 3 (= d/dt[x_1_1] (((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_1_1 / (z_1_1 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_1_1 / (z_1_1 + 0.5)))))) - (0.00005 * (1.0 - (z_1_1 / 20.0)))) * x_1_1) + (0.0 * x_1_1))))
(define-ode 3 (= d/dt[y_1_1] ((((0.00005 * (1.0 - (z_1_1 / 20.0))) * x_1_1) + (((0.0242 * (1.0 - (1.0 * (z_1_1 / 20.0)))) - 0.0168) * y_1_1)) + (0.0 * y_1_1))))
(define-ode 3 (= d/dt[z_1_1] (((0.0 - 1.0) * (z_1_1 / 62.5)) + (0.0 * z_1_1))))
(define-ode 3 (= d/dt[v_1_1] ((((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_1_1 / (z_1_1 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_1_1 / (z_1_1 + 0.5)))))) - (0.00005 * (1.0 - (z_1_1 / 20.0)))) * x_1_1) + (0.0 * x_1_1)) + ((((0.00005 * (1.0 - (z_1_1 / 20.0))) * x_1_1) + (((0.0242 * (1.0 - (1.0 * (z_1_1 / 20.0)))) - 0.0168) * y_1_1)) + (0.0 * y_1_1)))))
(define-ode 6 (= d/dt[x_2_2] (((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_2_2 / (z_2_2 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_2_2 / (z_2_2 + 0.5)))))) - (0.00005 * (1.0 - (z_2_2 / 20.0)))) * x_2_2) + (0.0 * x_2_2))))
(define-ode 6 (= d/dt[y_2_2] ((((0.00005 * (1.0 - (z_2_2 / 20.0))) * x_2_2) + (((0.0242 * (1.0 - (1.0 * (z_2_2 / 20.0)))) - 0.0168) * y_2_2)) + (0.0 * y_2_2))))
(define-ode 6 (= d/dt[z_2_2] (((20.0 - z_2_2) / 62.5) + (0.0 * z_2_2))))
(define-ode 6 (= d/dt[v_2_2] ((((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_2_2 / (z_2_2 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_2_2 / (z_2_2 + 0.5)))))) - (0.00005 * (1.0 - (z_2_2 / 20.0)))) * x_2_2) + (0.0 * x_2_2)) + ((((0.00005 * (1.0 - (z_2_2 / 20.0))) * x_2_2) + (((0.0242 * (1.0 - (1.0 * (z_2_2 / 20.0)))) - 0.0168) * y_2_2)) + (0.0 * y_2_2)))))
(define-ode 5 (= d/dt[x_2_1] (((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_2_1 / (z_2_1 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_2_1 / (z_2_1 + 0.5)))))) - (0.00005 * (1.0 - (z_2_1 / 20.0)))) * x_2_1) + (0.0 * x_2_1))))
(define-ode 5 (= d/dt[y_2_1] ((((0.00005 * (1.0 - (z_2_1 / 20.0))) * x_2_1) + (((0.0242 * (1.0 - (1.0 * (z_2_1 / 20.0)))) - 0.0168) * y_2_1)) + (0.0 * y_2_1))))
(define-ode 5 (= d/dt[z_2_1] (((0.0 - 1.0) * (z_2_1 / 62.5)) + (0.0 * z_2_1))))
(define-ode 5 (= d/dt[v_2_1] ((((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_2_1 / (z_2_1 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_2_1 / (z_2_1 + 0.5)))))) - (0.00005 * (1.0 - (z_2_1 / 20.0)))) * x_2_1) + (0.0 * x_2_1)) + ((((0.00005 * (1.0 - (z_2_1 / 20.0))) * x_2_1) + (((0.0242 * (1.0 - (1.0 * (z_2_1 / 20.0)))) - 0.0168) * y_2_1)) + (0.0 * y_2_1)))))
(define-ode 8 (= d/dt[x_3_2] (((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_3_2 / (z_3_2 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_3_2 / (z_3_2 + 0.5)))))) - (0.00005 * (1.0 - (z_3_2 / 20.0)))) * x_3_2) + (0.0 * x_3_2))))
(define-ode 8 (= d/dt[y_3_2] ((((0.00005 * (1.0 - (z_3_2 / 20.0))) * x_3_2) + (((0.0242 * (1.0 - (1.0 * (z_3_2 / 20.0)))) - 0.0168) * y_3_2)) + (0.0 * y_3_2))))
(define-ode 8 (= d/dt[z_3_2] (((20.0 - z_3_2) / 62.5) + (0.0 * z_3_2))))
(define-ode 8 (= d/dt[v_3_2] ((((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_3_2 / (z_3_2 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_3_2 / (z_3_2 + 0.5)))))) - (0.00005 * (1.0 - (z_3_2 / 20.0)))) * x_3_2) + (0.0 * x_3_2)) + ((((0.00005 * (1.0 - (z_3_2 / 20.0))) * x_3_2) + (((0.0242 * (1.0 - (1.0 * (z_3_2 / 20.0)))) - 0.0168) * y_3_2)) + (0.0 * y_3_2)))))
(define-ode 7 (= d/dt[x_3_1] (((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_3_1 / (z_3_1 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_3_1 / (z_3_1 + 0.5)))))) - (0.00005 * (1.0 - (z_3_1 / 20.0)))) * x_3_1) + (0.0 * x_3_1))))
(define-ode 7 (= d/dt[y_3_1] ((((0.00005 * (1.0 - (z_3_1 / 20.0))) * x_3_1) + (((0.0242 * (1.0 - (1.0 * (z_3_1 / 20.0)))) - 0.0168) * y_3_1)) + (0.0 * y_3_1))))
(define-ode 7 (= d/dt[z_3_1] (((0.0 - 1.0) * (z_3_1 / 62.5)) + (0.0 * z_3_1))))
(define-ode 7 (= d/dt[v_3_1] ((((((0.0204 * (0.0 + ((1.0 - 0.0) * (z_3_1 / (z_3_1 + 2.0))))) - (0.0076 * (0.8 + ((1.0 - 0.8) * (z_3_1 / (z_3_1 + 0.5)))))) - (0.00005 * (1.0 - (z_3_1 / 20.0)))) * x_3_1) + (0.0 * x_3_1)) + ((((0.00005 * (1.0 - (z_3_1 / 20.0))) * x_3_1) + (((0.0242 * (1.0 - (1.0 * (z_3_1 / 20.0)))) - 0.0168) * y_3_1)) + (0.0 * y_3_1)))))
(assert (<= 0.0 z_0_1_t))
(assert (<= z_0_1_t 100.0))
(assert (<= 0.0 z_0_1_0))
(assert (<= z_0_1_0 100.0))
(assert (<= 0.0 z_0_2_t))
(assert (<= z_0_2_t 100.0))
(assert (<= 0.0 z_0_2_0))
(assert (<= z_0_2_0 100.0))
(assert (<= 0.0 z_1_1_t))
(assert (<= z_1_1_t 100.0))
(assert (<= 0.0 z_1_1_0))
(assert (<= z_1_1_0 100.0))
(assert (<= 0.0 z_1_2_t))
(assert (<= z_1_2_t 100.0))
(assert (<= 0.0 z_1_2_0))
(assert (<= z_1_2_0 100.0))
(assert (<= 0.0 z_2_1_t))
(assert (<= z_2_1_t 100.0))
(assert (<= 0.0 z_2_1_0))
(assert (<= z_2_1_0 100.0))
(assert (<= 0.0 z_2_2_t))
(assert (<= z_2_2_t 100.0))
(assert (<= 0.0 z_2_2_0))
(assert (<= z_2_2_0 100.0))
(assert (<= 0.0 z_3_1_t))
(assert (<= z_3_1_t 100.0))
(assert (<= 0.0 z_3_1_0))
(assert (<= z_3_1_0 100.0))
(assert (<= 0.0 z_3_2_t))
(assert (<= z_3_2_t 100.0))
(assert (<= 0.0 z_3_2_0))
(assert (<= z_3_2_0 100.0))
(assert (<= 0.0 y_0_1_t))
(assert (<= y_0_1_t 100.0))
(assert (<= 0.0 y_0_1_0))
(assert (<= y_0_1_0 100.0))
(assert (<= 0.0 y_0_2_t))
(assert (<= y_0_2_t 100.0))
(assert (<= 0.0 y_0_2_0))
(assert (<= y_0_2_0 100.0))
(assert (<= 0.0 y_1_1_t))
(assert (<= y_1_1_t 100.0))
(assert (<= 0.0 y_1_1_0))
(assert (<= y_1_1_0 100.0))
(assert (<= 0.0 y_1_2_t))
(assert (<= y_1_2_t 100.0))
(assert (<= 0.0 y_1_2_0))
(assert (<= y_1_2_0 100.0))
(assert (<= 0.0 y_2_1_t))
(assert (<= y_2_1_t 100.0))
(assert (<= 0.0 y_2_1_0))
(assert (<= y_2_1_0 100.0))
(assert (<= 0.0 y_2_2_t))
(assert (<= y_2_2_t 100.0))
(assert (<= 0.0 y_2_2_0))
(assert (<= y_2_2_0 100.0))
(assert (<= 0.0 y_3_1_t))
(assert (<= y_3_1_t 100.0))
(assert (<= 0.0 y_3_1_0))
(assert (<= y_3_1_0 100.0))
(assert (<= 0.0 y_3_2_t))
(assert (<= y_3_2_t 100.0))
(assert (<= 0.0 y_3_2_0))
(assert (<= y_3_2_0 100.0))
(assert (<= 0.0 x_0_1_t))
(assert (<= x_0_1_t 100.0))
(assert (<= 0.0 x_0_1_0))
(assert (<= x_0_1_0 100.0))
(assert (<= 0.0 x_0_2_t))
(assert (<= x_0_2_t 100.0))
(assert (<= 0.0 x_0_2_0))
(assert (<= x_0_2_0 100.0))
(assert (<= 0.0 x_1_1_t))
(assert (<= x_1_1_t 100.0))
(assert (<= 0.0 x_1_1_0))
(assert (<= x_1_1_0 100.0))
(assert (<= 0.0 x_1_2_t))
(assert (<= x_1_2_t 100.0))
(assert (<= 0.0 x_1_2_0))
(assert (<= x_1_2_0 100.0))
(assert (<= 0.0 x_2_1_t))
(assert (<= x_2_1_t 100.0))
(assert (<= 0.0 x_2_1_0))
(assert (<= x_2_1_0 100.0))
(assert (<= 0.0 x_2_2_t))
(assert (<= x_2_2_t 100.0))
(assert (<= 0.0 x_2_2_0))
(assert (<= x_2_2_0 100.0))
(assert (<= 0.0 x_3_1_t))
(assert (<= x_3_1_t 100.0))
(assert (<= 0.0 x_3_1_0))
(assert (<= x_3_1_0 100.0))
(assert (<= 0.0 x_3_2_t))
(assert (<= x_3_2_t 100.0))
(assert (<= 0.0 x_3_2_0))
(assert (<= x_3_2_0 100.0))
(assert (<= 0.0 v_0_1_t))
(assert (<= v_0_1_t 100.0))
(assert (<= 0.0 v_0_1_0))
(assert (<= v_0_1_0 100.0))
(assert (<= 0.0 v_0_2_t))
(assert (<= v_0_2_t 100.0))
(assert (<= 0.0 v_0_2_0))
(assert (<= v_0_2_0 100.0))
(assert (<= 0.0 v_1_1_t))
(assert (<= v_1_1_t 100.0))
(assert (<= 0.0 v_1_1_0))
(assert (<= v_1_1_0 100.0))
(assert (<= 0.0 v_1_2_t))
(assert (<= v_1_2_t 100.0))
(assert (<= 0.0 v_1_2_0))
(assert (<= v_1_2_0 100.0))
(assert (<= 0.0 v_2_1_t))
(assert (<= v_2_1_t 100.0))
(assert (<= 0.0 v_2_1_0))
(assert (<= v_2_1_0 100.0))
(assert (<= 0.0 v_2_2_t))
(assert (<= v_2_2_t 100.0))
(assert (<= 0.0 v_2_2_0))
(assert (<= v_2_2_0 100.0))
(assert (<= 0.0 v_3_1_t))
(assert (<= v_3_1_t 100.0))
(assert (<= 0.0 v_3_1_0))
(assert (<= v_3_1_0 100.0))
(assert (<= 0.0 v_3_2_t))
(assert (<= v_3_2_t 100.0))
(assert (<= 0.0 v_3_2_0))
(assert (<= v_3_2_0 100.0))
(assert (<= 0.0 time_1))
(assert (<= time_1 1.0))
(assert (<= 0.0 time_4))
(assert (<= time_4 1.0))
(assert (<= 0.0 time_3))
(assert (<= time_3 1.0))
(assert (<= 0.0 time_6))
(assert (<= time_6 1.0))
(assert (<= 0.0 time_5))
(assert (<= time_5 1.0))
(assert (<= 0.0 time_8))
(assert (<= time_8 1.0))
(assert (<= 0.0 time_7))
(assert (<= time_7 1.0))
(assert (and (or (and (> (+ (+ (* (- (- (* 0.0204 (+ 0.0 (* (- 1.0 0.0) (/ z_2_2_t (+ z_2_2_t 2.0))))) (* 0.0076 (+ 0.8 (* (- 1.0 0.8) (/ z_2_2_t (+ z_2_2_t 0.5)))))) (* 0.00005 (- 1.0 (/ z_2_2_t 20.0)))) x_2_2_t) (* 0.0 x_2_2_t)) (+ (+ (* (* 0.00005 (- 1.0 (/ z_2_2_t 20.0))) x_2_2_t) (* (- (* 0.0242 (- 1.0 (* 1.0 (/ z_2_2_t 20.0)))) 0.0168) y_2_2_t)) (* 0.0 y_2_2_t))) 0.0) (< (+ x_2_2_t y_2_2_t) 10.0) (= v_3_1_0 v_2_2_t) (= z_3_1_0 z_2_2_t) (= y_3_1_0 y_2_2_t) (= x_3_1_0 x_2_2_t) (>= v_3_1_t 0.0) (>= z_3_1_t 0.0) (>= y_3_1_t 0.0) (>= x_3_1_t 0.0)) (and (> (+ (+ (* (- (- (* 0.0204 (+ 0.0 (* (- 1.0 0.0) (/ z_2_1_t (+ z_2_1_t 2.0))))) (* 0.0076 (+ 0.8 (* (- 1.0 0.8) (/ z_2_1_t (+ z_2_1_t 0.5)))))) (* 0.00005 (- 1.0 (/ z_2_1_t 20.0)))) x_2_1_t) (* 0.0 x_2_1_t)) (+ (+ (* (* 0.00005 (- 1.0 (/ z_2_1_t 20.0))) x_2_1_t) (* (- (* 0.0242 (- 1.0 (* 1.0 (/ z_2_1_t 20.0)))) 0.0168) y_2_1_t)) (* 0.0 y_2_1_t))) 0.0) (>= (+ x_2_1_t y_2_1_t) 10.0) (= v_3_2_0 v_2_1_t) (= z_3_2_0 z_2_1_t) (= y_3_2_0 y_2_1_t) (= x_3_2_0 x_2_1_t) (>= v_3_2_t 0.0) (>= z_3_2_t 0.0) (>= y_3_2_t 0.0) (>= x_3_2_t 0.0))) (or (and (> (+ (+ (* (- (- (* 0.0204 (+ 0.0 (* (- 1.0 0.0) (/ z_1_2_t (+ z_1_2_t 2.0))))) (* 0.0076 (+ 0.8 (* (- 1.0 0.8) (/ z_1_2_t (+ z_1_2_t 0.5)))))) (* 0.00005 (- 1.0 (/ z_1_2_t 20.0)))) x_1_2_t) (* 0.0 x_1_2_t)) (+ (+ (* (* 0.00005 (- 1.0 (/ z_1_2_t 20.0))) x_1_2_t) (* (- (* 0.0242 (- 1.0 (* 1.0 (/ z_1_2_t 20.0)))) 0.0168) y_1_2_t)) (* 0.0 y_1_2_t))) 0.0) (< (+ x_1_2_t y_1_2_t) 10.0) (= v_2_1_0 v_1_2_t) (= z_2_1_0 z_1_2_t) (= y_2_1_0 y_1_2_t) (= x_2_1_0 x_1_2_t) (>= v_2_1_t 0.0) (>= z_2_1_t 0.0) (>= y_2_1_t 0.0) (>= x_2_1_t 0.0)) (and (> (+ (+ (* (- (- (* 0.0204 (+ 0.0 (* (- 1.0 0.0) (/ z_1_1_t (+ z_1_1_t 2.0))))) (* 0.0076 (+ 0.8 (* (- 1.0 0.8) (/ z_1_1_t (+ z_1_1_t 0.5)))))) (* 0.00005 (- 1.0 (/ z_1_1_t 20.0)))) x_1_1_t) (* 0.0 x_1_1_t)) (+ (+ (* (* 0.00005 (- 1.0 (/ z_1_1_t 20.0))) x_1_1_t) (* (- (* 0.0242 (- 1.0 (* 1.0 (/ z_1_1_t 20.0)))) 0.0168) y_1_1_t)) (* 0.0 y_1_1_t))) 0.0) (>= (+ x_1_1_t y_1_1_t) 10.0) (= v_2_2_0 v_1_1_t) (= z_2_2_0 z_1_1_t) (= y_2_2_0 y_1_1_t) (= x_2_2_0 x_1_1_t) (>= v_2_2_t 0.0) (>= z_2_2_t 0.0) (>= y_2_2_t 0.0) (>= x_2_2_t 0.0))) (or (and (> (+ (+ (* (- (- (* 0.0204 (+ 0.0 (* (- 1.0 0.0) (/ z_0_2_t (+ z_0_2_t 2.0))))) (* 0.0076 (+ 0.8 (* (- 1.0 0.8) (/ z_0_2_t (+ z_0_2_t 0.5)))))) (* 0.00005 (- 1.0 (/ z_0_2_t 20.0)))) x_0_2_t) (* 0.0 x_0_2_t)) (+ (+ (* (* 0.00005 (- 1.0 (/ z_0_2_t 20.0))) x_0_2_t) (* (- (* 0.0242 (- 1.0 (* 1.0 (/ z_0_2_t 20.0)))) 0.0168) y_0_2_t)) (* 0.0 y_0_2_t))) 0.0) (< (+ x_0_2_t y_0_2_t) 10.0) (= v_1_1_0 v_0_2_t) (= z_1_1_0 z_0_2_t) (= y_1_1_0 y_0_2_t) (= x_1_1_0 x_0_2_t) (>= v_1_1_t 0.0) (>= z_1_1_t 0.0) (>= y_1_1_t 0.0) (>= x_1_1_t 0.0)) (and (> (+ (+ (* (- (- (* 0.0204 (+ 0.0 (* (- 1.0 0.0) (/ z_0_1_t (+ z_0_1_t 2.0))))) (* 0.0076 (+ 0.8 (* (- 1.0 0.8) (/ z_0_1_t (+ z_0_1_t 0.5)))))) (* 0.00005 (- 1.0 (/ z_0_1_t 20.0)))) x_0_1_t) (* 0.0 x_0_1_t)) (+ (+ (* (* 0.00005 (- 1.0 (/ z_0_1_t 20.0))) x_0_1_t) (* (- (* 0.0242 (- 1.0 (* 1.0 (/ z_0_1_t 20.0)))) 0.0168) y_0_1_t)) (* 0.0 y_0_1_t))) 0.0) (>= (+ x_0_1_t y_0_1_t) 10.0) (= v_1_2_0 v_0_1_t) (= z_1_2_0 z_0_1_t) (= y_1_2_0 y_0_1_t) (= x_1_2_0 x_0_1_t) (>= v_1_2_t 0.0) (>= z_1_2_t 0.0) (>= y_1_2_t 0.0) (>= x_1_2_t 0.0))) (= v_0_1_0 6.0) (= z_0_1_0 30.0) (= y_0_1_0 1.0) (= x_0_1_0 5.0) (>= v_0_1_t 0.0) (>= z_0_1_t 0.0) (>= y_0_1_t 0.0) (>= x_0_1_t 0.0) (>= v_3_1_t 0.0) (>= z_3_1_t 0.0) (>= y_3_1_t 0.0) (>= x_3_1_t 0.0)))
(check-sat)
(exit)
