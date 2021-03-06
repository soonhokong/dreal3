(set-logic QF_NRA)
(set-info :source |
SMT script generated by Ultimate Automizer [1,2].
Ultimate Automizer is a software verifier for C programs that implements an
automata-based approach [3].
The commands in this SMT scripts are used for a constraint-based synthesis
of invariants [4].

2016-04-30, Matthias Heizmann (heizmann@informatik.uni-freiburg.de)


[1] http://http://ultimate.informatik.uni-freiburg.de/automizer/
[2] Matthias Heizmann, Daniel Dietsch, Marius Greitschus, Jan Leike,
Betim Musa, Claus Schätzle, Andreas Podelski: Ultimate Automizer with
Two-track Proofs - (Competition Contribution). TACAS 2016: 950-953
[3] Matthias Heizmann, Jochen Hoenicke, Andreas Podelski: Software Model
Checking for People Who Love Automata. CAV 2013:36-52
[4] Michael Colon, Sriram Sankaranarayanan, Henny Sipma: Linear Invariant
Generation Using Non-linear Constraint Solving. CAV 2003: 420-432

|)
(set-info :smt-lib-version 2.5)
(set-info :category "industrial")
(set-info :status sat)
(declare-fun liipp_0_c () Real)
(declare-fun |liipp_0__(ite GateController_alarm 1 0)| () Real)
(declare-fun liipp_0__GateController_time () Real)
(declare-fun liipp_0__GateController_maxGatePosition () Real)
(declare-fun liipp_0__GateController_thousand () Real)
(declare-fun liipp_0__GateController_one () Real)
(declare-fun liipp_0__GateController_gate () Real)
(declare-fun liipp_1_c () Real)
(declare-fun |liipp_1__(ite GateController_alarm 1 0)| () Real)
(declare-fun liipp_1__GateController_time () Real)
(declare-fun liipp_1__GateController_maxGatePosition () Real)
(declare-fun liipp_1__GateController_thousand () Real)
(declare-fun liipp_1__GateController_one () Real)
(declare-fun liipp_1__GateController_gate () Real)
(declare-fun liipp_2_c () Real)
(declare-fun |liipp_2__(ite GateController_alarm 1 0)| () Real)
(declare-fun liipp_2__GateController_time () Real)
(declare-fun liipp_2__GateController_maxGatePosition () Real)
(declare-fun liipp_2__GateController_thousand () Real)
(declare-fun liipp_2__GateController_one () Real)
(declare-fun liipp_2__GateController_gate () Real)
(declare-fun liipp_3_c () Real)
(declare-fun |liipp_3__(ite GateController_alarm 1 0)| () Real)
(declare-fun liipp_3__GateController_time () Real)
(declare-fun liipp_3__GateController_maxGatePosition () Real)
(declare-fun liipp_3__GateController_thousand () Real)
(declare-fun liipp_3__GateController_one () Real)
(declare-fun liipp_3__GateController_gate () Real)
(declare-fun liipp_4_replace_0 () Real)
(declare-fun liipp_4_replace_1 () Real)
(declare-fun liipp_4_replace_2 () Real)
(declare-fun liipp_4_replace_3 () Real)
(declare-fun liipp_4_replace_4 () Real)
(declare-fun liipp_4_replace_5 () Real)
(declare-fun liipp_5_replace_0 () Real)
(declare-fun liipp_5_replace_1 () Real)
(declare-fun liipp_5_replace_2 () Real)
(declare-fun liipp_5_replace_3 () Real)
(declare-fun liipp_5_replace_4 () Real)
(declare-fun liipp_5_replace_5 () Real)
(declare-fun motzkin_205_0 () Real)
(declare-fun motzkin_205_1 () Real)
(declare-fun motzkin_205_2 () Real)
(declare-fun motzkin_205_3 () Real)
(declare-fun motzkin_205_4 () Real)
(declare-fun motzkin_205_5 () Real)
(declare-fun motzkin_205_6 () Real)
(declare-fun motzkin_205_7 () Real)
(declare-fun motzkin_205_8 () Real)
(declare-fun motzkin_205_9 () Real)
(declare-fun motzkin_205_10 () Real)
(declare-fun motzkin_205_11 () Real)
(declare-fun motzkin_206_0 () Real)
(declare-fun motzkin_206_1 () Real)
(declare-fun motzkin_206_2 () Real)
(declare-fun motzkin_206_3 () Real)
(declare-fun motzkin_206_4 () Real)
(declare-fun motzkin_206_5 () Real)
(declare-fun motzkin_206_6 () Real)
(declare-fun motzkin_206_7 () Real)
(declare-fun motzkin_206_8 () Real)
(declare-fun motzkin_206_9 () Real)
(declare-fun motzkin_206_10 () Real)
(declare-fun motzkin_206_11 () Real)
(assert (and (>= motzkin_205_0 0.0) (>= motzkin_205_1 0.0) (>= motzkin_205_2 0.0) (>= motzkin_205_3 0.0) (>= motzkin_205_4 0.0) (>= motzkin_205_5 0.0) (>= motzkin_205_6 0.0) (>= motzkin_205_7 0.0) (>= motzkin_205_8 0.0) (>= motzkin_205_9 0.0) (>= motzkin_205_10 0.0) (>= motzkin_205_11 0.0) (= (+ (* motzkin_205_0 (- 1.0)) motzkin_205_1 (* motzkin_205_11 (+ (* (- 1.0) liipp_2__GateController_gate) 0.0))) 0.0) (= (+ (* motzkin_205_2 (- 1.0)) motzkin_205_3 (* motzkin_205_11 (+ (* (- 1.0) liipp_2__GateController_thousand) 0.0))) 0.0) (= (+ (* motzkin_205_2 1000.0) (* motzkin_205_3 (- 1000.0)) (* motzkin_205_4 100.0) (* motzkin_205_5 (- 100.0)) (* motzkin_205_9 (- 1.0)) motzkin_205_10 (* motzkin_205_11 (+ (* (- 1.0) liipp_2__GateController_one) 0.0))) 0.0) (= (+ (* motzkin_205_4 (- 1.0)) motzkin_205_5 (* motzkin_205_11 (+ (* (- 1.0) liipp_2__GateController_maxGatePosition) 0.0))) 0.0) (= (+ (* motzkin_205_6 (- 1.0)) (* motzkin_205_11 (+ (* (- 1.0) |liipp_2__(ite GateController_alarm 1 0)|) 0.0))) 0.0) (= (+ motzkin_205_7 (* motzkin_205_8 (- 1.0)) (* motzkin_205_11 (+ (* (- 1.0) liipp_2__GateController_time) 0.0))) 0.0) (<= (+ motzkin_205_9 (* motzkin_205_10 (- 1.0)) (* motzkin_205_11 (+ (* (- 1.0) liipp_2_c) 0.0))) 0.0) (or (< (+ motzkin_205_9 (* motzkin_205_10 (- 1.0)) (* motzkin_205_11 (+ (* (- 1.0) liipp_2_c) 0.0))) 0.0) (> 0.0 0.0)) (>= motzkin_206_0 0.0) (>= motzkin_206_1 0.0) (>= motzkin_206_2 0.0) (>= motzkin_206_3 0.0) (>= motzkin_206_4 0.0) (>= motzkin_206_5 0.0) (>= motzkin_206_6 0.0) (>= motzkin_206_7 0.0) (>= motzkin_206_8 0.0) (>= motzkin_206_9 0.0) (>= motzkin_206_10 0.0) (>= motzkin_206_11 0.0) (= (+ (* motzkin_206_0 (- 1.0)) motzkin_206_1 (* motzkin_206_11 (+ (* (- 1.0) liipp_3__GateController_gate) 0.0))) 0.0) (= (+ (* motzkin_206_2 (- 1.0)) motzkin_206_3 (* motzkin_206_11 (+ (* (- 1.0) liipp_3__GateController_thousand) 0.0))) 0.0) (= (+ (* motzkin_206_2 1000.0) (* motzkin_206_3 (- 1000.0)) (* motzkin_206_4 100.0) (* motzkin_206_5 (- 100.0)) (* motzkin_206_9 (- 1.0)) motzkin_206_10 (* motzkin_206_11 (+ (* (- 1.0) liipp_3__GateController_one) 0.0))) 0.0) (= (+ (* motzkin_206_4 (- 1.0)) motzkin_206_5 (* motzkin_206_11 (+ (* (- 1.0) liipp_3__GateController_maxGatePosition) 0.0))) 0.0) (= (+ (* motzkin_206_6 (- 1.0)) (* motzkin_206_11 (+ (* (- 1.0) |liipp_3__(ite GateController_alarm 1 0)|) 0.0))) 0.0) (= (+ motzkin_206_7 (* motzkin_206_8 (- 1.0)) (* motzkin_206_11 (+ (* (- 1.0) liipp_3__GateController_time) 0.0))) 0.0) (<= (+ motzkin_206_9 (* motzkin_206_10 (- 1.0)) (* motzkin_206_11 (+ (* (- 1.0) liipp_3_c) 0.0))) 0.0) (or (< (+ motzkin_206_9 (* motzkin_206_10 (- 1.0))) 0.0) (> motzkin_206_11 0.0))))
(declare-fun liipp_7_replace_0 () Real)
(declare-fun liipp_7_replace_1 () Real)
(declare-fun liipp_7_replace_2 () Real)
(declare-fun liipp_7_replace_3 () Real)
(declare-fun motzkin_207_0 () Real)
(declare-fun motzkin_207_1 () Real)
(declare-fun motzkin_207_2 () Real)
(declare-fun motzkin_207_3 () Real)
(declare-fun motzkin_208_0 () Real)
(declare-fun motzkin_208_1 () Real)
(declare-fun motzkin_208_2 () Real)
(declare-fun motzkin_208_3 () Real)
(assert (and (>= motzkin_207_0 0.0) (>= motzkin_207_1 0.0) (>= motzkin_207_2 0.0) (>= motzkin_207_3 0.0) (= (+ motzkin_207_0 (* motzkin_207_1 (+ (* (- 1.0) liipp_0__GateController_time) 0.0)) (* motzkin_207_2 (+ (* 1.0 liipp_2__GateController_time) 0.0)) (* motzkin_207_3 (+ (* 1.0 liipp_3__GateController_time) 0.0))) 0.0) (= (+ (* motzkin_207_0 (- 1.0)) (* motzkin_207_1 (+ (* (- 1.0) liipp_0__GateController_thousand) 0.0)) (* motzkin_207_2 (+ (* 1.0 liipp_2__GateController_thousand) 0.0)) (* motzkin_207_3 (+ (* 1.0 liipp_3__GateController_thousand) 0.0))) 0.0) (= (+ (* motzkin_207_1 (+ (* (- 1.0) |liipp_0__(ite GateController_alarm 1 0)|) 0.0)) (* motzkin_207_2 (+ (* 1.0 |liipp_2__(ite GateController_alarm 1 0)|) 0.0)) (* motzkin_207_3 (+ (* 1.0 |liipp_3__(ite GateController_alarm 1 0)|) 0.0))) 0.0) (= (+ (* motzkin_207_1 (+ (* (- 1.0) liipp_0__GateController_maxGatePosition) 0.0)) (* motzkin_207_2 (+ (* 1.0 liipp_2__GateController_maxGatePosition) 0.0)) (* motzkin_207_3 (+ (* 1.0 liipp_3__GateController_maxGatePosition) 0.0))) 0.0) (= (+ (* motzkin_207_1 (+ (* (- 1.0) liipp_0__GateController_one) 0.0)) (* motzkin_207_2 (+ (* 1.0 liipp_2__GateController_one) 0.0)) (* motzkin_207_3 (+ (* 1.0 liipp_3__GateController_one) 0.0))) 0.0) (= (+ (* motzkin_207_1 (+ (* (- 1.0) liipp_0__GateController_gate) 0.0)) (* motzkin_207_2 (+ (* 1.0 liipp_2__GateController_gate) 0.0)) (* motzkin_207_3 (+ (* 1.0 liipp_3__GateController_gate) 0.0))) 0.0) (<= (+ (* motzkin_207_1 (+ (* (- 1.0) liipp_0_c) 0.0)) (* motzkin_207_2 (+ (* 1.0 liipp_2_c) 0.0)) (* motzkin_207_3 (+ (* 1.0 liipp_3_c) 0.0))) 0.0) (or (< (+ (* motzkin_207_1 (+ (* (- 1.0) liipp_0_c) 0.0)) (* motzkin_207_3 (+ (* 1.0 liipp_3_c) 0.0))) 0.0) (> motzkin_207_2 0.0)) (>= motzkin_208_0 0.0) (>= motzkin_208_1 0.0) (>= motzkin_208_2 0.0) (>= motzkin_208_3 0.0) (= (+ motzkin_208_0 (* motzkin_208_1 (+ (* (- 1.0) liipp_1__GateController_time) 0.0)) (* motzkin_208_2 (+ (* 1.0 liipp_2__GateController_time) 0.0)) (* motzkin_208_3 (+ (* 1.0 liipp_3__GateController_time) 0.0))) 0.0) (= (+ (* motzkin_208_0 (- 1.0)) (* motzkin_208_1 (+ (* (- 1.0) liipp_1__GateController_thousand) 0.0)) (* motzkin_208_2 (+ (* 1.0 liipp_2__GateController_thousand) 0.0)) (* motzkin_208_3 (+ (* 1.0 liipp_3__GateController_thousand) 0.0))) 0.0) (= (+ (* motzkin_208_1 (+ (* (- 1.0) |liipp_1__(ite GateController_alarm 1 0)|) 0.0)) (* motzkin_208_2 (+ (* 1.0 |liipp_2__(ite GateController_alarm 1 0)|) 0.0)) (* motzkin_208_3 (+ (* 1.0 |liipp_3__(ite GateController_alarm 1 0)|) 0.0))) 0.0) (= (+ (* motzkin_208_1 (+ (* (- 1.0) liipp_1__GateController_maxGatePosition) 0.0)) (* motzkin_208_2 (+ (* 1.0 liipp_2__GateController_maxGatePosition) 0.0)) (* motzkin_208_3 (+ (* 1.0 liipp_3__GateController_maxGatePosition) 0.0))) 0.0) (= (+ (* motzkin_208_1 (+ (* (- 1.0) liipp_1__GateController_one) 0.0)) (* motzkin_208_2 (+ (* 1.0 liipp_2__GateController_one) 0.0)) (* motzkin_208_3 (+ (* 1.0 liipp_3__GateController_one) 0.0))) 0.0) (= (+ (* motzkin_208_1 (+ (* (- 1.0) liipp_1__GateController_gate) 0.0)) (* motzkin_208_2 (+ (* 1.0 liipp_2__GateController_gate) 0.0)) (* motzkin_208_3 (+ (* 1.0 liipp_3__GateController_gate) 0.0))) 0.0) (<= (+ (* motzkin_208_1 (+ (* (- 1.0) liipp_1_c) 0.0)) (* motzkin_208_2 (+ (* 1.0 liipp_2_c) 0.0)) (* motzkin_208_3 (+ (* 1.0 liipp_3_c) 0.0))) 0.0) (or (< (* motzkin_208_3 (+ (* 1.0 liipp_3_c) 0.0)) 0.0) (> (+ motzkin_208_1 motzkin_208_2) 0.0))))
(declare-fun liipp_8_replace_0 () Real)
(declare-fun liipp_8_replace_1 () Real)
(declare-fun liipp_8_replace_2 () Real)
(declare-fun motzkin_209_0 () Real)
(declare-fun motzkin_209_1 () Real)
(declare-fun motzkin_209_2 () Real)
(declare-fun motzkin_209_3 () Real)
(declare-fun motzkin_210_0 () Real)
(declare-fun motzkin_210_1 () Real)
(declare-fun motzkin_210_2 () Real)
(declare-fun motzkin_210_3 () Real)
(assert (and (>= motzkin_209_0 0.0) (>= motzkin_209_1 0.0) (>= motzkin_209_2 0.0) (>= motzkin_209_3 0.0) (= (+ (* motzkin_209_0 (- 1.0)) (* motzkin_209_2 (+ (* 1.0 |liipp_0__(ite GateController_alarm 1 0)|) 0.0)) (* motzkin_209_3 (+ (* 1.0 |liipp_1__(ite GateController_alarm 1 0)|) 0.0))) 0.0) (= (+ (* motzkin_209_1 (- 1.0)) (* motzkin_209_2 (+ (* 1.0 liipp_0__GateController_gate) 0.0)) (* motzkin_209_3 (+ (* 1.0 liipp_1__GateController_gate) 0.0))) 0.0) (= (+ motzkin_209_1 (* motzkin_209_2 (+ (* 1.0 liipp_0__GateController_maxGatePosition) 0.0)) (* motzkin_209_3 (+ (* 1.0 liipp_1__GateController_maxGatePosition) 0.0))) 0.0) (= (+ (* motzkin_209_2 (+ (* 1.0 liipp_0__GateController_time) 0.0)) (* motzkin_209_3 (+ (* 1.0 liipp_1__GateController_time) 0.0))) 0.0) (= (+ (* motzkin_209_2 (+ (* 1.0 liipp_0__GateController_thousand) 0.0)) (* motzkin_209_3 (+ (* 1.0 liipp_1__GateController_thousand) 0.0))) 0.0) (= (+ (* motzkin_209_2 (+ (* 1.0 liipp_0__GateController_one) 0.0)) (* motzkin_209_3 (+ (* 1.0 liipp_1__GateController_one) 0.0))) 0.0) (<= (+ (* motzkin_209_1 (- 1.0)) (* motzkin_209_2 (+ (* 1.0 liipp_0_c) 0.0)) (* motzkin_209_3 (+ (* 1.0 liipp_1_c) 0.0))) 0.0) (or (< (+ (* motzkin_209_1 (- 1.0)) (* motzkin_209_3 (+ (* 1.0 liipp_1_c) 0.0))) 0.0) (> motzkin_209_2 0.0)) (>= motzkin_210_0 0.0) (>= motzkin_210_1 0.0) (>= motzkin_210_2 0.0) (>= motzkin_210_3 0.0) (= (+ (* motzkin_210_0 (- 1.0)) (* motzkin_210_2 (+ (* 1.0 |liipp_0__(ite GateController_alarm 1 0)|) 0.0)) (* motzkin_210_3 (+ (* 1.0 |liipp_1__(ite GateController_alarm 1 0)|) 0.0))) 0.0) (= (+ motzkin_210_1 (* motzkin_210_2 (+ (* 1.0 liipp_0__GateController_gate) 0.0)) (* motzkin_210_3 (+ (* 1.0 liipp_1__GateController_gate) 0.0))) 0.0) (= (+ (* motzkin_210_1 (- 1.0)) (* motzkin_210_2 (+ (* 1.0 liipp_0__GateController_maxGatePosition) 0.0)) (* motzkin_210_3 (+ (* 1.0 liipp_1__GateController_maxGatePosition) 0.0))) 0.0) (= (+ (* motzkin_210_2 (+ (* 1.0 liipp_0__GateController_time) 0.0)) (* motzkin_210_3 (+ (* 1.0 liipp_1__GateController_time) 0.0))) 0.0) (= (+ (* motzkin_210_2 (+ (* 1.0 liipp_0__GateController_thousand) 0.0)) (* motzkin_210_3 (+ (* 1.0 liipp_1__GateController_thousand) 0.0))) 0.0) (= (+ (* motzkin_210_2 (+ (* 1.0 liipp_0__GateController_one) 0.0)) (* motzkin_210_3 (+ (* 1.0 liipp_1__GateController_one) 0.0))) 0.0) (<= (+ (* motzkin_210_1 (- 1.0)) (* motzkin_210_2 (+ (* 1.0 liipp_0_c) 0.0)) (* motzkin_210_3 (+ (* 1.0 liipp_1_c) 0.0))) 0.0) (or (< (+ (* motzkin_210_1 (- 1.0)) (* motzkin_210_3 (+ (* 1.0 liipp_1_c) 0.0))) 0.0) (> motzkin_210_2 0.0))))
(assert (= |liipp_2__(ite GateController_alarm 1 0)| 0))
(assert (= liipp_2__GateController_gate 0))
(assert (= liipp_3__GateController_gate 0))
(assert (= liipp_2__GateController_one 0))
(assert (= liipp_1__GateController_time 0))
(assert (= |liipp_0__(ite GateController_alarm 1 0)| 0))
(assert (= liipp_3__GateController_one 0))
(assert (= liipp_0__GateController_thousand 0))
(assert (= liipp_0__GateController_maxGatePosition 0))
(assert (= liipp_0__GateController_time 0))
(assert (= liipp_0__GateController_one 0))
(assert (= |liipp_1__(ite GateController_alarm 1 0)| 0))
(assert (= liipp_1__GateController_maxGatePosition 0))
(assert (= liipp_1__GateController_thousand 0))
(assert (= |liipp_3__(ite GateController_alarm 1 0)| 0))
(assert (= liipp_2__GateController_maxGatePosition 0))
(assert (= liipp_0_c 0))
(assert (= liipp_1__GateController_one 0))
(assert (= liipp_1_c 0))
(assert (or (= liipp_1__GateController_gate 0) (= liipp_2__GateController_time 0) (= liipp_3_c 0) (= liipp_3__GateController_maxGatePosition 0)))
(assert (or (= liipp_3__GateController_maxGatePosition 0) (= liipp_2__GateController_time 0) (= liipp_3__GateController_time 0) (= liipp_3_c 0)))
(assert (or (= liipp_2__GateController_thousand 0) (= liipp_2_c 0) (= liipp_3__GateController_maxGatePosition 0) (= liipp_1__GateController_gate 0)))
(assert (or (= liipp_1__GateController_gate 0) (= liipp_3__GateController_maxGatePosition 0) (= liipp_3_c 0) (= liipp_2__GateController_thousand 0)))
(assert (or (= liipp_3__GateController_time 0) (= liipp_3__GateController_thousand 0) (= liipp_2_c 0) (= liipp_3_c 0)))
(assert (or (= liipp_1__GateController_gate 0) (= liipp_2_c 0) (= liipp_3__GateController_maxGatePosition 0) (= liipp_3__GateController_thousand 0)))
(assert (or (= liipp_2__GateController_time 0) (= liipp_2_c 0) (= liipp_1__GateController_gate 0) (= liipp_3__GateController_time 0)))
(assert (or (= liipp_2__GateController_time 0) (= liipp_2_c 0) (= liipp_3__GateController_time 0) (= liipp_3_c 0)))
(assert (or (= liipp_3__GateController_time 0) (= liipp_2__GateController_time 0) (= liipp_2_c 0) (= liipp_3_c 0)))
(assert (or (= liipp_2__GateController_thousand 0) (= liipp_2__GateController_time 0) (= liipp_2_c 0) (= liipp_3_c 0)))
(check-sat)
(exit)
