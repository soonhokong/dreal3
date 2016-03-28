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
(set-info :status unsat)
(declare-fun liipp_0_c () Real)
(declare-fun liipp_0__proc_a () Real)
(declare-fun liipp_1_c () Real)
(declare-fun liipp_1__proc_a () Real)
(declare-fun liipp_2_c () Real)
(declare-fun liipp_2__proc_a () Real)
(declare-fun liipp_3_c () Real)
(declare-fun liipp_3__proc_a () Real)
(declare-fun liipp_4_c () Real)
(declare-fun liipp_4__proc_a () Real)
(declare-fun liipp_5_c () Real)
(declare-fun liipp_5__proc_a () Real)
(declare-fun liipp_6_c () Real)
(declare-fun liipp_6__proc_a () Real)
(declare-fun liipp_7_c () Real)
(declare-fun liipp_7__proc_a () Real)
(declare-fun liipp_8_replace_0 () Real)
(declare-fun liipp_9_replace_0 () Real)
(declare-fun motzkin_5795_0 () Real)
(declare-fun motzkin_5795_1 () Real)
(declare-fun motzkin_5795_2 () Real)
(declare-fun motzkin_5795_3 () Real)
(declare-fun motzkin_5795_4 () Real)
(declare-fun motzkin_5795_5 () Real)
(declare-fun motzkin_5796_0 () Real)
(declare-fun motzkin_5796_1 () Real)
(declare-fun motzkin_5796_2 () Real)
(declare-fun motzkin_5796_3 () Real)
(declare-fun motzkin_5796_4 () Real)
(declare-fun motzkin_5796_5 () Real)
(declare-fun motzkin_5797_0 () Real)
(declare-fun motzkin_5797_1 () Real)
(declare-fun motzkin_5797_2 () Real)
(declare-fun motzkin_5797_3 () Real)
(declare-fun motzkin_5797_4 () Real)
(declare-fun motzkin_5797_5 () Real)
(declare-fun motzkin_5798_0 () Real)
(declare-fun motzkin_5798_1 () Real)
(declare-fun motzkin_5798_2 () Real)
(declare-fun motzkin_5798_3 () Real)
(declare-fun motzkin_5798_4 () Real)
(declare-fun motzkin_5798_5 () Real)
(declare-fun motzkin_5799_0 () Real)
(declare-fun motzkin_5799_1 () Real)
(declare-fun motzkin_5799_2 () Real)
(declare-fun motzkin_5799_3 () Real)
(declare-fun motzkin_5799_4 () Real)
(declare-fun motzkin_5799_5 () Real)
(declare-fun motzkin_5800_0 () Real)
(declare-fun motzkin_5800_1 () Real)
(declare-fun motzkin_5800_2 () Real)
(declare-fun motzkin_5800_3 () Real)
(declare-fun motzkin_5800_4 () Real)
(declare-fun motzkin_5800_5 () Real)
(declare-fun motzkin_5801_0 () Real)
(declare-fun motzkin_5801_1 () Real)
(declare-fun motzkin_5801_2 () Real)
(declare-fun motzkin_5801_3 () Real)
(declare-fun motzkin_5801_4 () Real)
(declare-fun motzkin_5801_5 () Real)
(declare-fun motzkin_5802_0 () Real)
(declare-fun motzkin_5802_1 () Real)
(declare-fun motzkin_5802_2 () Real)
(declare-fun motzkin_5802_3 () Real)
(declare-fun motzkin_5802_4 () Real)
(declare-fun motzkin_5802_5 () Real)
(declare-fun motzkin_5803_0 () Real)
(declare-fun motzkin_5803_1 () Real)
(declare-fun motzkin_5803_2 () Real)
(declare-fun motzkin_5803_3 () Real)
(declare-fun motzkin_5803_4 () Real)
(declare-fun motzkin_5803_5 () Real)
(declare-fun motzkin_5804_0 () Real)
(declare-fun motzkin_5804_1 () Real)
(declare-fun motzkin_5804_2 () Real)
(declare-fun motzkin_5804_3 () Real)
(declare-fun motzkin_5804_4 () Real)
(declare-fun motzkin_5804_5 () Real)
(declare-fun motzkin_5805_0 () Real)
(declare-fun motzkin_5805_1 () Real)
(declare-fun motzkin_5805_2 () Real)
(declare-fun motzkin_5805_3 () Real)
(declare-fun motzkin_5805_4 () Real)
(declare-fun motzkin_5805_5 () Real)
(declare-fun motzkin_5806_0 () Real)
(declare-fun motzkin_5806_1 () Real)
(declare-fun motzkin_5806_2 () Real)
(declare-fun motzkin_5806_3 () Real)
(declare-fun motzkin_5806_4 () Real)
(declare-fun motzkin_5806_5 () Real)
(declare-fun motzkin_5807_0 () Real)
(declare-fun motzkin_5807_1 () Real)
(declare-fun motzkin_5807_2 () Real)
(declare-fun motzkin_5807_3 () Real)
(declare-fun motzkin_5807_4 () Real)
(declare-fun motzkin_5807_5 () Real)
(declare-fun motzkin_5808_0 () Real)
(declare-fun motzkin_5808_1 () Real)
(declare-fun motzkin_5808_2 () Real)
(declare-fun motzkin_5808_3 () Real)
(declare-fun motzkin_5808_4 () Real)
(declare-fun motzkin_5808_5 () Real)
(declare-fun motzkin_5809_0 () Real)
(declare-fun motzkin_5809_1 () Real)
(declare-fun motzkin_5809_2 () Real)
(declare-fun motzkin_5809_3 () Real)
(declare-fun motzkin_5809_4 () Real)
(declare-fun motzkin_5809_5 () Real)
(declare-fun motzkin_5810_0 () Real)
(declare-fun motzkin_5810_1 () Real)
(declare-fun motzkin_5810_2 () Real)
(declare-fun motzkin_5810_3 () Real)
(declare-fun motzkin_5810_4 () Real)
(declare-fun motzkin_5810_5 () Real)
(assert (and (>= motzkin_5795_0 0.0) (>= motzkin_5795_1 0.0) (>= motzkin_5795_2 0.0) (>= motzkin_5795_3 0.0) (>= motzkin_5795_4 0.0) (>= motzkin_5795_5 0.0) (= (+ (* motzkin_5795_0 (- 1.0)) (* motzkin_5795_1 (- 1.0)) motzkin_5795_2 motzkin_5795_3) 0.0) (= (+ motzkin_5795_1 (* motzkin_5795_2 (- 1.0)) (* motzkin_5795_4 (+ (* (- 1.0) liipp_0__proc_a) 0.0)) (* motzkin_5795_5 (+ (* (- 1.0) liipp_4__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5795_1 (- 2.0)) (* motzkin_5795_2 2.0)) 0.0) (<= (+ (* motzkin_5795_4 (+ (* (- 1.0) liipp_0_c) 0.0)) (* motzkin_5795_5 (+ (* (- 1.0) liipp_4_c) 0.0))) 0.0) (or (< (+ (* motzkin_5795_4 (+ (* (- 1.0) liipp_0_c) 0.0)) (* motzkin_5795_5 (+ (* (- 1.0) liipp_4_c) 0.0))) 0.0) (> 0.0 0.0)) (>= motzkin_5796_0 0.0) (>= motzkin_5796_1 0.0) (>= motzkin_5796_2 0.0) (>= motzkin_5796_3 0.0) (>= motzkin_5796_4 0.0) (>= motzkin_5796_5 0.0) (= (+ (* motzkin_5796_0 (- 1.0)) (* motzkin_5796_1 (- 1.0)) motzkin_5796_2 motzkin_5796_3) 0.0) (= (+ motzkin_5796_1 (* motzkin_5796_2 (- 1.0)) (* motzkin_5796_4 (+ (* (- 1.0) liipp_0__proc_a) 0.0)) (* motzkin_5796_5 (+ (* (- 1.0) liipp_5__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5796_1 (- 2.0)) (* motzkin_5796_2 2.0)) 0.0) (<= (+ (* motzkin_5796_4 (+ (* (- 1.0) liipp_0_c) 0.0)) (* motzkin_5796_5 (+ (* (- 1.0) liipp_5_c) 0.0))) 0.0) (or (< (* motzkin_5796_4 (+ (* (- 1.0) liipp_0_c) 0.0)) 0.0) (> motzkin_5796_5 0.0)) (>= motzkin_5797_0 0.0) (>= motzkin_5797_1 0.0) (>= motzkin_5797_2 0.0) (>= motzkin_5797_3 0.0) (>= motzkin_5797_4 0.0) (>= motzkin_5797_5 0.0) (= (+ (* motzkin_5797_0 (- 1.0)) (* motzkin_5797_1 (- 1.0)) motzkin_5797_2 motzkin_5797_3) 0.0) (= (+ motzkin_5797_1 (* motzkin_5797_2 (- 1.0)) (* motzkin_5797_4 (+ (* (- 1.0) liipp_0__proc_a) 0.0)) (* motzkin_5797_5 (+ (* (- 1.0) liipp_6__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5797_1 (- 2.0)) (* motzkin_5797_2 2.0)) 0.0) (<= (+ (* motzkin_5797_4 (+ (* (- 1.0) liipp_0_c) 0.0)) (* motzkin_5797_5 (+ (* (- 1.0) liipp_6_c) 0.0))) 0.0) (or (< (+ (* motzkin_5797_4 (+ (* (- 1.0) liipp_0_c) 0.0)) (* motzkin_5797_5 (+ (* (- 1.0) liipp_6_c) 0.0))) 0.0) (> 0.0 0.0)) (>= motzkin_5798_0 0.0) (>= motzkin_5798_1 0.0) (>= motzkin_5798_2 0.0) (>= motzkin_5798_3 0.0) (>= motzkin_5798_4 0.0) (>= motzkin_5798_5 0.0) (= (+ (* motzkin_5798_0 (- 1.0)) (* motzkin_5798_1 (- 1.0)) motzkin_5798_2 motzkin_5798_3) 0.0) (= (+ motzkin_5798_1 (* motzkin_5798_2 (- 1.0)) (* motzkin_5798_4 (+ (* (- 1.0) liipp_0__proc_a) 0.0)) (* motzkin_5798_5 (+ (* (- 1.0) liipp_7__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5798_1 (- 2.0)) (* motzkin_5798_2 2.0)) 0.0) (<= (+ (* motzkin_5798_4 (+ (* (- 1.0) liipp_0_c) 0.0)) (* motzkin_5798_5 (+ (* (- 1.0) liipp_7_c) 0.0))) 0.0) (or (< (* motzkin_5798_4 (+ (* (- 1.0) liipp_0_c) 0.0)) 0.0) (> motzkin_5798_5 0.0)) (>= motzkin_5799_0 0.0) (>= motzkin_5799_1 0.0) (>= motzkin_5799_2 0.0) (>= motzkin_5799_3 0.0) (>= motzkin_5799_4 0.0) (>= motzkin_5799_5 0.0) (= (+ (* motzkin_5799_0 (- 1.0)) (* motzkin_5799_1 (- 1.0)) motzkin_5799_2 motzkin_5799_3) 0.0) (= (+ motzkin_5799_1 (* motzkin_5799_2 (- 1.0)) (* motzkin_5799_4 (+ (* (- 1.0) liipp_1__proc_a) 0.0)) (* motzkin_5799_5 (+ (* (- 1.0) liipp_4__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5799_1 (- 2.0)) (* motzkin_5799_2 2.0)) 0.0) (<= (+ (* motzkin_5799_4 (+ (* (- 1.0) liipp_1_c) 0.0)) (* motzkin_5799_5 (+ (* (- 1.0) liipp_4_c) 0.0))) 0.0) (or (< (* motzkin_5799_5 (+ (* (- 1.0) liipp_4_c) 0.0)) 0.0) (> motzkin_5799_4 0.0)) (>= motzkin_5800_0 0.0) (>= motzkin_5800_1 0.0) (>= motzkin_5800_2 0.0) (>= motzkin_5800_3 0.0) (>= motzkin_5800_4 0.0) (>= motzkin_5800_5 0.0) (= (+ (* motzkin_5800_0 (- 1.0)) (* motzkin_5800_1 (- 1.0)) motzkin_5800_2 motzkin_5800_3) 0.0) (= (+ motzkin_5800_1 (* motzkin_5800_2 (- 1.0)) (* motzkin_5800_4 (+ (* (- 1.0) liipp_1__proc_a) 0.0)) (* motzkin_5800_5 (+ (* (- 1.0) liipp_5__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5800_1 (- 2.0)) (* motzkin_5800_2 2.0)) 0.0) (<= (+ (* motzkin_5800_4 (+ (* (- 1.0) liipp_1_c) 0.0)) (* motzkin_5800_5 (+ (* (- 1.0) liipp_5_c) 0.0))) 0.0) (or (< 0.0 0.0) (> (+ motzkin_5800_4 motzkin_5800_5) 0.0)) (>= motzkin_5801_0 0.0) (>= motzkin_5801_1 0.0) (>= motzkin_5801_2 0.0) (>= motzkin_5801_3 0.0) (>= motzkin_5801_4 0.0) (>= motzkin_5801_5 0.0) (= (+ (* motzkin_5801_0 (- 1.0)) (* motzkin_5801_1 (- 1.0)) motzkin_5801_2 motzkin_5801_3) 0.0) (= (+ motzkin_5801_1 (* motzkin_5801_2 (- 1.0)) (* motzkin_5801_4 (+ (* (- 1.0) liipp_1__proc_a) 0.0)) (* motzkin_5801_5 (+ (* (- 1.0) liipp_6__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5801_1 (- 2.0)) (* motzkin_5801_2 2.0)) 0.0) (<= (+ (* motzkin_5801_4 (+ (* (- 1.0) liipp_1_c) 0.0)) (* motzkin_5801_5 (+ (* (- 1.0) liipp_6_c) 0.0))) 0.0) (or (< (* motzkin_5801_5 (+ (* (- 1.0) liipp_6_c) 0.0)) 0.0) (> motzkin_5801_4 0.0)) (>= motzkin_5802_0 0.0) (>= motzkin_5802_1 0.0) (>= motzkin_5802_2 0.0) (>= motzkin_5802_3 0.0) (>= motzkin_5802_4 0.0) (>= motzkin_5802_5 0.0) (= (+ (* motzkin_5802_0 (- 1.0)) (* motzkin_5802_1 (- 1.0)) motzkin_5802_2 motzkin_5802_3) 0.0) (= (+ motzkin_5802_1 (* motzkin_5802_2 (- 1.0)) (* motzkin_5802_4 (+ (* (- 1.0) liipp_1__proc_a) 0.0)) (* motzkin_5802_5 (+ (* (- 1.0) liipp_7__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5802_1 (- 2.0)) (* motzkin_5802_2 2.0)) 0.0) (<= (+ (* motzkin_5802_4 (+ (* (- 1.0) liipp_1_c) 0.0)) (* motzkin_5802_5 (+ (* (- 1.0) liipp_7_c) 0.0))) 0.0) (or (< 0.0 0.0) (> (+ motzkin_5802_4 motzkin_5802_5) 0.0)) (>= motzkin_5803_0 0.0) (>= motzkin_5803_1 0.0) (>= motzkin_5803_2 0.0) (>= motzkin_5803_3 0.0) (>= motzkin_5803_4 0.0) (>= motzkin_5803_5 0.0) (= (+ (* motzkin_5803_0 (- 1.0)) (* motzkin_5803_1 (- 1.0)) motzkin_5803_2 motzkin_5803_3) 0.0) (= (+ motzkin_5803_1 (* motzkin_5803_2 (- 1.0)) (* motzkin_5803_4 (+ (* (- 1.0) liipp_2__proc_a) 0.0)) (* motzkin_5803_5 (+ (* (- 1.0) liipp_4__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5803_1 (- 2.0)) (* motzkin_5803_2 2.0)) 0.0) (<= (+ (* motzkin_5803_4 (+ (* (- 1.0) liipp_2_c) 0.0)) (* motzkin_5803_5 (+ (* (- 1.0) liipp_4_c) 0.0))) 0.0) (or (< (+ (* motzkin_5803_4 (+ (* (- 1.0) liipp_2_c) 0.0)) (* motzkin_5803_5 (+ (* (- 1.0) liipp_4_c) 0.0))) 0.0) (> 0.0 0.0)) (>= motzkin_5804_0 0.0) (>= motzkin_5804_1 0.0) (>= motzkin_5804_2 0.0) (>= motzkin_5804_3 0.0) (>= motzkin_5804_4 0.0) (>= motzkin_5804_5 0.0) (= (+ (* motzkin_5804_0 (- 1.0)) (* motzkin_5804_1 (- 1.0)) motzkin_5804_2 motzkin_5804_3) 0.0) (= (+ motzkin_5804_1 (* motzkin_5804_2 (- 1.0)) (* motzkin_5804_4 (+ (* (- 1.0) liipp_2__proc_a) 0.0)) (* motzkin_5804_5 (+ (* (- 1.0) liipp_5__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5804_1 (- 2.0)) (* motzkin_5804_2 2.0)) 0.0) (<= (+ (* motzkin_5804_4 (+ (* (- 1.0) liipp_2_c) 0.0)) (* motzkin_5804_5 (+ (* (- 1.0) liipp_5_c) 0.0))) 0.0) (or (< (* motzkin_5804_4 (+ (* (- 1.0) liipp_2_c) 0.0)) 0.0) (> motzkin_5804_5 0.0)) (>= motzkin_5805_0 0.0) (>= motzkin_5805_1 0.0) (>= motzkin_5805_2 0.0) (>= motzkin_5805_3 0.0) (>= motzkin_5805_4 0.0) (>= motzkin_5805_5 0.0) (= (+ (* motzkin_5805_0 (- 1.0)) (* motzkin_5805_1 (- 1.0)) motzkin_5805_2 motzkin_5805_3) 0.0) (= (+ motzkin_5805_1 (* motzkin_5805_2 (- 1.0)) (* motzkin_5805_4 (+ (* (- 1.0) liipp_2__proc_a) 0.0)) (* motzkin_5805_5 (+ (* (- 1.0) liipp_6__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5805_1 (- 2.0)) (* motzkin_5805_2 2.0)) 0.0) (<= (+ (* motzkin_5805_4 (+ (* (- 1.0) liipp_2_c) 0.0)) (* motzkin_5805_5 (+ (* (- 1.0) liipp_6_c) 0.0))) 0.0) (or (< (+ (* motzkin_5805_4 (+ (* (- 1.0) liipp_2_c) 0.0)) (* motzkin_5805_5 (+ (* (- 1.0) liipp_6_c) 0.0))) 0.0) (> 0.0 0.0)) (>= motzkin_5806_0 0.0) (>= motzkin_5806_1 0.0) (>= motzkin_5806_2 0.0) (>= motzkin_5806_3 0.0) (>= motzkin_5806_4 0.0) (>= motzkin_5806_5 0.0) (= (+ (* motzkin_5806_0 (- 1.0)) (* motzkin_5806_1 (- 1.0)) motzkin_5806_2 motzkin_5806_3) 0.0) (= (+ motzkin_5806_1 (* motzkin_5806_2 (- 1.0)) (* motzkin_5806_4 (+ (* (- 1.0) liipp_2__proc_a) 0.0)) (* motzkin_5806_5 (+ (* (- 1.0) liipp_7__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5806_1 (- 2.0)) (* motzkin_5806_2 2.0)) 0.0) (<= (+ (* motzkin_5806_4 (+ (* (- 1.0) liipp_2_c) 0.0)) (* motzkin_5806_5 (+ (* (- 1.0) liipp_7_c) 0.0))) 0.0) (or (< (* motzkin_5806_4 (+ (* (- 1.0) liipp_2_c) 0.0)) 0.0) (> motzkin_5806_5 0.0)) (>= motzkin_5807_0 0.0) (>= motzkin_5807_1 0.0) (>= motzkin_5807_2 0.0) (>= motzkin_5807_3 0.0) (>= motzkin_5807_4 0.0) (>= motzkin_5807_5 0.0) (= (+ (* motzkin_5807_0 (- 1.0)) (* motzkin_5807_1 (- 1.0)) motzkin_5807_2 motzkin_5807_3) 0.0) (= (+ motzkin_5807_1 (* motzkin_5807_2 (- 1.0)) (* motzkin_5807_4 (+ (* (- 1.0) liipp_3__proc_a) 0.0)) (* motzkin_5807_5 (+ (* (- 1.0) liipp_4__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5807_1 (- 2.0)) (* motzkin_5807_2 2.0)) 0.0) (<= (+ (* motzkin_5807_4 (+ (* (- 1.0) liipp_3_c) 0.0)) (* motzkin_5807_5 (+ (* (- 1.0) liipp_4_c) 0.0))) 0.0) (or (< (* motzkin_5807_5 (+ (* (- 1.0) liipp_4_c) 0.0)) 0.0) (> motzkin_5807_4 0.0)) (>= motzkin_5808_0 0.0) (>= motzkin_5808_1 0.0) (>= motzkin_5808_2 0.0) (>= motzkin_5808_3 0.0) (>= motzkin_5808_4 0.0) (>= motzkin_5808_5 0.0) (= (+ (* motzkin_5808_0 (- 1.0)) (* motzkin_5808_1 (- 1.0)) motzkin_5808_2 motzkin_5808_3) 0.0) (= (+ motzkin_5808_1 (* motzkin_5808_2 (- 1.0)) (* motzkin_5808_4 (+ (* (- 1.0) liipp_3__proc_a) 0.0)) (* motzkin_5808_5 (+ (* (- 1.0) liipp_5__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5808_1 (- 2.0)) (* motzkin_5808_2 2.0)) 0.0) (<= (+ (* motzkin_5808_4 (+ (* (- 1.0) liipp_3_c) 0.0)) (* motzkin_5808_5 (+ (* (- 1.0) liipp_5_c) 0.0))) 0.0) (or (< 0.0 0.0) (> (+ motzkin_5808_4 motzkin_5808_5) 0.0)) (>= motzkin_5809_0 0.0) (>= motzkin_5809_1 0.0) (>= motzkin_5809_2 0.0) (>= motzkin_5809_3 0.0) (>= motzkin_5809_4 0.0) (>= motzkin_5809_5 0.0) (= (+ (* motzkin_5809_0 (- 1.0)) (* motzkin_5809_1 (- 1.0)) motzkin_5809_2 motzkin_5809_3) 0.0) (= (+ motzkin_5809_1 (* motzkin_5809_2 (- 1.0)) (* motzkin_5809_4 (+ (* (- 1.0) liipp_3__proc_a) 0.0)) (* motzkin_5809_5 (+ (* (- 1.0) liipp_6__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5809_1 (- 2.0)) (* motzkin_5809_2 2.0)) 0.0) (<= (+ (* motzkin_5809_4 (+ (* (- 1.0) liipp_3_c) 0.0)) (* motzkin_5809_5 (+ (* (- 1.0) liipp_6_c) 0.0))) 0.0) (or (< (* motzkin_5809_5 (+ (* (- 1.0) liipp_6_c) 0.0)) 0.0) (> motzkin_5809_4 0.0)) (>= motzkin_5810_0 0.0) (>= motzkin_5810_1 0.0) (>= motzkin_5810_2 0.0) (>= motzkin_5810_3 0.0) (>= motzkin_5810_4 0.0) (>= motzkin_5810_5 0.0) (= (+ (* motzkin_5810_0 (- 1.0)) (* motzkin_5810_1 (- 1.0)) motzkin_5810_2 motzkin_5810_3) 0.0) (= (+ motzkin_5810_1 (* motzkin_5810_2 (- 1.0)) (* motzkin_5810_4 (+ (* (- 1.0) liipp_3__proc_a) 0.0)) (* motzkin_5810_5 (+ (* (- 1.0) liipp_7__proc_a) 0.0))) 0.0) (= (+ (* motzkin_5810_1 (- 2.0)) (* motzkin_5810_2 2.0)) 0.0) (<= (+ (* motzkin_5810_4 (+ (* (- 1.0) liipp_3_c) 0.0)) (* motzkin_5810_5 (+ (* (- 1.0) liipp_7_c) 0.0))) 0.0) (or (< 0.0 0.0) (> (+ motzkin_5810_4 motzkin_5810_5) 0.0))))
(declare-fun motzkin_5811_0 () Real)
(declare-fun motzkin_5811_1 () Real)
(declare-fun motzkin_5811_2 () Real)
(declare-fun motzkin_5811_3 () Real)
(declare-fun motzkin_5811_4 () Real)
(declare-fun motzkin_5811_5 () Real)
(declare-fun motzkin_5812_0 () Real)
(declare-fun motzkin_5812_1 () Real)
(declare-fun motzkin_5812_2 () Real)
(declare-fun motzkin_5812_3 () Real)
(declare-fun motzkin_5812_4 () Real)
(declare-fun motzkin_5812_5 () Real)
(assert (and (>= motzkin_5811_0 0.0) (>= motzkin_5811_1 0.0) (>= motzkin_5811_2 0.0) (>= motzkin_5811_3 0.0) (>= motzkin_5811_4 0.0) (>= motzkin_5811_5 0.0) (= (+ (* motzkin_5811_0 (- 1.0)) motzkin_5811_1 (* motzkin_5811_2 (+ (* 1.0 liipp_0__proc_a) 0.0)) (* motzkin_5811_3 (+ (* 1.0 liipp_1__proc_a) 0.0)) (* motzkin_5811_4 (+ (* 1.0 liipp_2__proc_a) 0.0)) (* motzkin_5811_5 (+ (* 1.0 liipp_3__proc_a) 0.0))) 0.0) (<= (+ (* motzkin_5811_0 7.0) (* motzkin_5811_1 (- 7.0)) (* motzkin_5811_2 (+ (* 1.0 liipp_0_c) 0.0)) (* motzkin_5811_3 (+ (* 1.0 liipp_1_c) 0.0)) (* motzkin_5811_4 (+ (* 1.0 liipp_2_c) 0.0)) (* motzkin_5811_5 (+ (* 1.0 liipp_3_c) 0.0))) 0.0) (or (< (+ (* motzkin_5811_0 7.0) (* motzkin_5811_1 (- 7.0)) (* motzkin_5811_3 (+ (* 1.0 liipp_1_c) 0.0)) (* motzkin_5811_5 (+ (* 1.0 liipp_3_c) 0.0))) 0.0) (> (+ motzkin_5811_2 motzkin_5811_4) 0.0)) (>= motzkin_5812_0 0.0) (>= motzkin_5812_1 0.0) (>= motzkin_5812_2 0.0) (>= motzkin_5812_3 0.0) (>= motzkin_5812_4 0.0) (>= motzkin_5812_5 0.0) (= (+ (* motzkin_5812_0 (- 1.0)) motzkin_5812_1 (* motzkin_5812_2 (+ (* 1.0 liipp_4__proc_a) 0.0)) (* motzkin_5812_3 (+ (* 1.0 liipp_5__proc_a) 0.0)) (* motzkin_5812_4 (+ (* 1.0 liipp_6__proc_a) 0.0)) (* motzkin_5812_5 (+ (* 1.0 liipp_7__proc_a) 0.0))) 0.0) (<= (+ (* motzkin_5812_0 7.0) (* motzkin_5812_1 (- 7.0)) (* motzkin_5812_2 (+ (* 1.0 liipp_4_c) 0.0)) (* motzkin_5812_3 (+ (* 1.0 liipp_5_c) 0.0)) (* motzkin_5812_4 (+ (* 1.0 liipp_6_c) 0.0)) (* motzkin_5812_5 (+ (* 1.0 liipp_7_c) 0.0))) 0.0) (or (< (+ (* motzkin_5812_0 7.0) (* motzkin_5812_1 (- 7.0)) (* motzkin_5812_3 (+ (* 1.0 liipp_5_c) 0.0)) (* motzkin_5812_5 (+ (* 1.0 liipp_7_c) 0.0))) 0.0) (> (+ motzkin_5812_2 motzkin_5812_4) 0.0))))
(check-sat)
(exit)