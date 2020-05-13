% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh4m1OL
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
% 
% Output:
% Ja_transl [3x8]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:04
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = palh4m1OL_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(8,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1OL_jacobia_transl_sym_varpar: qJ has to be [8x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh4m1OL_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1OL_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1OL_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->24)
	unknown=NaN(3,8);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	Ja_transl = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->26)
	unknown=NaN(3,8);
	t1 = sin(qJ(1));
	t3 = cos(qJ(1));
	unknown(1,1) = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	unknown(2,2) = 0.0e0;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = 0.0e0;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	Ja_transl = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->9), mult. (22->18), div. (0->0), fcn. (22->4), ass. (0->32)
	unknown=NaN(3,8);
	t1 = sin(qJ(1));
	t2 = cos(qJ(2));
	t3 = t1 * t2;
	t5 = sin(qJ(2));
	t6 = t1 * t5;
	t8 = cos(qJ(1));
	t12 = t8 * t5;
	t14 = t8 * t2;
	unknown(1,1) = t1 * pkin(7) - t3 * r_i_i_C(1) + t6 * r_i_i_C(2) + t8 * r_i_i_C(3);
	unknown(1,2) = -t12 * r_i_i_C(1) - t14 * r_i_i_C(2);
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = -t8 * pkin(7) + t14 * r_i_i_C(1) - t12 * r_i_i_C(2) + t1 * r_i_i_C(3);
	unknown(2,2) = -t6 * r_i_i_C(1) - t3 * r_i_i_C(2);
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t2 * r_i_i_C(1) - t5 * r_i_i_C(2);
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	Ja_transl = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (14->14), mult. (33->23), div. (0->0), fcn. (36->4), ass. (0->32)
	unknown=NaN(3,8);
	t1 = sin(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t5 = cos(qJ(1));
	t7 = cos(qJ(2));
	t8 = t1 * t7;
	t13 = t5 * t7;
	t15 = t5 * t2;
	unknown(1,1) = t1 * pkin(7) - t3 * r_i_i_C(1) - t5 * r_i_i_C(2) - t8 * r_i_i_C(3) - t8 * qJ(3);
	unknown(1,2) = t13 * r_i_i_C(1) - t15 * r_i_i_C(3) - t15 * qJ(3);
	unknown(1,3) = t13;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = -t5 * pkin(7) + t15 * r_i_i_C(1) - t1 * r_i_i_C(2) + t13 * r_i_i_C(3) + t13 * qJ(3);
	unknown(2,2) = t8 * r_i_i_C(1) - t3 * r_i_i_C(3) - t3 * qJ(3);
	unknown(2,3) = t8;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t2 * r_i_i_C(1) + t7 * r_i_i_C(3) + t7 * qJ(3);
	unknown(3,3) = t2;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	Ja_transl = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (33->24), mult. (87->35), div. (0->0), fcn. (106->6), ass. (0->44)
	unknown=NaN(3,8);
	t1 = sin(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t4 = sin(qJ(4));
	t6 = cos(qJ(2));
	t7 = t1 * t6;
	t8 = cos(qJ(4));
	t10 = t3 * t4 - t7 * t8;
	t14 = t3 * t8 + t7 * t4;
	t16 = cos(qJ(1));
	t21 = t16 * t2;
	t23 = t16 * t6;
	t25 = -t21 * t8 - t23 * t4;
	t26 = t25 * r_i_i_C(1);
	t29 = t21 * t4 - t23 * t8;
	t30 = t29 * r_i_i_C(2);
	t40 = -t14 * r_i_i_C(1);
	t41 = t10 * r_i_i_C(2);
	t48 = (-t2 * t4 + t6 * t8) * r_i_i_C(1);
	t52 = (-t2 * t8 - t6 * t4) * r_i_i_C(2);
	unknown(1,1) = t1 * pkin(7) + t10 * r_i_i_C(1) + t14 * r_i_i_C(2) + t16 * r_i_i_C(3) - t7 * qJ(3);
	unknown(1,2) = -t21 * qJ(3) + t26 + t30;
	unknown(1,3) = t23;
	unknown(1,4) = t26 + t30;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = -t16 * pkin(7) - t29 * r_i_i_C(1) + t25 * r_i_i_C(2) + t1 * r_i_i_C(3) + t23 * qJ(3);
	unknown(2,2) = -t3 * qJ(3) + t40 + t41;
	unknown(2,3) = t7;
	unknown(2,4) = t40 + t41;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t6 * qJ(3) + t48 + t52;
	unknown(3,3) = t2;
	unknown(3,4) = t48 + t52;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	Ja_transl = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (102->44), mult. (157->49), div. (0->0), fcn. (182->8), ass. (0->55)
	unknown=NaN(3,8);
	t1 = sin(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t4 = qJ(4) + qJ(5);
	t5 = sin(t4);
	t7 = cos(qJ(2));
	t8 = t1 * t7;
	t9 = cos(t4);
	t11 = t3 * t5 - t8 * t9;
	t15 = t3 * t9 + t8 * t5;
	t17 = cos(qJ(1));
	t19 = sin(qJ(4));
	t20 = t19 * pkin(3);
	t22 = cos(qJ(4));
	t23 = t22 * pkin(3);
	t28 = t17 * t2;
	t30 = t17 * t7;
	t32 = -t28 * t9 - t30 * t5;
	t33 = t32 * r_i_i_C(1);
	t36 = t28 * t5 - t30 * t9;
	t37 = t36 * r_i_i_C(2);
	t38 = t30 * t20;
	t39 = t28 * t23;
	t52 = -t15 * r_i_i_C(1);
	t53 = t11 * r_i_i_C(2);
	t54 = t8 * t20;
	t55 = t3 * t23;
	t63 = (-t2 * t5 + t7 * t9) * r_i_i_C(1);
	t67 = (-t2 * t9 - t7 * t5) * r_i_i_C(2);
	t69 = t2 * t19 * pkin(3);
	t71 = t7 * t22 * pkin(3);
	unknown(1,1) = t1 * pkin(7) + t11 * r_i_i_C(1) + t15 * r_i_i_C(2) + t17 * r_i_i_C(3) - t8 * qJ(3) + t3 * t20 - t8 * t23;
	unknown(1,2) = -t28 * qJ(3) + t33 + t37 - t38 - t39;
	unknown(1,3) = t30;
	unknown(1,4) = t33 + t37 - t39 - t38;
	unknown(1,5) = t33 + t37;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = -t17 * pkin(7) - t36 * r_i_i_C(1) + t32 * r_i_i_C(2) + t1 * r_i_i_C(3) + t30 * qJ(3) - t28 * t20 + t30 * t23;
	unknown(2,2) = -t3 * qJ(3) + t52 + t53 - t54 - t55;
	unknown(2,3) = t8;
	unknown(2,4) = t52 + t53 - t55 - t54;
	unknown(2,5) = t52 + t53;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t7 * qJ(3) + t63 + t67 - t69 + t71;
	unknown(3,3) = t2;
	unknown(3,4) = t63 + t67 + t71 - t69;
	unknown(3,5) = t63 + t67;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	Ja_transl = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:11
	% EndTime: 2020-04-11 23:04:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (217->74), mult. (315->81), div. (0->0), fcn. (376->10), ass. (0->69)
	unknown=NaN(3,8);
	t1 = sin(qJ(1));
	t2 = sin(qJ(2));
	t3 = t1 * t2;
	t4 = qJ(4) + qJ(5);
	t5 = cos(t4);
	t7 = cos(qJ(2));
	t8 = t1 * t7;
	t9 = sin(t4);
	t11 = t3 * t5 + t8 * t9;
	t12 = cos(qJ(6));
	t14 = cos(qJ(1));
	t15 = sin(qJ(6));
	t16 = t14 * t15;
	t20 = t14 * t12;
	t25 = t3 * t9 - t8 * t5;
	t28 = sin(qJ(4));
	t29 = t28 * pkin(3);
	t31 = cos(qJ(4));
	t32 = t31 * pkin(3);
	t37 = t14 * t7;
	t39 = t14 * t2;
	t41 = -t37 * t5 + t39 * t9;
	t43 = t41 * t12 * r_i_i_C(1);
	t45 = t41 * t15 * r_i_i_C(2);
	t48 = -t37 * t9 - t39 * t5;
	t49 = t48 * r_i_i_C(3);
	t50 = t48 * pkin(4);
	t51 = t37 * t29;
	t52 = t39 * t32;
	t59 = -t1 * t12 + t48 * t15;
	t63 = t1 * t15 + t48 * t12;
	t76 = t25 * t12 * r_i_i_C(1);
	t78 = t25 * t15 * r_i_i_C(2);
	t79 = -t11 * r_i_i_C(3);
	t80 = -t11 * pkin(4);
	t81 = t8 * t29;
	t82 = t3 * t32;
	t96 = -t2 * t5 - t7 * t9;
	t98 = t96 * t12 * r_i_i_C(1);
	t100 = t96 * t15 * r_i_i_C(2);
	t103 = -t2 * t9 + t7 * t5;
	t104 = t103 * r_i_i_C(3);
	t105 = t103 * pkin(4);
	t107 = t2 * t28 * pkin(3);
	t109 = t7 * t31 * pkin(3);
	unknown(1,1) = (-t11 * t12 - t16) * r_i_i_C(1) + (t11 * t15 - t20) * r_i_i_C(2) + t25 * r_i_i_C(3) + t25 * pkin(4) + t3 * t29 - t8 * t32 - t8 * qJ(3) + t1 * pkin(7);
	unknown(1,2) = -t39 * qJ(3) - t43 + t45 + t49 + t50 - t51 - t52;
	unknown(1,3) = t37;
	unknown(1,4) = -t43 + t45 + t49 + t50 - t52 - t51;
	unknown(1,5) = -t43 + t45 + t49 + t50;
	unknown(1,6) = t59 * r_i_i_C(1) + t63 * r_i_i_C(2);
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = -t41 * pkin(4) - t14 * pkin(7) - t63 * r_i_i_C(1) + t59 * r_i_i_C(2) - t41 * r_i_i_C(3) + t37 * qJ(3) - t39 * t29 + t37 * t32;
	unknown(2,2) = -t3 * qJ(3) - t76 + t78 + t79 + t80 - t81 - t82;
	unknown(2,3) = t8;
	unknown(2,4) = -t76 + t78 + t79 + t80 - t82 - t81;
	unknown(2,5) = -t76 + t78 + t79 + t80;
	unknown(2,6) = (-t11 * t15 + t20) * r_i_i_C(1) + (-t11 * t12 - t16) * r_i_i_C(2);
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t7 * qJ(3) + t100 + t104 + t105 - t107 + t109 - t98;
	unknown(3,3) = t2;
	unknown(3,4) = -t98 + t100 + t104 + t105 + t109 - t107;
	unknown(3,5) = -t98 + t100 + t104 + t105;
	unknown(3,6) = t103 * t15 * r_i_i_C(1) + t103 * t12 * r_i_i_C(2);
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	Ja_transl = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:12
	% EndTime: 2020-04-11 23:04:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->9), mult. (22->18), div. (0->0), fcn. (22->4), ass. (0->32)
	unknown=NaN(3,8);
	t1 = sin(qJ(1));
	t2 = sin(qJ(7));
	t3 = t1 * t2;
	t5 = cos(qJ(7));
	t6 = t1 * t5;
	t8 = cos(qJ(1));
	t12 = t8 * t5;
	t14 = t8 * t2;
	unknown(1,1) = -t1 * pkin(6) + t3 * r_i_i_C(1) + t6 * r_i_i_C(2) + t8 * r_i_i_C(3);
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = -t12 * r_i_i_C(1) + t14 * r_i_i_C(2);
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t8 * pkin(6) - t14 * r_i_i_C(1) - t12 * r_i_i_C(2) + t1 * r_i_i_C(3);
	unknown(2,2) = 0.0e0;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = -t6 * r_i_i_C(1) + t3 * r_i_i_C(2);
	unknown(2,8) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = 0.0e0;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = -t2 * r_i_i_C(1) - t5 * r_i_i_C(2);
	unknown(3,8) = 0.0e0;
	Ja_transl = unknown;
else
	Ja_transl=NaN(3,8);
end