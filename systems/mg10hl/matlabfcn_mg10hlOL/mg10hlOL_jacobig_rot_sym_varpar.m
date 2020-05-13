% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% mg10hlOL
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% Jg_rot [3x13]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:06
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = mg10hlOL_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),uint8(0),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'mg10hlOL_jacobig_rot_sym_varpar: qJ has to be [13x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'mg10hlOL_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'mg10hlOL_jacobig_rot_sym_varpar: pkin has to be [16x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->39)
	unknown=NaN(3,13);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->39)
	unknown=NaN(3,13);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->41)
	unknown=NaN(3,13);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t2;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->41)
	unknown=NaN(3,13);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t2;
	unknown(2,3) = -t2;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:47
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (6->2), ass. (0->41)
	unknown=NaN(3,13);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = t1;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t2;
	unknown(2,3) = -t2;
	unknown(2,4) = -t2;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->41)
	unknown=NaN(3,13);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = t1;
	unknown(1,5) = t1;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t2;
	unknown(2,3) = -t2;
	unknown(2,4) = -t2;
	unknown(2,5) = -t2;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:49
	% EndTime: 2020-04-11 13:05:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (85->24), mult. (130->42), div. (0->0), fcn. (198->14), ass. (0->68)
	unknown=NaN(3,13);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	t3 = sin(qJ(2));
	t4 = t2 * t3;
	t5 = cos(pkin(14));
	t6 = cos(qJ(3));
	t8 = sin(pkin(14));
	t9 = sin(qJ(3));
	t11 = -t5 * t6 + t8 * t9;
	t13 = cos(qJ(2));
	t14 = t2 * t13;
	t17 = -t5 * t9 - t8 * t6;
	t19 = -t4 * t11 - t14 * t17;
	t20 = cos(pkin(16));
	t21 = -qJ(4) + pkin(15);
	t22 = cos(t21);
	t24 = sin(pkin(16));
	t25 = sin(t21);
	t27 = -t20 * t22 - t24 * t25;
	t31 = -t14 * t11 + t4 * t17;
	t34 = t20 * t25 - t24 * t22;
	t37 = cos(qJ(5));
	t42 = sin(qJ(5));
	t45 = t1 * t3;
	t47 = t1 * t13;
	t49 = -t45 * t11 - t47 * t17;
	t53 = -t47 * t11 + t45 * t17;
	t64 = t13 * t11 - t3 * t17;
	t68 = -t3 * t11 - t13 * t17;
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = t1;
	unknown(1,5) = t1;
	unknown(1,6) = (-(t19 * t27 + t31 * t34) * t37 - (-t19 * t34 + t31 * t27) * t42);
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t2;
	unknown(2,3) = -t2;
	unknown(2,4) = -t2;
	unknown(2,5) = -t2;
	unknown(2,6) = (-(t49 * t27 + t53 * t34) * t37 - (t53 * t27 - t49 * t34) * t42);
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = (-(t64 * t27 + t68 * t34) * t37 - (t68 * t27 - t64 * t34) * t42);
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:50
	% EndTime: 2020-04-11 13:05:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (169->30), mult. (265->53), div. (0->0), fcn. (395->16), ass. (0->76)
	unknown=NaN(3,13);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	t3 = sin(qJ(2));
	t4 = t2 * t3;
	t5 = cos(pkin(14));
	t6 = cos(qJ(3));
	t8 = sin(pkin(14));
	t9 = sin(qJ(3));
	t11 = -t5 * t6 + t8 * t9;
	t13 = cos(qJ(2));
	t14 = t2 * t13;
	t17 = -t5 * t9 - t8 * t6;
	t19 = -t4 * t11 - t14 * t17;
	t20 = cos(pkin(16));
	t21 = -qJ(4) + pkin(15);
	t22 = cos(t21);
	t24 = sin(pkin(16));
	t25 = sin(t21);
	t27 = -t20 * t22 - t24 * t25;
	t31 = -t14 * t11 + t4 * t17;
	t34 = t20 * t25 - t24 * t22;
	t36 = t19 * t27 + t31 * t34;
	t37 = cos(qJ(5));
	t41 = -t19 * t34 + t31 * t27;
	t42 = sin(qJ(5));
	t48 = sin(qJ(6));
	t50 = cos(qJ(6));
	t53 = t1 * t3;
	t55 = t1 * t13;
	t57 = -t53 * t11 - t55 * t17;
	t61 = -t55 * t11 + t53 * t17;
	t63 = t57 * t27 + t61 * t34;
	t67 = t61 * t27 - t57 * t34;
	t78 = t13 * t11 - t3 * t17;
	t82 = -t3 * t11 - t13 * t17;
	t84 = t78 * t27 + t82 * t34;
	t88 = t82 * t27 - t78 * t34;
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = t1;
	unknown(1,5) = t1;
	unknown(1,6) = (-t36 * t37 - t41 * t42);
	unknown(1,7) = (-(t36 * t42 - t41 * t37) * t48 + t1 * t50);
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t2;
	unknown(2,3) = -t2;
	unknown(2,4) = -t2;
	unknown(2,5) = -t2;
	unknown(2,6) = (-t63 * t37 - t67 * t42);
	unknown(2,7) = (-(-t67 * t37 + t63 * t42) * t48 - t2 * t50);
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = (-t84 * t37 - t88 * t42);
	unknown(3,7) = -((-t88 * t37 + t84 * t42) * t48);
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobig_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:50
	% EndTime: 2020-04-11 13:05:50
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (336->35), mult. (536->64), div. (0->0), fcn. (788->18), ass. (0->84)
	unknown=NaN(3,13);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	t3 = sin(qJ(2));
	t4 = t2 * t3;
	t5 = cos(pkin(14));
	t6 = cos(qJ(3));
	t8 = sin(pkin(14));
	t9 = sin(qJ(3));
	t11 = -t5 * t6 + t8 * t9;
	t13 = cos(qJ(2));
	t14 = t2 * t13;
	t17 = -t5 * t9 - t8 * t6;
	t19 = -t4 * t11 - t14 * t17;
	t20 = cos(pkin(16));
	t21 = -qJ(4) + pkin(15);
	t22 = cos(t21);
	t24 = sin(pkin(16));
	t25 = sin(t21);
	t27 = -t20 * t22 - t24 * t25;
	t31 = -t14 * t11 + t4 * t17;
	t34 = t20 * t25 - t24 * t22;
	t36 = t19 * t27 + t31 * t34;
	t37 = cos(qJ(5));
	t41 = -t19 * t34 + t31 * t27;
	t42 = sin(qJ(5));
	t44 = -t36 * t37 - t41 * t42;
	t47 = t36 * t42 - t41 * t37;
	t48 = sin(qJ(6));
	t50 = cos(qJ(6));
	t56 = sin(qJ(7));
	t58 = cos(qJ(7));
	t61 = t1 * t3;
	t63 = t1 * t13;
	t65 = -t61 * t11 - t63 * t17;
	t69 = -t63 * t11 + t61 * t17;
	t71 = t65 * t27 + t69 * t34;
	t75 = t69 * t27 - t65 * t34;
	t77 = -t71 * t37 - t75 * t42;
	t80 = -t75 * t37 + t71 * t42;
	t92 = t13 * t11 - t3 * t17;
	t96 = -t3 * t11 - t13 * t17;
	t98 = t92 * t27 + t96 * t34;
	t102 = t96 * t27 - t92 * t34;
	t104 = -t102 * t42 - t98 * t37;
	t107 = -t102 * t37 + t98 * t42;
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = t1;
	unknown(1,5) = t1;
	unknown(1,6) = t44;
	unknown(1,7) = (t1 * t50 - t47 * t48);
	unknown(1,8) = ((t1 * t48 + t47 * t50) * t56 + t44 * t58);
	unknown(1,9) = 0;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t2;
	unknown(2,3) = -t2;
	unknown(2,4) = -t2;
	unknown(2,5) = -t2;
	unknown(2,6) = t77;
	unknown(2,7) = (-t2 * t50 - t80 * t48);
	unknown(2,8) = ((-t2 * t48 + t80 * t50) * t56 + t77 * t58);
	unknown(2,9) = 0;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = t104;
	unknown(3,7) = -(t107 * t48);
	unknown(3,8) = (t107 * t50 * t56 + t104 * t58);
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobig_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->41)
	unknown=NaN(3,13);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = t1;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = t1;
	unknown(1,10) = 0;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t2;
	unknown(2,3) = -t2;
	unknown(2,4) = -t2;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = -t2;
	unknown(2,10) = 0;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobig_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->41)
	unknown=NaN(3,13);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = t1;
	unknown(1,11) = 0;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = -t2;
	unknown(2,11) = 0;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 11
	%% Symbolic Calculation
	% From jacobig_rot_11_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:49
	% EndTime: 2020-04-11 13:05:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (0->0), div. (0->0), fcn. (10->2), ass. (0->41)
	unknown=NaN(3,13);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = t1;
	unknown(1,3) = t1;
	unknown(1,4) = t1;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = t1;
	unknown(1,10) = 0;
	unknown(1,11) = t1;
	unknown(1,12) = 0;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = -t2;
	unknown(2,3) = -t2;
	unknown(2,4) = -t2;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = -t2;
	unknown(2,10) = 0;
	unknown(2,11) = -t2;
	unknown(2,12) = 0;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
elseif link_index == 12
	%% Symbolic Calculation
	% From jacobig_rot_12_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 13:05:48
	% EndTime: 2020-04-11 13:05:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->41)
	unknown=NaN(3,13);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(1,6) = 0;
	unknown(1,7) = 0;
	unknown(1,8) = 0;
	unknown(1,9) = 0;
	unknown(1,10) = t1;
	unknown(1,11) = 0;
	unknown(1,12) = t1;
	unknown(1,13) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(2,6) = 0;
	unknown(2,7) = 0;
	unknown(2,8) = 0;
	unknown(2,9) = 0;
	unknown(2,10) = -t2;
	unknown(2,11) = 0;
	unknown(2,12) = -t2;
	unknown(2,13) = 0;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	unknown(3,9) = 0;
	unknown(3,10) = 0;
	unknown(3,11) = 0;
	unknown(3,12) = 0;
	unknown(3,13) = 0;
	Jg_rot = unknown;
else
	Jg_rot=NaN(3,13);
end