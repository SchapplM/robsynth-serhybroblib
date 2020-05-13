% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh4m1OL
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh4m1OL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
% 
% Output:
% Ja_rot [3x8]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:04
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh4m1OL_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1OL_jacobia_rot_sym_varpar: qJ has to be [8x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1OL_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1OL_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.01s
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
	Ja_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (6->3), div. (5->2), fcn. (6->2), ass. (0->30)
	unknown=NaN(3,8);
	t1 = sin(qJ(1));
	t2 = t1 ^ 2;
	t3 = cos(qJ(1));
	t4 = t3 ^ 2;
	t6 = t2 / t4;
	t8 = 0.1e1 / (0.1e1 + t6);
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
	unknown(3,1) = (t6 * t8 + t8);
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(3,6) = 0;
	unknown(3,7) = 0;
	unknown(3,8) = 0;
	Ja_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:09
	% EndTime: 2020-04-11 23:04:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->24)
	unknown=NaN(3,8);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	Ja_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (60->16), mult. (145->53), div. (49->9), fcn. (228->7), ass. (0->50)
	unknown=NaN(3,8);
	t1 = cos(qJ(1));
	t2 = cos(qJ(2));
	t3 = t1 * t2;
	t4 = sin(qJ(2));
	t5 = 0.1e1 / t4;
	t6 = sin(qJ(1));
	t7 = t6 ^ 2;
	t8 = t2 ^ 2;
	t10 = t4 ^ 2;
	t11 = 0.1e1 / t10;
	t14 = 0.1e1 / (t7 * t8 * t11 + 0.1e1);
	t21 = t8 * t6 * t11 * t14 + t6 * t14;
	t22 = t6 * t2;
	t23 = atan2(-t22, t4);
	t24 = cos(t23);
	t26 = sin(t23);
	t27 = t26 * t6;
	t29 = -t27 * t2 + t24 * t4;
	t31 = t1 ^ 2;
	t33 = t29 ^ 2;
	t34 = 0.1e1 / t33;
	t37 = 0.1e1 / (t31 * t8 * t34 + 0.1e1);
	t38 = 0.1e1 / t29 * t37;
	t52 = t2 * t34 * t37;
	t67 = 0.1e1 / t31;
	t71 = 0.1e1 / (t7 * t67 * t11 + 0.1e1);
	unknown(1,1) = -t3 * t5 * t14;
	unknown(1,2) = t21;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = -t22 * t38 - (t1 * t8 * t5 * t14 * t24 * t6 - t26 * t1 * t2 + t3 * t14 * t26) * t1 * t52;
	unknown(2,2) = -t1 * t4 * t38 - (-t21 * t24 * t22 - t21 * t26 * t4 + t24 * t2 + t27 * t4) * t1 * t52;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = t5 * t7 * t67 * t71 + t5 * t71;
	unknown(3,2) = -0.1e1 / t1 * t2 * t6 * t11 * t71;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(3,6) = 0.0e0;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	Ja_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->24)
	unknown=NaN(3,8);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	Ja_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:10
	% EndTime: 2020-04-11 23:04:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->24)
	unknown=NaN(3,8);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	Ja_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:11
	% EndTime: 2020-04-11 23:04:11
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (994->26), mult. (1174->77), div. (115->9), fcn. (1719->11), ass. (0->70)
	unknown=NaN(3,8);
	t1 = cos(qJ(1));
	t2 = cos(qJ(2));
	t3 = t1 * t2;
	t4 = qJ(4) + qJ(5);
	t5 = cos(t4);
	t7 = sin(qJ(2));
	t8 = t1 * t7;
	t9 = sin(t4);
	t11 = -t3 * t5 + t8 * t9;
	t14 = t2 * t9 + t7 * t5;
	t15 = 0.1e1 / t14;
	t16 = t11 * t15;
	t17 = sin(qJ(1));
	t18 = t17 * t7;
	t20 = t17 * t2;
	t22 = t18 * t9 - t20 * t5;
	t23 = t22 ^ 2;
	t24 = t14 ^ 2;
	t25 = 0.1e1 / t24;
	t28 = 0.1e1 / (t23 * t25 + 0.1e1);
	t32 = t18 * t5 + t20 * t9;
	t37 = t2 * t5 - t7 * t9;
	t41 = -t37 * t22 * t25 * t28 + t32 * t15 * t28;
	t42 = atan2(t22, t14);
	t43 = cos(t42);
	t45 = sin(t42);
	t47 = t43 * t14 + t45 * t22;
	t48 = 0.1e1 / t47;
	t50 = t11 ^ 2;
	t51 = t47 ^ 2;
	t52 = 0.1e1 / t51;
	t55 = 0.1e1 / (t50 * t52 + 0.1e1);
	t65 = t52 * t55;
	t70 = -t3 * t9 - t8 * t5;
	t82 = t70 * t48 * t55 + (-t41 * t45 * t14 + t41 * t43 * t22 + t45 * t32 + t43 * t37) * t11 * t65;
	t83 = sin(qJ(6));
	t85 = cos(qJ(6));
	t90 = -t17 * t83 - t70 * t85;
	t91 = 0.1e1 / t90;
	t95 = t17 * t85 - t70 * t83;
	t96 = t95 ^ 2;
	t97 = t90 ^ 2;
	t98 = 0.1e1 / t97;
	t101 = 0.1e1 / (t96 * t98 + 0.1e1);
	t107 = t98 * t101;
	t117 = t11 * t85 * t95 * t98 * t101 - t11 * t83 * t91 * t101;
	unknown(1,1) = t16 * t28;
	unknown(1,2) = t41;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = t41;
	unknown(1,5) = t41;
	unknown(1,6) = 0.0e0;
	unknown(1,7) = 0.0e0;
	unknown(1,8) = 0.0e0;
	unknown(2,1) = t22 * t48 * t55 + (t16 * t28 * t43 * t22 - t11 * t28 * t45 + t45 * t11) * t11 * t65;
	unknown(2,2) = t82;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = t82;
	unknown(2,5) = t82;
	unknown(2,6) = 0.0e0;
	unknown(2,7) = 0.0e0;
	unknown(2,8) = 0.0e0;
	unknown(3,1) = (t1 * t85 - t32 * t83) * t91 * t101 - (-t1 * t83 - t32 * t85) * t95 * t107;
	unknown(3,2) = t117;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = t117;
	unknown(3,5) = t117;
	unknown(3,6) = t95 ^ 2 * t107 + t101;
	unknown(3,7) = 0.0e0;
	unknown(3,8) = 0.0e0;
	Ja_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 23:04:12
	% EndTime: 2020-04-11 23:04:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->24)
	unknown=NaN(3,8);
	unknown(1,1) = NaN;
	unknown(1,2) = NaN;
	unknown(1,3) = NaN;
	unknown(1,4) = NaN;
	unknown(1,5) = NaN;
	unknown(1,6) = NaN;
	unknown(1,7) = NaN;
	unknown(1,8) = NaN;
	unknown(2,1) = NaN;
	unknown(2,2) = NaN;
	unknown(2,3) = NaN;
	unknown(2,4) = NaN;
	unknown(2,5) = NaN;
	unknown(2,6) = NaN;
	unknown(2,7) = NaN;
	unknown(2,8) = NaN;
	unknown(3,1) = NaN;
	unknown(3,2) = NaN;
	unknown(3,3) = NaN;
	unknown(3,4) = NaN;
	unknown(3,5) = NaN;
	unknown(3,6) = NaN;
	unknown(3,7) = NaN;
	unknown(3,8) = NaN;
	Ja_rot = unknown;
else
	Ja_rot=NaN(3,8);
end