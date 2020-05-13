% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh1m2TE_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2TE_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_jacobiR_rot_sym_varpar: pkin has to be [22x1] (double)');
JR_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:38
	% EndTime: 2020-05-01 20:48:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:38
	% EndTime: 2020-05-01 20:48:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0; t9, 0, 0, 0; 0, 0, 0, 0; -t9, 0, 0, 0; -t8, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:38
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t12 = cos(qJ(1));
	t9 = sin(qJ(2));
	t15 = t12 * t9;
	t10 = sin(qJ(1));
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t8 = t10 * t9;
	t1 = [t8, -t13, 0, 0; -t15, -t14, 0, 0; 0, -t9, 0, 0; t14, t15, 0, 0; -t13, t8, 0, 0; 0, -t11, 0, 0; t12, 0, 0, 0; t10, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->14), mult. (44->8), div. (0->0), fcn. (78->6), ass. (0->13)
	t27 = sin(qJ(3));
	t28 = sin(qJ(2));
	t30 = cos(qJ(3));
	t31 = cos(qJ(2));
	t25 = -t28 * t27 + t30 * t31;
	t29 = sin(qJ(1));
	t34 = t29 * t25;
	t26 = t31 * t27 + t30 * t28;
	t23 = t29 * t26;
	t32 = cos(qJ(1));
	t33 = t32 * t25;
	t24 = t32 * t26;
	t1 = [-t34, -t24, -t24, 0; t33, -t23, -t23, 0; 0, t25, t25, 0; t23, -t33, -t33, 0; -t24, -t34, -t34, 0; 0, -t26, -t26, 0; t32, 0, 0, 0; t29, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->7), mult. (28->12), div. (0->0), fcn. (46->8), ass. (0->14)
	t39 = cos(pkin(18));
	t38 = cos(qJ(1));
	t37 = sin(pkin(18));
	t36 = sin(qJ(1));
	t35 = cos(pkin(20));
	t34 = sin(pkin(20));
	t33 = pkin(22) + pkin(21);
	t32 = cos(t33);
	t31 = sin(t33);
	t30 = t37 * t34 + t39 * t35;
	t29 = -t39 * t34 + t37 * t35;
	t28 = -t29 * t32 + t31 * t30;
	t27 = t29 * t31 + t32 * t30;
	t1 = [t36 * t27, 0, 0, 0; -t38 * t27, 0, 0, 0; 0, 0, 0, 0; t36 * t28, 0, 0, 0; -t38 * t28, 0, 0, 0; 0, 0, 0, 0; t38, 0, 0, 0; t36, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (71->14), mult. (116->22), div. (0->0), fcn. (172->10), ass. (0->22)
	t101 = cos(pkin(18));
	t94 = sin(pkin(20));
	t95 = cos(pkin(20));
	t98 = sin(pkin(18));
	t89 = -t101 * t94 + t98 * t95;
	t90 = t101 * t95 + t98 * t94;
	t93 = pkin(22) + pkin(21);
	t91 = sin(t93);
	t92 = cos(t93);
	t104 = t89 * t91 + t90 * t92;
	t97 = sin(qJ(1));
	t103 = t104 * t97;
	t100 = cos(qJ(1));
	t102 = t100 * t104;
	t99 = cos(qJ(4));
	t96 = sin(qJ(4));
	t88 = -t89 * t92 + t91 * t90;
	t87 = t99 * t102 - t97 * t96;
	t86 = t96 * t102 + t97 * t99;
	t85 = t100 * t96 + t99 * t103;
	t84 = -t100 * t99 + t96 * t103;
	t1 = [t85, 0, 0, t86; -t87, 0, 0, t84; 0, 0, 0, -t88 * t96; -t84, 0, 0, t87; t86, 0, 0, t85; 0, 0, 0, -t88 * t99; -t97 * t88, 0, 0, 0; t100 * t88, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (32->7), mult. (68->12), div. (0->0), fcn. (110->8), ass. (0->17)
	t34 = sin(pkin(18));
	t37 = cos(pkin(17));
	t40 = sin(pkin(17));
	t41 = cos(pkin(18));
	t29 = t34 * t37 - t41 * t40;
	t30 = t40 * t34 + t37 * t41;
	t32 = sin(qJ(2));
	t35 = cos(qJ(2));
	t25 = -t29 * t35 + t30 * t32;
	t33 = sin(qJ(1));
	t45 = t33 * t25;
	t36 = cos(qJ(1));
	t44 = t36 * t25;
	t38 = -t32 * t29 - t30 * t35;
	t43 = t33 * t38;
	t42 = t36 * t38;
	t1 = [t45, t42, 0, 0; -t44, t43, 0, 0; 0, -t25, 0, 0; -t43, t44, 0, 0; t42, t45, 0, 0; 0, t38, 0, 0; t36, 0, 0, 0; t33, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->4), mult. (12->8), div. (0->0), fcn. (22->6), ass. (0->9)
	t22 = cos(pkin(18));
	t21 = cos(qJ(1));
	t20 = sin(pkin(18));
	t19 = sin(qJ(1));
	t18 = cos(pkin(22));
	t17 = sin(pkin(22));
	t16 = t22 * t17 - t20 * t18;
	t15 = t20 * t17 + t18 * t22;
	t1 = [t19 * t15, 0, 0, 0; -t21 * t15, 0, 0, 0; 0, 0, 0, 0; t19 * t16, 0, 0, 0; -t21 * t16, 0, 0, 0; 0, 0, 0, 0; t21, 0, 0, 0; t19, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (55->17), mult. (108->22), div. (0->0), fcn. (174->8), ass. (0->19)
	t56 = sin(pkin(19));
	t57 = cos(pkin(19));
	t59 = sin(qJ(2));
	t62 = cos(qJ(2));
	t52 = t56 * t59 - t57 * t62;
	t54 = t56 * t62 + t57 * t59;
	t58 = sin(qJ(3));
	t61 = cos(qJ(3));
	t43 = t52 * t61 + t58 * t54;
	t63 = cos(qJ(1));
	t73 = t63 * t43;
	t60 = sin(qJ(1));
	t67 = -t52 * t58 + t54 * t61;
	t68 = t60 * t67;
	t72 = t63 * t67;
	t51 = t61 * t56 + t58 * t57;
	t53 = t58 * t56 - t61 * t57;
	t64 = t60 * (t59 * t51 + t62 * t53);
	t1 = [t60 * t43, -t72, -t72, 0; -t73, -t68, -t68, 0; 0, -t43, -t43, 0; t68, t73, t73, 0; t63 * (-t51 * t62 + t59 * t53), t64, t64, 0; 0, -t67, -t67, 0; t63, 0, 0, 0; t60, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiR_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t26 = sin(qJ(1));
	t27 = cos(qJ(2));
	t31 = t26 * t27;
	t25 = sin(qJ(2));
	t28 = cos(qJ(1));
	t30 = t28 * t25;
	t29 = t28 * t27;
	t24 = t26 * t25;
	t1 = [t24, -t29, 0, 0; -t30, -t31, 0, 0; 0, -t25, 0, 0; t31, t30, 0, 0; -t29, t24, 0, 0; 0, -t27, 0, 0; t28, 0, 0, 0; t26, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiR_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t26 = qJ(3) + qJ(2);
	t24 = sin(t26);
	t27 = sin(qJ(1));
	t32 = t27 * t24;
	t25 = cos(t26);
	t31 = t27 * t25;
	t28 = cos(qJ(1));
	t30 = t28 * t24;
	t29 = t28 * t25;
	t1 = [-t31, -t30, -t30, 0; t29, -t32, -t32, 0; 0, t25, t25, 0; t32, -t29, -t29, 0; -t30, -t31, -t31, 0; 0, -t24, -t24, 0; t28, 0, 0, 0; t27, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
end