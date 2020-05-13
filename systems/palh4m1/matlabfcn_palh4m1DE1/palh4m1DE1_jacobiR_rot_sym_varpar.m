% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh4m1DE1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AD,CB,CE,EP,HC,OT,TA,TD]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 22:26
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh4m1DE1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1DE1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1DE1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1DE1_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->45)
	unknown=NaN(9,5);
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = 0;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = 0;
	unknown(3,1) = 0;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	unknown(4,1) = 0;
	unknown(4,2) = 0;
	unknown(4,3) = 0;
	unknown(4,4) = 0;
	unknown(4,5) = 0;
	unknown(5,1) = 0;
	unknown(5,2) = 0;
	unknown(5,3) = 0;
	unknown(5,4) = 0;
	unknown(5,5) = 0;
	unknown(6,1) = 0;
	unknown(6,2) = 0;
	unknown(6,3) = 0;
	unknown(6,4) = 0;
	unknown(6,5) = 0;
	unknown(7,1) = 0;
	unknown(7,2) = 0;
	unknown(7,3) = 0;
	unknown(7,4) = 0;
	unknown(7,5) = 0;
	unknown(8,1) = 0;
	unknown(8,2) = 0;
	unknown(8,3) = 0;
	unknown(8,4) = 0;
	unknown(8,5) = 0;
	unknown(9,1) = 0;
	unknown(9,2) = 0;
	unknown(9,3) = 0;
	unknown(9,4) = 0;
	unknown(9,5) = 0;
	JR_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->47)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = -t1;
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = 0.0e0;
	unknown(2,1) = t2;
	unknown(2,2) = 0.0e0;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = 0.0e0;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = 0.0e0;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = 0.0e0;
	unknown(4,1) = -t2;
	unknown(4,2) = 0.0e0;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(5,1) = -t1;
	unknown(5,2) = 0.0e0;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = 0.0e0;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(7,1) = 0.0e0;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(8,1) = 0.0e0;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:50
	% EndTime: 2020-04-11 22:24:50
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (5404->62), mult. (5732->131), div. (264->9), fcn. (1576->8), ass. (0->95)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = cos(qJ(5));
	t3 = sin(qJ(5));
	t4 = pkin(1) * t3;
	t6 = 0.2e1 * pkin(2) * t4;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t15 * t11);
	t19 = pkin(1) * t17 * t2;
	t20 = t4 - pkin(2);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t27 = -t25 * t20 - t19;
	t29 = 0.1e1 / t22;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t34 = t2 * pkin(1);
	t35 = t25 * t34;
	t36 = -t17 * t20 + t35;
	t37 = t36 ^ 2;
	t38 = 0.1e1 / t23;
	t40 = t30 ^ 2;
	t41 = 0.1e1 / t40;
	t43 = t27 ^ 2;
	t47 = sqrt(t41 * t38 * t37 + t41 * t38 * t43);
	t49 = 0.1e1 / t47 * t31 * t29;
	t51 = cos(qJ(1));
	t52 = 0.1e1 / t17;
	t53 = -t52 * t20;
	t58 = -0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t68 = 0.1e1 / t27;
	t71 = 0.1e1 / t43;
	t74 = 0.1e1 / (t71 * t37 + 0.1e1);
	t75 = t74 * t30 * t22;
	t77 = t52 * t2;
	t91 = t74 * t71 * t30;
	t93 = t75 * t68 * (t31 * t29 * (t58 * t53 / 0.2e1 + 0.2e1 * t22 * t34) - t31 * t38 * t36) - t91 * t22 * t36 * (t31 * t29 * (-t58 * pkin(1) * t77 / 0.2e1 - 0.2e1 * t22 * t20) - t31 * t38 * t27);
	t94 = t93 * t51;
	t100 = pkin(1) * pkin(2);
	t102 = t15 * pkin(2) * t34 + t100 * t2 * t11;
	t106 = t2 ^ 2;
	t115 = pkin(2) * t34;
	t140 = t75 * t68 * (t31 * t29 * (-0.2e1 * pkin(2) * t106 * t7 + t102 * t53 - t25 * t4 - t19) + 0.2e1 * t115 * t41 * t29 * t36) - t91 * t22 * t36 * (t31 * t29 * (-t102 * pkin(1) * t77 + pkin(1) * t17 * t3 + 0.2e1 * t100 * t2 * t20 - t35) + 0.2e1 * t115 * t41 * t29 * t27);
	t141 = t140 * t51;
	t146 = t93 * t1;
	t149 = t140 * t1;
	unknown(1,1) = -t49 * t27 * t1;
	unknown(1,2) = -t49 * t36 * t94;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = -t49 * t36 * t141;
	unknown(2,1) = t49 * t27 * t51;
	unknown(2,2) = -t49 * t36 * t146;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = -t49 * t36 * t149;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t49 * t27 * t93;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = t49 * t27 * t140;
	unknown(4,1) = t49 * t36 * t1;
	unknown(4,2) = -t49 * t27 * t94;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = -t49 * t27 * t141;
	unknown(5,1) = -t49 * t36 * t51;
	unknown(5,2) = -t49 * t27 * t146;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = -t49 * t27 * t149;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = -t49 * t36 * t93;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = -t49 * t36 * t140;
	unknown(7,1) = t51;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(8,1) = t1;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:51
	% EndTime: 2020-04-11 22:24:51
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (5400->57), mult. (5732->131), div. (264->9), fcn. (1576->8), ass. (0->95)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = pkin(1) * t2;
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * pkin(2) * t3;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t15 * t11);
	t19 = cos(qJ(5));
	t20 = t19 * pkin(1);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t25 * t20;
	t27 = t17 * t4 + t26;
	t29 = 0.1e1 / t22;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t33 = t27 ^ 2;
	t34 = 0.1e1 / t23;
	t36 = t30 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t47 = sqrt(t37 * t34 * t33 + t37 * t34 * t43);
	t49 = 0.1e1 / t47 * t31 * t29;
	t51 = cos(qJ(1));
	t52 = 0.1e1 / t17;
	t53 = t52 * t4;
	t58 = -0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t68 = 0.1e1 / t42;
	t71 = 0.1e1 / t43;
	t74 = 0.1e1 / (t71 * t33 + 0.1e1);
	t75 = t74 * t30 * t22;
	t77 = t52 * t19;
	t91 = t74 * t71 * t30;
	t93 = t75 * t68 * (t31 * t29 * (t58 * t53 / 0.2e1 + 0.2e1 * t22 * t20) - t31 * t34 * t27) - t91 * t22 * t27 * (t31 * t29 * (-t58 * pkin(1) * t77 / 0.2e1 + 0.2e1 * t22 * t4) - t31 * t34 * t42);
	t94 = t93 * t51;
	t100 = pkin(1) * pkin(2);
	t102 = t15 * pkin(2) * t20 + t100 * t19 * t11;
	t106 = t19 ^ 2;
	t115 = pkin(2) * t20;
	t140 = t75 * t68 * (t31 * t29 * (-0.2e1 * pkin(2) * t106 * t7 + t102 * t53 - t25 * t3 - t40) + 0.2e1 * t115 * t37 * t29 * t27) - t91 * t22 * t27 * (t31 * t29 * (-t102 * pkin(1) * t77 + pkin(1) * t17 * t2 - 0.2e1 * t100 * t19 * t4 - t26) + 0.2e1 * t115 * t37 * t29 * t42);
	t141 = t140 * t51;
	t146 = t93 * t1;
	t149 = t140 * t1;
	unknown(1,1) = -t49 * t27 * t1;
	unknown(1,2) = t49 * t42 * t94;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t49 * t42 * t141;
	unknown(2,1) = t49 * t27 * t51;
	unknown(2,2) = t49 * t42 * t146;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t49 * t42 * t149;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t49 * t27 * t93;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = t49 * t27 * t140;
	unknown(4,1) = -t51;
	unknown(4,2) = 0.0e0;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = 0.0e0;
	unknown(5,1) = -t1;
	unknown(5,2) = 0.0e0;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = 0.0e0;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = 0.0e0;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = 0.0e0;
	unknown(7,1) = -t49 * t42 * t1;
	unknown(7,2) = -t49 * t27 * t94;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = -t49 * t27 * t141;
	unknown(8,1) = t49 * t42 * t51;
	unknown(8,2) = -t49 * t27 * t146;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = -t49 * t27 * t149;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = t49 * t42 * t93;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = t49 * t42 * t140;
	JR_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:51
	% EndTime: 2020-04-11 22:24:51
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (17716->99), mult. (18468->240), div. (1192->15), fcn. (4738->10), ass. (0->148)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = pkin(1) * t2;
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * pkin(2) * t3;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t16 = t15 * t11;
	t17 = sqrt(-t16);
	t19 = cos(qJ(5));
	t20 = t19 * pkin(1);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t25 * t20;
	t27 = t17 * t4 + t26;
	t29 = 0.1e1 / t23;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t32 = t31 * t29;
	t33 = t32 * t27 * t1;
	t34 = t27 ^ 2;
	t36 = t30 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t47 = sqrt(t37 * t29 * t34 + t37 * t29 * t43);
	t48 = 0.1e1 / t47;
	t49 = 0.1e1 / pkin(3);
	t51 = 0.1e1 / t24;
	t54 = t6 - t7 - t21 + t23 + t24;
	t55 = t54 ^ 2;
	t59 = sqrt(-t16 * t51 * t29 + t51 * t29 * t55);
	t60 = 0.1e1 / t59;
	t61 = t60 * t17;
	t62 = t61 * t49 * t48;
	t65 = t32 * t42 * t1;
	t67 = t60 * t49;
	t68 = t67 * t54 * t48;
	t71 = cos(qJ(1));
	t72 = 0.1e1 / t17;
	t73 = t72 * t4;
	t78 = -0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t83 = 0.1e1 / t22;
	t86 = t29 * t27;
	t89 = 0.1e1 / t42;
	t92 = 0.1e1 / t43;
	t95 = 0.1e1 / (t92 * t34 + 0.1e1);
	t96 = t95 * t30 * t22;
	t98 = t72 * t19;
	t106 = t29 * t42;
	t112 = t95 * t92 * t30;
	t114 = t96 * t89 * (t31 * t83 * (t78 * t73 / 0.2e1 + 0.2e1 * t22 * t20) - t31 * t86) - t112 * t22 * t27 * (t31 * t83 * (-t78 * pkin(1) * t98 / 0.2e1 + 0.2e1 * t22 * t4) - t31 * t106);
	t115 = t114 * t71;
	t116 = t106 * t115;
	t117 = t48 * t31;
	t119 = t60 * t17 * t49;
	t120 = t119 * t117;
	t123 = t32 * t27 * t71;
	t131 = 0.1e1 / t54;
	t134 = 0.1e1 / t55;
	t137 = 0.1e1 / (-t134 * t16 + 0.1e1);
	t148 = t137 * t134 * t17;
	t150 = t137 * pkin(3) * t22 * t131 * (-t17 * t49 * t29 + t78 * t72 * t49 * t83 / 0.2e1) - t148 * pkin(3) * t22 * (0.2e1 * t49 * t83 * t22 - t49 * t29 * t54);
	t151 = t150 * t48;
	t153 = t60 * t49 * t54;
	t154 = t153 * t151;
	t156 = t86 * t115;
	t157 = t153 * t117;
	t160 = t32 * t42 * t71;
	t161 = t119 * t151;
	t167 = pkin(1) * pkin(2);
	t169 = t15 * pkin(2) * t20 + t167 * t19 * t11;
	t173 = t19 ^ 2;
	t182 = pkin(2) * t20;
	t207 = t96 * t89 * (t31 * t83 * (-0.2e1 * pkin(2) * t173 * t7 + t169 * t73 - t25 * t3 - t40) + 0.2e1 * t182 * t37 * t83 * t27) - t112 * t22 * t27 * (t31 * t83 * (-t169 * pkin(1) * t98 + pkin(1) * t17 * t2 - 0.2e1 * t167 * t19 * t4 - t26) + 0.2e1 * t182 * t37 * t83 * t42);
	t208 = t207 * t71;
	t209 = t106 * t208;
	t217 = t137 * t131 * t169 * t72 - 0.2e1 * t148 * t182;
	t218 = t217 * t48;
	t219 = t153 * t218;
	t221 = t86 * t208;
	t223 = t119 * t218;
	t229 = t114 * t1;
	t230 = t106 * t229;
	t233 = t86 * t229;
	t237 = t207 * t1;
	t238 = t106 * t237;
	t241 = t86 * t237;
	t246 = t32 * t27 * t114;
	t248 = t117 * t106;
	t250 = t67 * t54 * t150;
	t253 = t32 * t42 * t114;
	t255 = t117 * t86;
	t257 = t61 * t49 * t150;
	t261 = t32 * t27 * t207;
	t264 = t67 * t54 * t217;
	t267 = t32 * t42 * t207;
	t270 = t61 * t49 * t217;
	unknown(1,1) = t62 * t33 - t68 * t65;
	unknown(1,2) = -t120 * t116 - t154 * t123 - t157 * t156 - t161 * t160;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = -t120 * t209 - t219 * t123 - t157 * t221 - t223 * t160;
	unknown(2,1) = -t62 * t123 + t68 * t160;
	unknown(2,2) = -t120 * t230 - t154 * t33 - t157 * t233 - t161 * t65;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = -t120 * t238 - t157 * t241 - t219 * t33 - t223 * t65;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t62 * t246 + t250 * t248 + t68 * t253 - t257 * t255;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = t264 * t248 - t270 * t255 - t62 * t261 + t68 * t267;
	unknown(4,1) = t68 * t33 + t62 * t65;
	unknown(4,2) = -t157 * t116 + t120 * t156 + t161 * t123 - t154 * t160;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = t120 * t221 + t223 * t123 - t157 * t209 - t219 * t160;
	unknown(5,1) = -t68 * t123 - t62 * t160;
	unknown(5,2) = t120 * t233 - t154 * t65 - t157 * t230 + t161 * t33;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = t120 * t241 - t157 * t238 - t219 * t65 + t223 * t33;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = -t68 * t246 - t257 * t248 - t250 * t255 - t62 * t253;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = -t270 * t248 - t264 * t255 - t68 * t261 - t62 * t267;
	unknown(7,1) = t71;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(8,1) = t1;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:51
	% EndTime: 2020-04-11 22:24:52
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (38274->122), mult. (40080->284), div. (2648->15), fcn. (10338->12), ass. (0->172)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = pkin(1) * t2;
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * pkin(2) * t3;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t16 = t15 * t11;
	t17 = sqrt(-t16);
	t19 = cos(qJ(5));
	t20 = t19 * pkin(1);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t25 * t20;
	t27 = t17 * t4 + t26;
	t29 = 0.1e1 / t23;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t32 = t31 * t29;
	t33 = t32 * t27 * t1;
	t34 = t27 ^ 2;
	t36 = t30 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t47 = sqrt(t37 * t29 * t34 + t37 * t29 * t43);
	t48 = 0.1e1 / t47;
	t49 = 0.1e1 / pkin(3);
	t51 = 0.1e1 / t24;
	t54 = t6 - t7 - t21 + t23 + t24;
	t55 = t54 ^ 2;
	t59 = sqrt(-t16 * t51 * t29 + t51 * t29 * t55);
	t60 = 0.1e1 / t59;
	t61 = t60 * t17;
	t62 = t61 * t49 * t48;
	t65 = t32 * t42 * t1;
	t67 = t60 * t49;
	t68 = t67 * t54 * t48;
	t70 = t62 * t33 - t68 * t65;
	t71 = cos(qJ(3));
	t75 = t68 * t33 + t62 * t65;
	t76 = sin(qJ(3));
	t79 = cos(qJ(1));
	t80 = 0.1e1 / t17;
	t81 = t80 * t4;
	t86 = -0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t91 = 0.1e1 / t22;
	t94 = t29 * t27;
	t95 = t31 * t94;
	t97 = 0.1e1 / t42;
	t100 = 0.1e1 / t43;
	t103 = 0.1e1 / (t100 * t34 + 0.1e1);
	t104 = t103 * t30 * t22;
	t106 = t80 * t19;
	t114 = t29 * t42;
	t115 = t31 * t114;
	t120 = t103 * t100 * t30;
	t122 = t104 * t97 * (t31 * t91 * (t86 * t81 / 0.2e1 + 0.2e1 * t22 * t20) - t95) - t120 * t22 * t27 * (t31 * t91 * (-t86 * pkin(1) * t106 / 0.2e1 + 0.2e1 * t22 * t4) - t115);
	t123 = t122 * t79;
	t124 = t114 * t123;
	t125 = t48 * t31;
	t127 = t60 * t17 * t49;
	t128 = t127 * t125;
	t131 = t32 * t27 * t79;
	t139 = 0.1e1 / t54;
	t142 = 0.1e1 / t55;
	t145 = 0.1e1 / (-t142 * t16 + 0.1e1);
	t156 = t145 * t142 * t17;
	t158 = t145 * pkin(3) * t22 * t139 * (-t17 * t49 * t29 + t86 * t80 * t49 * t91 / 0.2e1) - t156 * pkin(3) * t22 * (0.2e1 * t49 * t91 * t22 - t49 * t29 * t54);
	t159 = t158 * t48;
	t161 = t60 * t49 * t54;
	t162 = t161 * t159;
	t164 = t94 * t123;
	t165 = t161 * t125;
	t168 = t32 * t42 * t79;
	t169 = t127 * t159;
	t171 = -t128 * t124 - t162 * t131 - t165 * t164 - t169 * t168;
	t177 = -t165 * t124 + t128 * t164 + t169 * t131 - t162 * t168;
	t182 = -t62 * t131 + t68 * t168;
	t186 = -t68 * t131 - t62 * t168;
	t188 = -t76 * t182 + t71 * t186;
	t192 = pkin(1) * pkin(2);
	t194 = t15 * pkin(2) * t20 + t192 * t19 * t11;
	t198 = t19 ^ 2;
	t207 = pkin(2) * t20;
	t232 = t104 * t97 * (t31 * t91 * (-0.2e1 * pkin(2) * t198 * t7 + t194 * t81 - t25 * t3 - t40) + 0.2e1 * t207 * t37 * t91 * t27) - t120 * t22 * t27 * (t31 * t91 * (-t194 * pkin(1) * t106 + pkin(1) * t17 * t2 - 0.2e1 * t192 * t19 * t4 - t26) + 0.2e1 * t207 * t37 * t91 * t42);
	t233 = t232 * t79;
	t234 = t114 * t233;
	t242 = t145 * t139 * t194 * t80 - 0.2e1 * t156 * t207;
	t243 = t242 * t48;
	t244 = t161 * t243;
	t246 = t94 * t233;
	t248 = t127 * t243;
	t250 = -t128 * t234 - t244 * t131 - t165 * t246 - t248 * t168;
	t256 = t128 * t246 + t248 * t131 - t165 * t234 - t244 * t168;
	t261 = t71 * t182 + t76 * t186;
	t262 = t122 * t1;
	t263 = t114 * t262;
	t266 = t94 * t262;
	t269 = -t128 * t263 - t162 * t33 - t165 * t266 - t169 * t65;
	t275 = t128 * t266 - t162 * t65 - t165 * t263 + t169 * t33;
	t281 = t232 * t1;
	t282 = t114 * t281;
	t285 = t94 * t281;
	t288 = -t128 * t282 - t165 * t285 - t244 * t33 - t248 * t65;
	t294 = t128 * t285 - t165 * t282 - t244 * t65 + t248 * t33;
	t298 = t32 * t27 * t122;
	t300 = t125 * t114;
	t302 = t67 * t54 * t158;
	t305 = t32 * t42 * t122;
	t307 = t125 * t94;
	t309 = t61 * t49 * t158;
	t311 = -t62 * t298 + t302 * t300 + t68 * t305 - t309 * t307;
	t317 = -t68 * t298 - t309 * t300 - t302 * t307 - t62 * t305;
	t322 = t62 * t115 + t68 * t95;
	t326 = t68 * t115 - t62 * t95;
	t330 = t32 * t27 * t232;
	t333 = t67 * t54 * t242;
	t336 = t32 * t42 * t232;
	t339 = t61 * t49 * t242;
	t341 = t333 * t300 - t339 * t307 - t62 * t330 + t68 * t336;
	t347 = -t339 * t300 - t333 * t307 - t68 * t330 - t62 * t336;
	unknown(1,1) = t71 * t70 + t76 * t75;
	unknown(1,2) = t71 * t171 + t76 * t177;
	unknown(1,3) = t188;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t71 * t250 + t76 * t256;
	unknown(2,1) = t261;
	unknown(2,2) = t71 * t269 + t76 * t275;
	unknown(2,3) = t76 * t70 - t71 * t75;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t71 * t288 + t76 * t294;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t71 * t311 + t76 * t317;
	unknown(3,3) = -t76 * t322 + t71 * t326;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = t71 * t341 + t76 * t347;
	unknown(4,1) = -t76 * t70 + t71 * t75;
	unknown(4,2) = -t76 * t171 + t71 * t177;
	unknown(4,3) = -t261;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = -t76 * t250 + t71 * t256;
	unknown(5,1) = t188;
	unknown(5,2) = -t76 * t269 + t71 * t275;
	unknown(5,3) = t71 * t70 + t76 * t75;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = -t76 * t288 + t71 * t294;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = -t76 * t311 + t71 * t317;
	unknown(6,3) = -t71 * t322 - t76 * t326;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = -t76 * t341 + t71 * t347;
	unknown(7,1) = t79;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(8,1) = t1;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	JR_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:55
	% EndTime: 2020-04-11 22:24:56
	% DurationCPUTime: 0.90s
	% Computational Cost: add. (60254->138), mult. (63268->314), div. (4236->15), fcn. (16380->14), ass. (0->189)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = pkin(1) * t2;
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * pkin(2) * t3;
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t16 = t15 * t11;
	t17 = sqrt(-t16);
	t19 = cos(qJ(5));
	t20 = t19 * pkin(1);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t25 * t20;
	t27 = t17 * t4 + t26;
	t29 = 0.1e1 / t23;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t32 = t31 * t29;
	t33 = t32 * t27 * t1;
	t34 = t27 ^ 2;
	t36 = t30 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t47 = sqrt(t37 * t29 * t34 + t37 * t29 * t43);
	t48 = 0.1e1 / t47;
	t49 = 0.1e1 / pkin(3);
	t51 = 0.1e1 / t24;
	t54 = t6 - t7 - t21 + t23 + t24;
	t55 = t54 ^ 2;
	t59 = sqrt(-t16 * t51 * t29 + t51 * t29 * t55);
	t60 = 0.1e1 / t59;
	t61 = t60 * t17;
	t62 = t61 * t49 * t48;
	t65 = t32 * t42 * t1;
	t67 = t60 * t49;
	t68 = t67 * t54 * t48;
	t70 = t62 * t33 - t68 * t65;
	t71 = sin(qJ(3));
	t75 = t68 * t33 + t62 * t65;
	t76 = cos(qJ(3));
	t78 = -t71 * t70 + t76 * t75;
	t79 = cos(qJ(4));
	t81 = cos(qJ(1));
	t82 = sin(qJ(4));
	t83 = t82 * t81;
	t85 = 0.1e1 / t17;
	t86 = t85 * t4;
	t91 = -0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t96 = 0.1e1 / t22;
	t99 = t29 * t27;
	t100 = t31 * t99;
	t102 = 0.1e1 / t42;
	t105 = 0.1e1 / t43;
	t108 = 0.1e1 / (t105 * t34 + 0.1e1);
	t109 = t108 * t30 * t22;
	t111 = t85 * t19;
	t119 = t29 * t42;
	t120 = t31 * t119;
	t125 = t108 * t105 * t30;
	t127 = t109 * t102 * (t31 * t96 * (t91 * t86 / 0.2e1 + 0.2e1 * t22 * t20) - t100) - t125 * t22 * t27 * (t31 * t96 * (-t91 * pkin(1) * t111 / 0.2e1 + 0.2e1 * t22 * t4) - t120);
	t128 = t127 * t81;
	t129 = t119 * t128;
	t130 = t48 * t31;
	t132 = t60 * t17 * t49;
	t133 = t132 * t130;
	t136 = t32 * t27 * t81;
	t144 = 0.1e1 / t54;
	t147 = 0.1e1 / t55;
	t150 = 0.1e1 / (-t147 * t16 + 0.1e1);
	t161 = t150 * t147 * t17;
	t163 = t150 * pkin(3) * t22 * t144 * (-t17 * t49 * t29 + t91 * t85 * t49 * t96 / 0.2e1) - t161 * pkin(3) * t22 * (0.2e1 * t49 * t96 * t22 - t49 * t29 * t54);
	t164 = t163 * t48;
	t166 = t60 * t49 * t54;
	t167 = t166 * t164;
	t169 = t99 * t128;
	t170 = t166 * t130;
	t173 = t32 * t42 * t81;
	t174 = t132 * t164;
	t176 = -t133 * t129 - t167 * t136 - t170 * t169 - t174 * t173;
	t182 = -t170 * t129 + t133 * t169 + t174 * t136 - t167 * t173;
	t184 = -t71 * t176 + t76 * t182;
	t188 = -t62 * t136 + t68 * t173;
	t192 = -t68 * t136 - t62 * t173;
	t194 = -t76 * t188 - t71 * t192;
	t198 = -t71 * t188 + t76 * t192;
	t201 = -t79 * t1 + t82 * t198;
	t205 = pkin(1) * pkin(2);
	t207 = t15 * pkin(2) * t20 + t205 * t19 * t11;
	t211 = t19 ^ 2;
	t220 = pkin(2) * t20;
	t245 = t109 * t102 * (t31 * t96 * (-0.2e1 * pkin(2) * t211 * t7 + t207 * t86 - t25 * t3 - t40) + 0.2e1 * t220 * t37 * t96 * t27) - t125 * t22 * t27 * (t31 * t96 * (-t207 * pkin(1) * t111 + pkin(1) * t17 * t2 - 0.2e1 * t205 * t19 * t4 - t26) + 0.2e1 * t220 * t37 * t96 * t42);
	t246 = t245 * t81;
	t247 = t119 * t246;
	t255 = t150 * t144 * t207 * t85 - 0.2e1 * t161 * t220;
	t256 = t255 * t48;
	t257 = t166 * t256;
	t259 = t99 * t246;
	t261 = t132 * t256;
	t263 = -t133 * t247 - t257 * t136 - t170 * t259 - t261 * t173;
	t269 = t133 * t259 + t261 * t136 - t170 * t247 - t257 * t173;
	t271 = -t71 * t263 + t76 * t269;
	t275 = -t82 * t1 - t79 * t198;
	t276 = t127 * t1;
	t277 = t119 * t276;
	t280 = t99 * t276;
	t283 = -t133 * t277 - t167 * t33 - t170 * t280 - t174 * t65;
	t289 = t133 * t280 - t167 * t65 - t170 * t277 + t174 * t33;
	t291 = -t71 * t283 + t76 * t289;
	t295 = t76 * t70 + t71 * t75;
	t299 = t71 * t70 - t76 * t75;
	t301 = t79 * t81;
	t303 = t245 * t1;
	t304 = t119 * t303;
	t307 = t99 * t303;
	t310 = -t133 * t304 - t170 * t307 - t257 * t33 - t261 * t65;
	t316 = t133 * t307 - t170 * t304 - t257 * t65 + t261 * t33;
	t318 = -t71 * t310 + t76 * t316;
	t321 = t32 * t27 * t127;
	t323 = t130 * t119;
	t325 = t67 * t54 * t163;
	t328 = t32 * t42 * t127;
	t330 = t130 * t99;
	t332 = t61 * t49 * t163;
	t334 = -t62 * t321 + t325 * t323 + t68 * t328 - t332 * t330;
	t340 = -t68 * t321 - t332 * t323 - t325 * t330 - t62 * t328;
	t342 = -t71 * t334 + t76 * t340;
	t346 = t68 * t100 + t62 * t120;
	t350 = -t62 * t100 + t68 * t120;
	t352 = -t76 * t346 - t71 * t350;
	t356 = -t71 * t346 + t76 * t350;
	t359 = t32 * t27 * t245;
	t362 = t67 * t54 * t255;
	t365 = t32 * t42 * t245;
	t368 = t61 * t49 * t255;
	t370 = t362 * t323 - t368 * t330 - t62 * t359 + t68 * t365;
	t376 = -t368 * t323 - t362 * t330 - t68 * t359 - t62 * t365;
	t378 = -t71 * t370 + t76 * t376;
	unknown(1,1) = -t79 * t78 - t83;
	unknown(1,2) = -t79 * t184;
	unknown(1,3) = -t79 * t194;
	unknown(1,4) = t201;
	unknown(1,5) = -t79 * t271;
	unknown(2,1) = t275;
	unknown(2,2) = -t79 * t291;
	unknown(2,3) = -t79 * t295;
	unknown(2,4) = t82 * t299 + t301;
	unknown(2,5) = -t79 * t318;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t79 * t342;
	unknown(3,3) = -t79 * t352;
	unknown(3,4) = t82 * t356;
	unknown(3,5) = -t79 * t378;
	unknown(4,1) = t82 * t78 - t301;
	unknown(4,2) = t82 * t184;
	unknown(4,3) = t82 * t194;
	unknown(4,4) = -t275;
	unknown(4,5) = t82 * t271;
	unknown(5,1) = t201;
	unknown(5,2) = t82 * t291;
	unknown(5,3) = t82 * t295;
	unknown(5,4) = t79 * t299 - t83;
	unknown(5,5) = t82 * t318;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t82 * t342;
	unknown(6,3) = t82 * t352;
	unknown(6,4) = t79 * t356;
	unknown(6,5) = t82 * t378;
	unknown(7,1) = t76 * t70 + t71 * t75;
	unknown(7,2) = t76 * t176 + t71 * t182;
	unknown(7,3) = t198;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = t76 * t263 + t71 * t269;
	unknown(8,1) = -t194;
	unknown(8,2) = t76 * t283 + t71 * t289;
	unknown(8,3) = t299;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = t76 * t310 + t71 * t316;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = t76 * t334 + t71 * t340;
	unknown(9,3) = t356;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = t76 * t370 + t71 * t376;
	JR_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->53)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = t1 * t2;
	t4 = cos(qJ(1));
	t5 = cos(qJ(5));
	t6 = t4 * t5;
	t7 = t4 * t2;
	t8 = t1 * t5;
	unknown(1,1) = t3;
	unknown(1,2) = 0.0e0;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = -t6;
	unknown(2,1) = -t7;
	unknown(2,2) = 0.0e0;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = -t8;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = 0.0e0;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = -t2;
	unknown(4,1) = t8;
	unknown(4,2) = 0.0e0;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = t7;
	unknown(5,1) = -t6;
	unknown(5,2) = 0.0e0;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = t3;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = 0.0e0;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = -t5;
	unknown(7,1) = t4;
	unknown(7,2) = 0.0e0;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = 0.0e0;
	unknown(8,1) = t1;
	unknown(8,2) = 0.0e0;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = 0.0e0;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = 0.0e0;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = 0.0e0;
	JR_rot = unknown;
else
	JR_rot=NaN(9,5);
end