% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh4m1DE2
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
% Datum: 2020-04-11 22:54
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh4m1DE2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1DE2_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1DE2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1DE2_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
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
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
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
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:33
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
	% StartTime: 2020-04-11 22:53:33
	% EndTime: 2020-04-11 22:53:33
	% DurationCPUTime: 0.20s
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
	% StartTime: 2020-04-11 22:53:33
	% EndTime: 2020-04-11 22:53:33
	% DurationCPUTime: 0.30s
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
	% StartTime: 2020-04-11 22:53:34
	% EndTime: 2020-04-11 22:53:34
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (18446->103), mult. (18896->236), div. (1256->13), fcn. (5114->11), ass. (0->150)
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
	t29 = 0.1e1 / t22;
	t30 = t29 * t27 * t1;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t33 = t27 ^ 2;
	t34 = 0.1e1 / t23;
	t36 = t31 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t47 = sqrt(t37 * t34 * t33 + t37 * t34 * t43);
	t48 = 0.1e1 / t47;
	t49 = t48 * t32;
	t50 = 0.1e1 / pkin(3);
	t51 = t50 * t29;
	t53 = t6 - t7 - t21 + t23 + t24;
	t56 = atan2(t17 * t51, t50 * t29 * t53);
	t57 = t56 + qJ(3);
	t58 = sin(t57);
	t59 = t58 * t49;
	t62 = t29 * t42 * t1;
	t63 = cos(t57);
	t64 = t63 * t49;
	t66 = t59 * t30 - t64 * t62;
	t67 = cos(qJ(1));
	t68 = 0.1e1 / t17;
	t69 = t68 * t4;
	t74 = -0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t84 = 0.1e1 / t42;
	t87 = 0.1e1 / t43;
	t90 = 0.1e1 / (t87 * t33 + 0.1e1);
	t91 = t90 * t31 * t22;
	t93 = t68 * t19;
	t107 = t90 * t87 * t31;
	t109 = t91 * t84 * (t32 * t29 * (t74 * t69 / 0.2e1 + 0.2e1 * t22 * t20) - t32 * t34 * t27) - t107 * t22 * t27 * (t32 * t29 * (-t74 * pkin(1) * t93 / 0.2e1 + 0.2e1 * t22 * t4) - t32 * t34 * t42);
	t110 = t109 * t67;
	t111 = t42 * t110;
	t112 = t32 * t29;
	t114 = t58 * t48 * t112;
	t117 = t29 * t27 * t67;
	t124 = 0.1e1 / t53;
	t127 = t53 ^ 2;
	t128 = 0.1e1 / t127;
	t131 = 0.1e1 / (-t128 * t16 + 0.1e1);
	t142 = t131 * t128 * t17;
	t144 = t131 * pkin(3) * t22 * t124 * (-t17 * t50 * t34 + t74 * t68 * t51 / 0.2e1) - t142 * pkin(3) * t22 * (0.2e1 * t50 * t29 * t22 - t50 * t34 * t53);
	t146 = t63 * t144 * t49;
	t148 = t27 * t110;
	t150 = t63 * t48 * t112;
	t153 = t29 * t42 * t67;
	t155 = t58 * t144 * t49;
	t160 = -t64 * t117 - t59 * t153;
	t164 = pkin(1) * pkin(2);
	t166 = t15 * pkin(2) * t20 + t164 * t19 * t11;
	t170 = t19 ^ 2;
	t177 = t29 * t27;
	t179 = pkin(2) * t20;
	t196 = t29 * t42;
	t204 = t91 * t84 * (t32 * t29 * (-0.2e1 * pkin(2) * t170 * t7 + t166 * t69 - t25 * t3 - t40) + 0.2e1 * t179 * t37 * t177) - t107 * t22 * t27 * (t32 * t29 * (-t166 * pkin(1) * t93 + pkin(1) * t17 * t2 - 0.2e1 * t164 * t19 * t4 - t26) + 0.2e1 * t179 * t37 * t196);
	t205 = t204 * t67;
	t206 = t42 * t205;
	t214 = t131 * t124 * t166 * t68 - 0.2e1 * t142 * t179;
	t216 = t63 * t214 * t49;
	t218 = t27 * t205;
	t221 = t58 * t214 * t49;
	t226 = -t59 * t117 + t64 * t153;
	t227 = t109 * t1;
	t228 = t42 * t227;
	t231 = t27 * t227;
	t237 = -t64 * t30 - t59 * t62;
	t238 = t204 * t1;
	t239 = t42 * t238;
	t242 = t27 * t238;
	t247 = t29 * t27 * t109;
	t249 = t32 * t196;
	t250 = t144 * t48;
	t251 = t63 * t250;
	t254 = t29 * t42 * t109;
	t256 = t32 * t177;
	t257 = t58 * t250;
	t264 = t29 * t27 * t204;
	t266 = t214 * t48;
	t267 = t63 * t266;
	t270 = t29 * t42 * t204;
	t272 = t58 * t266;
	unknown(1,1) = t66;
	unknown(1,2) = -t114 * t111 - t146 * t117 - t150 * t148 - t155 * t153;
	unknown(1,3) = t160;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = -t114 * t206 - t216 * t117 - t150 * t218 - t221 * t153;
	unknown(2,1) = t226;
	unknown(2,2) = -t114 * t228 - t146 * t30 - t150 * t231 - t155 * t62;
	unknown(2,3) = t237;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = -t114 * t239 - t150 * t242 - t216 * t30 - t221 * t62;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t59 * t247 + t251 * t249 + t64 * t254 - t257 * t256;
	unknown(3,3) = -t59 * t177 + t64 * t196;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = t267 * t249 - t272 * t256 - t59 * t264 + t64 * t270;
	unknown(4,1) = -t237;
	unknown(4,2) = -t150 * t111 + t114 * t148 + t155 * t117 - t146 * t153;
	unknown(4,3) = -t226;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = t114 * t218 + t221 * t117 - t150 * t206 - t216 * t153;
	unknown(5,1) = t160;
	unknown(5,2) = t114 * t231 - t146 * t62 - t150 * t228 + t155 * t30;
	unknown(5,3) = t66;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = t114 * t242 - t150 * t239 - t216 * t62 + t221 * t30;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = -t64 * t247 - t257 * t249 - t251 * t256 - t59 * t254;
	unknown(6,3) = -t64 * t177 - t59 * t196;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = -t272 * t249 - t267 * t256 - t64 * t264 - t59 * t270;
	unknown(7,1) = t67;
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
	% StartTime: 2020-04-11 22:53:36
	% EndTime: 2020-04-11 22:53:36
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (28976->119), mult. (29732->266), div. (2004->13), fcn. (8116->13), ass. (0->164)
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
	t29 = 0.1e1 / t22;
	t30 = t29 * t27 * t1;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t33 = t27 ^ 2;
	t34 = 0.1e1 / t23;
	t36 = t31 ^ 2;
	t37 = 0.1e1 / t36;
	t40 = pkin(1) * t17 * t19;
	t42 = t25 * t4 - t40;
	t43 = t42 ^ 2;
	t47 = sqrt(t37 * t34 * t33 + t37 * t34 * t43);
	t48 = 0.1e1 / t47;
	t49 = t48 * t32;
	t50 = 0.1e1 / pkin(3);
	t51 = t50 * t29;
	t53 = t6 - t7 - t21 + t23 + t24;
	t56 = atan2(t17 * t51, t50 * t29 * t53);
	t57 = t56 + qJ(3);
	t58 = cos(t57);
	t59 = t58 * t49;
	t62 = t29 * t42 * t1;
	t63 = sin(t57);
	t64 = t63 * t49;
	t66 = t59 * t30 + t64 * t62;
	t67 = cos(qJ(4));
	t69 = cos(qJ(1));
	t70 = sin(qJ(4));
	t71 = t70 * t69;
	t73 = 0.1e1 / t17;
	t74 = t73 * t4;
	t79 = -0.2e1 * t15 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t11;
	t89 = 0.1e1 / t42;
	t92 = 0.1e1 / t43;
	t95 = 0.1e1 / (t92 * t33 + 0.1e1);
	t96 = t95 * t31 * t22;
	t98 = t73 * t19;
	t112 = t95 * t92 * t31;
	t114 = t96 * t89 * (t32 * t29 * (t79 * t74 / 0.2e1 + 0.2e1 * t22 * t20) - t32 * t34 * t27) - t112 * t22 * t27 * (t32 * t29 * (-t79 * pkin(1) * t98 / 0.2e1 + 0.2e1 * t22 * t4) - t32 * t34 * t42);
	t115 = t114 * t69;
	t116 = t42 * t115;
	t117 = t32 * t29;
	t119 = t58 * t48 * t117;
	t122 = t29 * t27 * t69;
	t129 = 0.1e1 / t53;
	t132 = t53 ^ 2;
	t133 = 0.1e1 / t132;
	t136 = 0.1e1 / (-t133 * t16 + 0.1e1);
	t147 = t136 * t133 * t17;
	t149 = t136 * pkin(3) * t22 * t129 * (-t17 * t50 * t34 + t79 * t73 * t51 / 0.2e1) - t147 * pkin(3) * t22 * (0.2e1 * t50 * t29 * t22 - t50 * t34 * t53);
	t151 = t63 * t149 * t49;
	t153 = t27 * t115;
	t155 = t63 * t48 * t117;
	t158 = t29 * t42 * t69;
	t160 = t58 * t149 * t49;
	t162 = -t119 * t116 + t151 * t122 + t155 * t153 - t160 * t158;
	t166 = t64 * t122 - t59 * t158;
	t170 = -t59 * t122 - t64 * t158;
	t173 = -t67 * t1 + t70 * t170;
	t177 = pkin(1) * pkin(2);
	t179 = t15 * pkin(2) * t20 + t177 * t19 * t11;
	t183 = t19 ^ 2;
	t190 = t29 * t27;
	t192 = pkin(2) * t20;
	t209 = t29 * t42;
	t217 = t96 * t89 * (t32 * t29 * (-0.2e1 * pkin(2) * t183 * t7 + t179 * t74 - t25 * t3 - t40) + 0.2e1 * t192 * t37 * t190) - t112 * t22 * t27 * (t32 * t29 * (pkin(1) * t17 * t2 - t179 * pkin(1) * t98 - 0.2e1 * t177 * t19 * t4 - t26) + 0.2e1 * t192 * t37 * t209);
	t218 = t217 * t69;
	t219 = t42 * t218;
	t227 = t136 * t129 * t179 * t73 - 0.2e1 * t147 * t192;
	t229 = t63 * t227 * t49;
	t231 = t27 * t218;
	t234 = t58 * t227 * t49;
	t236 = -t119 * t219 + t229 * t122 + t155 * t231 - t234 * t158;
	t240 = -t70 * t1 - t67 * t170;
	t241 = t114 * t1;
	t242 = t42 * t241;
	t245 = t27 * t241;
	t248 = -t119 * t242 + t151 * t30 + t155 * t245 - t160 * t62;
	t252 = t64 * t30 - t59 * t62;
	t255 = t67 * t69;
	t257 = t217 * t1;
	t258 = t42 * t257;
	t261 = t27 * t257;
	t264 = -t119 * t258 + t155 * t261 + t229 * t30 - t234 * t62;
	t267 = t29 * t27 * t114;
	t269 = t32 * t209;
	t270 = t149 * t48;
	t271 = t63 * t270;
	t274 = t29 * t42 * t114;
	t276 = t32 * t190;
	t277 = t58 * t270;
	t279 = -t59 * t267 - t271 * t269 - t64 * t274 - t277 * t276;
	t283 = -t59 * t190 - t64 * t209;
	t287 = -t64 * t190 + t59 * t209;
	t290 = t29 * t27 * t217;
	t292 = t227 * t48;
	t293 = t63 * t292;
	t296 = t29 * t42 * t217;
	t298 = t58 * t292;
	t300 = -t293 * t269 - t298 * t276 - t59 * t290 - t64 * t296;
	unknown(1,1) = -t67 * t66 - t71;
	unknown(1,2) = -t67 * t162;
	unknown(1,3) = -t67 * t166;
	unknown(1,4) = t173;
	unknown(1,5) = -t67 * t236;
	unknown(2,1) = t240;
	unknown(2,2) = -t67 * t248;
	unknown(2,3) = -t67 * t252;
	unknown(2,4) = -t70 * t66 + t255;
	unknown(2,5) = -t67 * t264;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t67 * t279;
	unknown(3,3) = -t67 * t283;
	unknown(3,4) = t70 * t287;
	unknown(3,5) = -t67 * t300;
	unknown(4,1) = t70 * t66 - t255;
	unknown(4,2) = t70 * t162;
	unknown(4,3) = t70 * t166;
	unknown(4,4) = -t240;
	unknown(4,5) = t70 * t236;
	unknown(5,1) = t173;
	unknown(5,2) = t70 * t248;
	unknown(5,3) = t70 * t252;
	unknown(5,4) = -t67 * t66 - t71;
	unknown(5,5) = t70 * t264;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t70 * t279;
	unknown(6,3) = t70 * t283;
	unknown(6,4) = t67 * t287;
	unknown(6,5) = t70 * t300;
	unknown(7,1) = t252;
	unknown(7,2) = -t155 * t116 - t119 * t153 - t160 * t122 - t151 * t158;
	unknown(7,3) = t170;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = -t119 * t231 - t234 * t122 - t155 * t219 - t229 * t158;
	unknown(8,1) = -t166;
	unknown(8,2) = -t119 * t245 - t151 * t62 - t155 * t242 - t160 * t30;
	unknown(8,3) = -t66;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = -t119 * t261 - t155 * t258 - t229 * t62 - t234 * t30;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = -t64 * t267 + t277 * t269 - t271 * t276 + t59 * t274;
	unknown(9,3) = t287;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = t298 * t269 - t293 * t276 - t64 * t290 + t59 * t296;
	JR_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
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