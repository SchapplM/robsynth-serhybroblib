% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh4m1TE
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
% Datum: 2020-04-11 21:48
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh4m1TE_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1TE_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1TE_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1TE_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
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
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
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
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:19
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (1086->54), mult. (1180->132), div. (56->5), fcn. (321->6), ass. (0->91)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = cos(qJ(5));
	t3 = sin(qJ(5));
	t4 = t3 * pkin(1);
	t6 = 0.2e1 * t4 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = t2 * t17 * pkin(1);
	t20 = t4 - pkin(2);
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t27 = -t20 * t25 - t19;
	t28 = t1 * t27;
	t29 = 0.1e1 / t22;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t32 = t29 * t31;
	t35 = cos(qJ(1));
	t36 = 0.1e1 / t17;
	t37 = t2 * t36;
	t42 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t47 = -t37 * pkin(1) * t42 / 0.2e1 - 0.2e1 * t20 * t22;
	t50 = t35 * t27;
	t51 = 0.1e1 / t23;
	t52 = t51 * t31;
	t57 = pkin(1) * t2;
	t61 = pkin(1) * pkin(2);
	t63 = pkin(2) * t15 * t57 + t11 * t2 * t61;
	t67 = t57 * t25;
	t71 = pkin(1) * t17 * t3 - pkin(1) * t37 * t63 + 0.2e1 * t2 * t20 * t61 - t67;
	t76 = t30 ^ 2;
	t77 = 0.1e1 / t76;
	t79 = t77 * t2 * t61;
	t94 = -t20 * t36;
	t98 = t94 * t42 / 0.2e1 + 0.2e1 * t57 * t22;
	t102 = -t17 * t20 + t67;
	t109 = t2 ^ 2;
	t113 = -0.2e1 * pkin(2) * t109 * t7 - t25 * t4 + t63 * t94 - t19;
	t119 = t57 * pkin(2);
	t122 = t1 * t102;
	t127 = t35 * t102;
	unknown(1,1) = -t28 * t32 / 0.2e1;
	unknown(1,2) = t35 * t47 * t32 / 0.2e1 - t50 * t52 / 0.2e1;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t35 * t71 * t32 / 0.2e1 + t50 * t29 * t79;
	unknown(2,1) = t50 * t32 / 0.2e1;
	unknown(2,2) = t1 * t47 * t32 / 0.2e1 - t28 * t52 / 0.2e1;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t1 * t71 * t32 / 0.2e1 + t28 * t29 * t79;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t102 * t51 * t31 / 0.2e1 + t98 * t29 * t31 / 0.2e1;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = t113 * t29 * t31 / 0.2e1 + t102 * t29 * t77 * t119;
	unknown(4,1) = t122 * t32 / 0.2e1;
	unknown(4,2) = -t35 * t98 * t32 / 0.2e1 + t127 * t52 / 0.2e1;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = -t35 * t113 * t32 / 0.2e1 - t127 * t29 * t79;
	unknown(5,1) = -t127 * t32 / 0.2e1;
	unknown(5,2) = -t1 * t98 * t32 / 0.2e1 + t122 * t52 / 0.2e1;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = -t1 * t113 * t32 / 0.2e1 - t122 * t29 * t79;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = -t27 * t51 * t31 / 0.2e1 + t47 * t29 * t31 / 0.2e1;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = t71 * t29 * t31 / 0.2e1 + t27 * t29 * t77 * t119;
	unknown(7,1) = t35;
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
	% StartTime: 2020-04-11 21:48:19
	% EndTime: 2020-04-11 21:48:19
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (1088->55), mult. (1180->132), div. (56->5), fcn. (321->6), ass. (0->91)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = t2 * pkin(1);
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * t3 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = cos(qJ(5));
	t20 = pkin(1) * t19;
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t20 * t25;
	t27 = t4 * t17 + t26;
	t28 = t1 * t27;
	t29 = 0.1e1 / t22;
	t30 = -t6 + t7 + t21;
	t31 = 0.1e1 / t30;
	t32 = t29 * t31;
	t35 = cos(qJ(1));
	t36 = 0.1e1 / t17;
	t37 = t4 * t36;
	t42 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t46 = t37 * t42 / 0.2e1 + 0.2e1 * t20 * t22;
	t49 = t35 * t27;
	t50 = 0.1e1 / t23;
	t51 = t50 * t31;
	t55 = t19 * t17 * pkin(1);
	t59 = pkin(1) * pkin(2);
	t61 = t20 * pkin(2) * t15 + t11 * t19 * t59;
	t65 = t19 ^ 2;
	t69 = -0.2e1 * t7 * t65 * pkin(2) - t3 * t25 + t37 * t61 - t55;
	t74 = t30 ^ 2;
	t75 = 0.1e1 / t74;
	t77 = t75 * t19 * t59;
	t92 = t19 * t36;
	t97 = -t92 * pkin(1) * t42 / 0.2e1 + 0.2e1 * t4 * t22;
	t101 = t4 * t25 - t55;
	t113 = t2 * t17 * pkin(1) - t92 * pkin(1) * t61 - 0.2e1 * t4 * t19 * t59 - t26;
	t119 = t20 * pkin(2);
	t122 = t1 * t101;
	t127 = t35 * t101;
	unknown(1,1) = -t28 * t32 / 0.2e1;
	unknown(1,2) = t35 * t46 * t32 / 0.2e1 - t49 * t51 / 0.2e1;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t35 * t69 * t32 / 0.2e1 + t49 * t29 * t77;
	unknown(2,1) = t49 * t32 / 0.2e1;
	unknown(2,2) = t1 * t46 * t32 / 0.2e1 - t28 * t51 / 0.2e1;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t1 * t69 * t32 / 0.2e1 + t28 * t29 * t77;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t101 * t50 * t31 / 0.2e1 - t97 * t29 * t31 / 0.2e1;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = -t113 * t29 * t31 / 0.2e1 - t101 * t29 * t75 * t119;
	unknown(4,1) = -t35;
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
	unknown(7,1) = -t122 * t32 / 0.2e1;
	unknown(7,2) = t35 * t97 * t32 / 0.2e1 - t127 * t51 / 0.2e1;
	unknown(7,3) = 0.0e0;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = t35 * t113 * t32 / 0.2e1 + t127 * t29 * t77;
	unknown(8,1) = t127 * t32 / 0.2e1;
	unknown(8,2) = t1 * t97 * t32 / 0.2e1 - t122 * t51 / 0.2e1;
	unknown(8,3) = 0.0e0;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = t1 * t113 * t32 / 0.2e1 + t122 * t29 * t77;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = -t27 * t50 * t31 / 0.2e1 + t46 * t29 * t31 / 0.2e1;
	unknown(9,3) = 0.0e0;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = t69 * t29 * t31 / 0.2e1 + t27 * t29 * t75 * t119;
	JR_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:19
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (3914->109), mult. (4172->288), div. (240->7), fcn. (1078->6), ass. (0->138)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = t2 * pkin(1);
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * t3 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = cos(qJ(5));
	t20 = pkin(1) * t19;
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t20 * t25;
	t27 = t4 * t17 + t26;
	t28 = t1 * t27;
	t29 = 0.1e1 / t23;
	t30 = t28 * t29;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t33 = 0.1e1 / pkin(3);
	t34 = t32 * t33;
	t35 = t34 * t17;
	t37 = t19 * t17;
	t38 = t37 * pkin(1);
	t40 = t4 * t25 - t38;
	t41 = t1 * t40;
	t42 = t41 * t29;
	t43 = t6 - t7 - t21 + t23 + t24;
	t45 = t32 * t43 * t33;
	t48 = cos(qJ(1));
	t49 = 0.1e1 / t17;
	t50 = t4 * t49;
	t55 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t59 = t50 * t55 / 0.2e1 + 0.2e1 * t20 * t22;
	t61 = t48 * t59 * t29;
	t64 = t48 * t27;
	t66 = 0.1e1 / t23 / t22;
	t67 = t64 * t66;
	t70 = t64 * t29;
	t72 = t34 * t49 * t55;
	t75 = t19 * t49;
	t80 = -t75 * pkin(1) * t55 / 0.2e1 + 0.2e1 * t4 * t22;
	t82 = t48 * t80 * t29;
	t85 = t48 * t40;
	t86 = t85 * t66;
	t89 = t85 * t29;
	t91 = 0.2e1 * t32 * t22 * t33;
	t98 = pkin(1) * pkin(2);
	t100 = t20 * pkin(2) * t15 + t11 * t19 * t98;
	t104 = t19 ^ 2;
	t108 = -0.2e1 * t7 * t104 * pkin(2) + t50 * t100 - t3 * t25 - t38;
	t110 = t48 * t108 * t29;
	t113 = t31 ^ 2;
	t114 = 0.1e1 / t113;
	t115 = t29 * t114;
	t116 = t64 * t115;
	t118 = t20 * pkin(2);
	t119 = t33 * t17 * t118;
	t123 = 0.2e1 * t34 * t49 * t100;
	t134 = -t75 * pkin(1) * t100 + t2 * t17 * pkin(1) - 0.2e1 * t4 * t19 * t98 - t26;
	t136 = t48 * t134 * t29;
	t139 = t85 * t115;
	t141 = t43 * t33 * t118;
	t144 = t29 * t32;
	t147 = t20 * pkin(2) * t33;
	t155 = t1 * t59 * t29;
	t158 = t28 * t66;
	t164 = t1 * t80 * t29;
	t167 = t41 * t66;
	t174 = t1 * t108 * t29;
	t177 = t28 * t115;
	t183 = t1 * t134 * t29;
	t186 = t41 * t115;
	t193 = t80 * t29;
	t196 = t40 * t66;
	t199 = t40 * t29;
	t200 = t199 * t32;
	t201 = t33 * t49;
	t202 = t201 * t55;
	t205 = t59 * t29;
	t208 = t27 * t66;
	t211 = t27 * t29;
	t215 = t134 * t29;
	t218 = t114 * t33;
	t220 = t37 * t98;
	t223 = 0.2e1 * t201 * t100;
	t226 = t108 * t29;
	t229 = t114 * t43;
	t233 = t211 * t32;
	unknown(1,1) = t30 * t35 / 0.4e1 - t42 * t45 / 0.4e1;
	unknown(1,2) = -t61 * t35 / 0.4e1 + t67 * t35 / 0.2e1 - t70 * t72 / 0.8e1 + t82 * t45 / 0.4e1 - t86 * t45 / 0.2e1 + t89 * t91 / 0.4e1;
	unknown(1,3) = 0.0e0;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = -t110 * t35 / 0.4e1 - t116 * t119 / 0.2e1 - t70 * t123 / 0.8e1 + t136 * t45 / 0.4e1 + t139 * t141 / 0.2e1 + t85 * t144 * t147 / 0.2e1;
	unknown(2,1) = -t70 * t35 / 0.4e1 + t89 * t45 / 0.4e1;
	unknown(2,2) = -t155 * t35 / 0.4e1 + t158 * t35 / 0.2e1 - t30 * t72 / 0.8e1 + t164 * t45 / 0.4e1 - t167 * t45 / 0.2e1 + t42 * t91 / 0.4e1;
	unknown(2,3) = 0.0e0;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = -t174 * t35 / 0.4e1 - t177 * t119 / 0.2e1 - t30 * t123 / 0.8e1 + t183 * t45 / 0.4e1 + t186 * t141 / 0.2e1 + t41 * t144 * t147 / 0.2e1;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t193 * t35 / 0.4e1 - t196 * t35 / 0.2e1 + t200 * t202 / 0.8e1 + t205 * t45 / 0.4e1 - t208 * t45 / 0.2e1 + t211 * t91 / 0.4e1;
	unknown(3,3) = 0.0e0;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = t215 * t35 / 0.4e1 + t199 * t218 * t220 / 0.2e1 + t200 * t223 / 0.8e1 + t226 * t45 / 0.4e1 + t211 * t229 * t147 / 0.2e1 + t233 * t147 / 0.2e1;
	unknown(4,1) = t30 * t45 / 0.4e1 + t42 * t35 / 0.4e1;
	unknown(4,2) = -t61 * t45 / 0.4e1 + t67 * t45 / 0.2e1 - t70 * t91 / 0.4e1 - t82 * t35 / 0.4e1 + t86 * t35 / 0.2e1 - t89 * t72 / 0.8e1;
	unknown(4,3) = 0.0e0;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = -t110 * t45 / 0.4e1 - t116 * t141 / 0.2e1 - t64 * t144 * t147 / 0.2e1 - t136 * t35 / 0.4e1 - t139 * t119 / 0.2e1 - t89 * t123 / 0.8e1;
	unknown(5,1) = -t89 * t35 / 0.4e1 - t70 * t45 / 0.4e1;
	unknown(5,2) = -t155 * t45 / 0.4e1 + t158 * t45 / 0.2e1 - t30 * t91 / 0.4e1 - t164 * t35 / 0.4e1 + t167 * t35 / 0.2e1 - t72 * t42 / 0.8e1;
	unknown(5,3) = 0.0e0;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = -t174 * t45 / 0.4e1 - t177 * t141 / 0.2e1 - t28 * t144 * t147 / 0.2e1 - t183 * t35 / 0.4e1 - t186 * t119 / 0.2e1 - t42 * t123 / 0.8e1;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t193 * t45 / 0.4e1 - t196 * t45 / 0.2e1 + t199 * t91 / 0.4e1 - t205 * t35 / 0.4e1 + t208 * t35 / 0.2e1 - t233 * t202 / 0.8e1;
	unknown(6,3) = 0.0e0;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = t215 * t45 / 0.4e1 + t199 * t229 * t147 / 0.2e1 + t200 * t147 / 0.2e1 - t226 * t35 / 0.4e1 - t211 * t218 * t220 / 0.2e1 - t233 * t223 / 0.8e1;
	unknown(7,1) = t48;
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
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:19
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (8774->132), mult. (9376->340), div. (552->7), fcn. (2466->8), ass. (0->160)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = t2 * pkin(1);
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * t3 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = cos(qJ(5));
	t20 = pkin(1) * t19;
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t20 * t25;
	t27 = t4 * t17 + t26;
	t28 = t1 * t27;
	t29 = 0.1e1 / t23;
	t30 = t28 * t29;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t33 = 0.1e1 / pkin(3);
	t34 = t32 * t33;
	t35 = t34 * t17;
	t37 = t19 * t17;
	t38 = t37 * pkin(1);
	t40 = t4 * t25 - t38;
	t41 = t1 * t40;
	t42 = t41 * t29;
	t43 = t6 - t7 - t21 + t23 + t24;
	t45 = t32 * t43 * t33;
	t47 = t30 * t35 - t42 * t45;
	t48 = cos(qJ(3));
	t52 = t30 * t45 + t42 * t35;
	t53 = sin(qJ(3));
	t56 = cos(qJ(1));
	t57 = 0.1e1 / t17;
	t58 = t4 * t57;
	t63 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t67 = t58 * t63 / 0.2e1 + 0.2e1 * t20 * t22;
	t69 = t56 * t67 * t29;
	t72 = t56 * t27;
	t74 = 0.1e1 / t23 / t22;
	t75 = t72 * t74;
	t78 = t72 * t29;
	t80 = t34 * t57 * t63;
	t83 = t19 * t57;
	t88 = -t83 * pkin(1) * t63 / 0.2e1 + 0.2e1 * t4 * t22;
	t90 = t56 * t88 * t29;
	t93 = t56 * t40;
	t94 = t93 * t74;
	t97 = t93 * t29;
	t99 = 0.2e1 * t32 * t22 * t33;
	t102 = -t69 * t35 / 0.4e1 + t75 * t35 / 0.2e1 - t78 * t80 / 0.8e1 + t90 * t45 / 0.4e1 - t94 * t45 / 0.2e1 + t97 * t99 / 0.4e1;
	t116 = -t69 * t45 / 0.4e1 + t75 * t45 / 0.2e1 - t78 * t99 / 0.4e1 - t90 * t35 / 0.4e1 + t94 * t35 / 0.2e1 - t97 * t80 / 0.8e1;
	t121 = -t78 * t35 + t97 * t45;
	t125 = -t97 * t35 - t78 * t45;
	t127 = -t121 * t53 / 0.4e1 + t125 * t48 / 0.4e1;
	t131 = pkin(1) * pkin(2);
	t133 = t20 * pkin(2) * t15 + t11 * t19 * t131;
	t137 = t19 ^ 2;
	t141 = -0.2e1 * t7 * t137 * pkin(2) + t58 * t133 - t3 * t25 - t38;
	t143 = t56 * t141 * t29;
	t146 = t31 ^ 2;
	t147 = 0.1e1 / t146;
	t148 = t29 * t147;
	t149 = t72 * t148;
	t151 = t20 * pkin(2);
	t152 = t33 * t17 * t151;
	t156 = 0.2e1 * t34 * t57 * t133;
	t167 = -t83 * pkin(1) * t133 + t2 * t17 * pkin(1) - 0.2e1 * t4 * t19 * t131 - t26;
	t169 = t56 * t167 * t29;
	t172 = t93 * t148;
	t174 = t43 * t33 * t151;
	t177 = t29 * t32;
	t180 = t20 * pkin(2) * t33;
	t183 = -t143 * t35 / 0.4e1 - t149 * t152 / 0.2e1 - t78 * t156 / 0.8e1 + t169 * t45 / 0.4e1 + t172 * t174 / 0.2e1 + t93 * t177 * t180 / 0.2e1;
	t198 = -t143 * t45 / 0.4e1 - t149 * t174 / 0.2e1 - t72 * t177 * t180 / 0.2e1 - t169 * t35 / 0.4e1 - t172 * t152 / 0.2e1 - t97 * t156 / 0.8e1;
	t203 = t121 * t48 / 0.4e1 + t125 * t53 / 0.4e1;
	t205 = t1 * t67 * t29;
	t208 = t28 * t74;
	t214 = t1 * t88 * t29;
	t217 = t41 * t74;
	t222 = -t205 * t35 / 0.4e1 + t208 * t35 / 0.2e1 - t30 * t80 / 0.8e1 + t214 * t45 / 0.4e1 - t217 * t45 / 0.2e1 + t42 * t99 / 0.4e1;
	t236 = -t205 * t45 / 0.4e1 + t208 * t45 / 0.2e1 - t30 * t99 / 0.4e1 - t214 * t35 / 0.4e1 + t217 * t35 / 0.2e1 - t42 * t80 / 0.8e1;
	t243 = t1 * t141 * t29;
	t246 = t28 * t148;
	t252 = t1 * t167 * t29;
	t255 = t41 * t148;
	t261 = -t243 * t35 / 0.4e1 - t246 * t152 / 0.2e1 - t30 * t156 / 0.8e1 + t252 * t45 / 0.4e1 + t255 * t174 / 0.2e1 + t41 * t177 * t180 / 0.2e1;
	t276 = -t243 * t45 / 0.4e1 - t246 * t174 / 0.2e1 - t28 * t177 * t180 / 0.2e1 - t252 * t35 / 0.4e1 - t255 * t152 / 0.2e1 - t42 * t156 / 0.8e1;
	t279 = t88 * t29;
	t282 = t40 * t74;
	t285 = t40 * t29;
	t286 = t285 * t32;
	t287 = t33 * t57;
	t288 = t287 * t63;
	t291 = t67 * t29;
	t294 = t27 * t74;
	t297 = t27 * t29;
	t300 = t279 * t35 / 0.4e1 - t282 * t35 / 0.2e1 + t286 * t288 / 0.8e1 + t291 * t45 / 0.4e1 - t294 * t45 / 0.2e1 + t297 * t99 / 0.4e1;
	t312 = t297 * t32;
	t315 = t279 * t45 / 0.4e1 - t282 * t45 / 0.2e1 + t285 * t99 / 0.4e1 - t291 * t35 / 0.4e1 + t294 * t35 / 0.2e1 - t312 * t288 / 0.8e1;
	t320 = t285 * t35 + t297 * t45;
	t324 = t285 * t45 - t297 * t35;
	t327 = t167 * t29;
	t330 = t147 * t33;
	t332 = t37 * t131;
	t335 = 0.2e1 * t287 * t133;
	t338 = t141 * t29;
	t341 = t147 * t43;
	t347 = t327 * t35 / 0.4e1 + t285 * t330 * t332 / 0.2e1 + t286 * t335 / 0.8e1 + t338 * t45 / 0.4e1 + t297 * t341 * t180 / 0.2e1 + t312 * t180 / 0.2e1;
	t363 = t327 * t45 / 0.4e1 + t285 * t341 * t180 / 0.2e1 + t286 * t180 / 0.2e1 - t338 * t35 / 0.4e1 - t297 * t330 * t332 / 0.2e1 - t312 * t335 / 0.8e1;
	unknown(1,1) = t47 * t48 / 0.4e1 + t52 * t53 / 0.4e1;
	unknown(1,2) = t102 * t48 + t116 * t53;
	unknown(1,3) = t127;
	unknown(1,4) = 0.0e0;
	unknown(1,5) = t183 * t48 + t198 * t53;
	unknown(2,1) = t203;
	unknown(2,2) = t222 * t48 + t236 * t53;
	unknown(2,3) = t47 * t53 / 0.4e1 - t52 * t48 / 0.4e1;
	unknown(2,4) = 0.0e0;
	unknown(2,5) = t261 * t48 + t276 * t53;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = t300 * t48 + t315 * t53;
	unknown(3,3) = -t320 * t53 / 0.4e1 + t324 * t48 / 0.4e1;
	unknown(3,4) = 0.0e0;
	unknown(3,5) = t347 * t48 + t363 * t53;
	unknown(4,1) = -t47 * t53 / 0.4e1 + t52 * t48 / 0.4e1;
	unknown(4,2) = -t102 * t53 + t116 * t48;
	unknown(4,3) = -t203;
	unknown(4,4) = 0.0e0;
	unknown(4,5) = -t183 * t53 + t198 * t48;
	unknown(5,1) = t127;
	unknown(5,2) = -t222 * t53 + t236 * t48;
	unknown(5,3) = t47 * t48 / 0.4e1 + t52 * t53 / 0.4e1;
	unknown(5,4) = 0.0e0;
	unknown(5,5) = -t261 * t53 + t276 * t48;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = -t300 * t53 + t315 * t48;
	unknown(6,3) = -t320 * t48 / 0.4e1 - t324 * t53 / 0.4e1;
	unknown(6,4) = 0.0e0;
	unknown(6,5) = -t347 * t53 + t363 * t48;
	unknown(7,1) = t56;
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
	% StartTime: 2020-04-11 21:48:22
	% EndTime: 2020-04-11 21:48:22
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (14108->148), mult. (15100->370), div. (900->7), fcn. (4020->10), ass. (0->177)
	unknown=NaN(9,5);
	t1 = sin(qJ(1));
	t2 = sin(qJ(5));
	t3 = t2 * pkin(1);
	t4 = -t3 + pkin(2);
	t6 = 0.2e1 * t3 * pkin(2);
	t7 = pkin(1) ^ 2;
	t11 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t15 = -t6 + t7 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t17 = sqrt(-t11 * t15);
	t19 = cos(qJ(5));
	t20 = pkin(1) * t19;
	t21 = pkin(2) ^ 2;
	t22 = qJ(2) + pkin(6);
	t23 = t22 ^ 2;
	t24 = pkin(3) ^ 2;
	t25 = -t6 + t7 + t21 + t23 - t24;
	t26 = t20 * t25;
	t27 = t4 * t17 + t26;
	t28 = t1 * t27;
	t29 = 0.1e1 / t23;
	t30 = t28 * t29;
	t31 = -t6 + t7 + t21;
	t32 = 0.1e1 / t31;
	t33 = 0.1e1 / pkin(3);
	t34 = t32 * t33;
	t35 = t34 * t17;
	t37 = t19 * t17;
	t38 = t37 * pkin(1);
	t40 = t4 * t25 - t38;
	t41 = t1 * t40;
	t42 = t41 * t29;
	t43 = t6 - t7 - t21 + t23 + t24;
	t45 = t32 * t43 * t33;
	t47 = t30 * t35 - t42 * t45;
	t48 = sin(qJ(3));
	t52 = t30 * t45 + t42 * t35;
	t53 = cos(qJ(3));
	t55 = -t47 * t48 / 0.4e1 + t52 * t53 / 0.4e1;
	t56 = cos(qJ(4));
	t58 = cos(qJ(1));
	t59 = sin(qJ(4));
	t60 = t58 * t59;
	t62 = 0.1e1 / t17;
	t63 = t4 * t62;
	t68 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t15 - 0.2e1 * t11 * (-qJ(2) - pkin(6) + pkin(3));
	t72 = t63 * t68 / 0.2e1 + 0.2e1 * t20 * t22;
	t74 = t58 * t72 * t29;
	t77 = t58 * t27;
	t79 = 0.1e1 / t23 / t22;
	t80 = t77 * t79;
	t83 = t77 * t29;
	t85 = t34 * t62 * t68;
	t88 = t19 * t62;
	t93 = -t88 * pkin(1) * t68 / 0.2e1 + 0.2e1 * t4 * t22;
	t95 = t58 * t93 * t29;
	t98 = t58 * t40;
	t99 = t98 * t79;
	t102 = t98 * t29;
	t104 = 0.2e1 * t32 * t22 * t33;
	t107 = -t74 * t35 / 0.4e1 + t80 * t35 / 0.2e1 - t83 * t85 / 0.8e1 + t95 * t45 / 0.4e1 - t99 * t45 / 0.2e1 + t102 * t104 / 0.4e1;
	t121 = -t74 * t45 / 0.4e1 + t80 * t45 / 0.2e1 - t83 * t104 / 0.4e1 - t95 * t35 / 0.4e1 + t99 * t35 / 0.2e1 - t102 * t85 / 0.8e1;
	t123 = -t107 * t48 + t121 * t53;
	t127 = t102 * t45 - t83 * t35;
	t131 = -t102 * t35 - t83 * t45;
	t133 = -t127 * t53 / 0.4e1 - t131 * t48 / 0.4e1;
	t137 = -t127 * t48 / 0.4e1 + t131 * t53 / 0.4e1;
	t140 = -t1 * t56 + t137 * t59;
	t144 = pkin(1) * pkin(2);
	t146 = pkin(2) * t15 * t20 + t11 * t144 * t19;
	t150 = t19 ^ 2;
	t154 = -0.2e1 * pkin(2) * t150 * t7 + t63 * t146 - t3 * t25 - t38;
	t156 = t58 * t154 * t29;
	t159 = t31 ^ 2;
	t160 = 0.1e1 / t159;
	t161 = t29 * t160;
	t162 = t77 * t161;
	t164 = t20 * pkin(2);
	t165 = t33 * t17 * t164;
	t169 = 0.2e1 * t34 * t62 * t146;
	t180 = -pkin(1) * t146 * t88 + pkin(1) * t17 * t2 - 0.2e1 * t4 * t19 * t144 - t26;
	t182 = t58 * t180 * t29;
	t185 = t98 * t161;
	t187 = t43 * t33 * t164;
	t190 = t29 * t32;
	t193 = t20 * pkin(2) * t33;
	t196 = -t156 * t35 / 0.4e1 - t162 * t165 / 0.2e1 - t83 * t169 / 0.8e1 + t182 * t45 / 0.4e1 + t185 * t187 / 0.2e1 + t98 * t190 * t193 / 0.2e1;
	t211 = -t156 * t45 / 0.4e1 - t162 * t187 / 0.2e1 - t77 * t190 * t193 / 0.2e1 - t182 * t35 / 0.4e1 - t185 * t165 / 0.2e1 - t102 * t169 / 0.8e1;
	t213 = -t196 * t48 + t211 * t53;
	t217 = -t1 * t59 - t137 * t56;
	t219 = t1 * t72 * t29;
	t222 = t28 * t79;
	t228 = t1 * t93 * t29;
	t231 = t41 * t79;
	t236 = -t219 * t35 / 0.4e1 + t222 * t35 / 0.2e1 - t30 * t85 / 0.8e1 + t228 * t45 / 0.4e1 - t231 * t45 / 0.2e1 + t42 * t104 / 0.4e1;
	t250 = -t219 * t45 / 0.4e1 + t222 * t45 / 0.2e1 - t30 * t104 / 0.4e1 - t228 * t35 / 0.4e1 + t231 * t35 / 0.2e1 - t42 * t85 / 0.8e1;
	t252 = -t236 * t48 + t250 * t53;
	t256 = t47 * t53 / 0.4e1 + t52 * t48 / 0.4e1;
	t260 = t47 * t48 / 0.4e1 - t52 * t53 / 0.4e1;
	t262 = t58 * t56;
	t265 = t1 * t154 * t29;
	t268 = t28 * t161;
	t274 = t1 * t180 * t29;
	t277 = t41 * t161;
	t283 = -t265 * t35 / 0.4e1 - t268 * t165 / 0.2e1 - t30 * t169 / 0.8e1 + t274 * t45 / 0.4e1 + t277 * t187 / 0.2e1 + t41 * t190 * t193 / 0.2e1;
	t298 = -t265 * t45 / 0.4e1 - t268 * t187 / 0.2e1 - t28 * t190 * t193 / 0.2e1 - t274 * t35 / 0.4e1 - t277 * t165 / 0.2e1 - t42 * t169 / 0.8e1;
	t300 = -t283 * t48 + t298 * t53;
	t302 = t93 * t29;
	t305 = t40 * t79;
	t308 = t40 * t29;
	t309 = t308 * t32;
	t310 = t33 * t62;
	t311 = t310 * t68;
	t314 = t72 * t29;
	t317 = t27 * t79;
	t320 = t27 * t29;
	t323 = t302 * t35 / 0.4e1 - t305 * t35 / 0.2e1 + t309 * t311 / 0.8e1 + t314 * t45 / 0.4e1 - t317 * t45 / 0.2e1 + t320 * t104 / 0.4e1;
	t335 = t320 * t32;
	t338 = t302 * t45 / 0.4e1 - t305 * t45 / 0.2e1 + t308 * t104 / 0.4e1 - t314 * t35 / 0.4e1 + t317 * t35 / 0.2e1 - t335 * t311 / 0.8e1;
	t340 = -t323 * t48 + t338 * t53;
	t344 = t308 * t35 + t320 * t45;
	t348 = t308 * t45 - t320 * t35;
	t350 = -t344 * t53 / 0.4e1 - t348 * t48 / 0.4e1;
	t354 = -t344 * t48 / 0.4e1 + t348 * t53 / 0.4e1;
	t356 = t180 * t29;
	t359 = t160 * t33;
	t361 = t37 * t144;
	t364 = 0.2e1 * t310 * t146;
	t367 = t154 * t29;
	t370 = t160 * t43;
	t376 = t356 * t35 / 0.4e1 + t308 * t359 * t361 / 0.2e1 + t309 * t364 / 0.8e1 + t367 * t45 / 0.4e1 + t320 * t370 * t193 / 0.2e1 + t335 * t193 / 0.2e1;
	t392 = t356 * t45 / 0.4e1 + t308 * t370 * t193 / 0.2e1 + t309 * t193 / 0.2e1 - t367 * t35 / 0.4e1 - t320 * t359 * t361 / 0.2e1 - t335 * t364 / 0.8e1;
	t394 = -t376 * t48 + t392 * t53;
	unknown(1,1) = -t55 * t56 - t60;
	unknown(1,2) = -t123 * t56;
	unknown(1,3) = -t133 * t56;
	unknown(1,4) = t140;
	unknown(1,5) = -t213 * t56;
	unknown(2,1) = t217;
	unknown(2,2) = -t252 * t56;
	unknown(2,3) = -t256 * t56;
	unknown(2,4) = t260 * t59 + t262;
	unknown(2,5) = -t300 * t56;
	unknown(3,1) = 0.0e0;
	unknown(3,2) = -t340 * t56;
	unknown(3,3) = -t350 * t56;
	unknown(3,4) = t354 * t59;
	unknown(3,5) = -t394 * t56;
	unknown(4,1) = t55 * t59 - t262;
	unknown(4,2) = t123 * t59;
	unknown(4,3) = t133 * t59;
	unknown(4,4) = -t217;
	unknown(4,5) = t213 * t59;
	unknown(5,1) = t140;
	unknown(5,2) = t252 * t59;
	unknown(5,3) = t256 * t59;
	unknown(5,4) = t260 * t56 - t60;
	unknown(5,5) = t300 * t59;
	unknown(6,1) = 0.0e0;
	unknown(6,2) = t340 * t59;
	unknown(6,3) = t350 * t59;
	unknown(6,4) = t354 * t56;
	unknown(6,5) = t394 * t59;
	unknown(7,1) = t47 * t53 / 0.4e1 + t52 * t48 / 0.4e1;
	unknown(7,2) = t107 * t53 + t121 * t48;
	unknown(7,3) = t137;
	unknown(7,4) = 0.0e0;
	unknown(7,5) = t196 * t53 + t211 * t48;
	unknown(8,1) = -t133;
	unknown(8,2) = t236 * t53 + t250 * t48;
	unknown(8,3) = t260;
	unknown(8,4) = 0.0e0;
	unknown(8,5) = t283 * t53 + t298 * t48;
	unknown(9,1) = 0.0e0;
	unknown(9,2) = t323 * t53 + t338 * t48;
	unknown(9,3) = t354;
	unknown(9,4) = 0.0e0;
	unknown(9,5) = t376 * t53 + t392 * t48;
	JR_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
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