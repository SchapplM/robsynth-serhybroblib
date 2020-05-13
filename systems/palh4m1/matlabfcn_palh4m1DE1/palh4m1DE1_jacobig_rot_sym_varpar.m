% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh4m1DE1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
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
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 22:26
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = palh4m1DE1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1DE1_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1DE1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1DE1_jacobig_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->15)
	unknown=NaN(3,5);
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
	Jg_rot = unknown;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->15)
	unknown=NaN(3,5);
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
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	Jg_rot = unknown;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:50
	% EndTime: 2020-04-11 22:24:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1362->50), mult. (1416->93), div. (56->8), fcn. (386->6), ass. (0->59)
	unknown=NaN(3,5);
	t1 = sin(qJ(5));
	t2 = pkin(1) * t1;
	t3 = -t2 + pkin(2);
	t5 = 0.2e1 * pkin(2) * t2;
	t6 = pkin(1) ^ 2;
	t10 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t14 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t16 = sqrt(-t14 * t10);
	t17 = 0.1e1 / t16;
	t18 = t17 * t3;
	t23 = -0.2e1 * t14 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t10;
	t26 = cos(qJ(5));
	t27 = t26 * pkin(1);
	t28 = qJ(2) + pkin(6);
	t31 = 0.1e1 / t28;
	t33 = pkin(2) ^ 2;
	t34 = -t5 + t6 + t33;
	t35 = 0.1e1 / t34;
	t38 = t28 ^ 2;
	t39 = pkin(3) ^ 2;
	t40 = -t5 + t6 + t33 + t38 - t39;
	t41 = t40 * t27;
	t42 = t16 * t3 + t41;
	t43 = 0.1e1 / t38;
	t48 = pkin(1) * t16 * t26;
	t50 = t40 * t3 - t48;
	t51 = 0.1e1 / t50;
	t54 = t42 ^ 2;
	t55 = t50 ^ 2;
	t56 = 0.1e1 / t55;
	t59 = 0.1e1 / (t56 * t54 + 0.1e1);
	t60 = t59 * t34 * t28;
	t62 = t17 * t26;
	t76 = t59 * t56 * t34;
	t78 = t60 * t51 * (t35 * t31 * (t23 * t18 / 0.2e1 + 0.2e1 * t28 * t27) - t35 * t43 * t42) - t76 * t28 * t42 * (t35 * t31 * (-t23 * pkin(1) * t62 / 0.2e1 + 0.2e1 * t28 * t3) - t35 * t43 * t50);
	t79 = sin(qJ(1));
	t84 = pkin(1) * pkin(2);
	t86 = t14 * pkin(2) * t27 + t84 * t26 * t10;
	t90 = t26 ^ 2;
	t98 = t34 ^ 2;
	t99 = 0.1e1 / t98;
	t101 = pkin(2) * t27;
	t126 = t60 * t51 * (t35 * t31 * (-0.2e1 * pkin(2) * t90 * t6 + t86 * t18 - t40 * t2 - t48) + 0.2e1 * t101 * t99 * t31 * t42) - t76 * t28 * t42 * (t35 * t31 * (pkin(1) * t16 * t1 - t86 * pkin(1) * t62 - 0.2e1 * t84 * t26 * t3 - t41) + 0.2e1 * t101 * t99 * t31 * t50);
	t128 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t79 * t78);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = (t79 * t126);
	unknown(2,1) = 0;
	unknown(2,2) = -(t128 * t78);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = -(t128 * t126);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	Jg_rot = unknown;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:50
	% EndTime: 2020-04-11 22:24:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1362->50), mult. (1416->93), div. (56->8), fcn. (386->6), ass. (0->59)
	unknown=NaN(3,5);
	t1 = sin(qJ(5));
	t2 = pkin(1) * t1;
	t3 = -t2 + pkin(2);
	t5 = 0.2e1 * pkin(2) * t2;
	t6 = pkin(1) ^ 2;
	t10 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t14 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t16 = sqrt(-t14 * t10);
	t17 = 0.1e1 / t16;
	t18 = t17 * t3;
	t23 = -0.2e1 * t14 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t10;
	t26 = cos(qJ(5));
	t27 = t26 * pkin(1);
	t28 = qJ(2) + pkin(6);
	t31 = 0.1e1 / t28;
	t33 = pkin(2) ^ 2;
	t34 = -t5 + t6 + t33;
	t35 = 0.1e1 / t34;
	t38 = t28 ^ 2;
	t39 = pkin(3) ^ 2;
	t40 = -t5 + t6 + t33 + t38 - t39;
	t41 = t40 * t27;
	t42 = t16 * t3 + t41;
	t43 = 0.1e1 / t38;
	t48 = pkin(1) * t16 * t26;
	t50 = t40 * t3 - t48;
	t51 = 0.1e1 / t50;
	t54 = t42 ^ 2;
	t55 = t50 ^ 2;
	t56 = 0.1e1 / t55;
	t59 = 0.1e1 / (t56 * t54 + 0.1e1);
	t60 = t59 * t34 * t28;
	t62 = t17 * t26;
	t76 = t59 * t56 * t34;
	t78 = t60 * t51 * (t35 * t31 * (t23 * t18 / 0.2e1 + 0.2e1 * t28 * t27) - t35 * t43 * t42) - t76 * t28 * t42 * (t35 * t31 * (-t23 * pkin(1) * t62 / 0.2e1 + 0.2e1 * t28 * t3) - t35 * t43 * t50);
	t79 = sin(qJ(1));
	t84 = pkin(1) * pkin(2);
	t86 = t14 * pkin(2) * t27 + t84 * t26 * t10;
	t90 = t26 ^ 2;
	t98 = t34 ^ 2;
	t99 = 0.1e1 / t98;
	t101 = pkin(2) * t27;
	t126 = t60 * t51 * (t35 * t31 * (-0.2e1 * pkin(2) * t90 * t6 + t86 * t18 - t40 * t2 - t48) + 0.2e1 * t101 * t99 * t31 * t42) - t76 * t28 * t42 * (t35 * t31 * (pkin(1) * t16 * t1 - t86 * pkin(1) * t62 - 0.2e1 * t84 * t26 * t3 - t41) + 0.2e1 * t101 * t99 * t31 * t50);
	t128 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t79 * t78);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = (t79 * t126);
	unknown(2,1) = 0;
	unknown(2,2) = -(t128 * t78);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = -(t128 * t126);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	Jg_rot = unknown;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:50
	% EndTime: 2020-04-11 22:24:50
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1858->61), mult. (1872->124), div. (96->12), fcn. (478->6), ass. (0->69)
	unknown=NaN(3,5);
	t1 = sin(qJ(5));
	t2 = pkin(1) * t1;
	t3 = -t2 + pkin(2);
	t5 = 0.2e1 * pkin(2) * t2;
	t6 = pkin(1) ^ 2;
	t10 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t14 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t15 = t14 * t10;
	t16 = sqrt(-t15);
	t17 = 0.1e1 / t16;
	t18 = t17 * t3;
	t23 = -0.2e1 * t14 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t10;
	t26 = cos(qJ(5));
	t27 = t26 * pkin(1);
	t28 = qJ(2) + pkin(6);
	t31 = 0.1e1 / t28;
	t33 = pkin(2) ^ 2;
	t34 = -t5 + t6 + t33;
	t35 = 0.1e1 / t34;
	t38 = t28 ^ 2;
	t39 = pkin(3) ^ 2;
	t40 = -t5 + t6 + t33 + t38 - t39;
	t41 = t40 * t27;
	t42 = t16 * t3 + t41;
	t43 = 0.1e1 / t38;
	t48 = pkin(1) * t16 * t26;
	t50 = t40 * t3 - t48;
	t51 = 0.1e1 / t50;
	t54 = t42 ^ 2;
	t55 = t50 ^ 2;
	t56 = 0.1e1 / t55;
	t59 = 0.1e1 / (t56 * t54 + 0.1e1);
	t60 = t59 * t34 * t28;
	t62 = t17 * t26;
	t76 = t59 * t56 * t34;
	t78 = t60 * t51 * (t35 * t31 * (t23 * t18 / 0.2e1 + 0.2e1 * t28 * t27) - t35 * t43 * t42) - t76 * t28 * t42 * (t35 * t31 * (-t23 * pkin(1) * t62 / 0.2e1 + 0.2e1 * t28 * t3) - t35 * t43 * t50);
	t79 = sin(qJ(1));
	t81 = 0.1e1 / pkin(3);
	t89 = t5 - t6 - t33 + t38 + t39;
	t90 = 0.1e1 / t89;
	t93 = t89 ^ 2;
	t94 = 0.1e1 / t93;
	t97 = 0.1e1 / (-t94 * t15 + 0.1e1);
	t108 = t97 * t94 * t16;
	t110 = t97 * pkin(3) * t28 * t90 * (-t16 * t81 * t43 + t23 * t17 * t81 * t31 / 0.2e1) - t108 * pkin(3) * t28 * (0.2e1 * t81 * t31 * t28 - t81 * t43 * t89);
	t116 = pkin(1) * pkin(2);
	t118 = t14 * pkin(2) * t27 + t116 * t26 * t10;
	t122 = t26 ^ 2;
	t130 = t34 ^ 2;
	t131 = 0.1e1 / t130;
	t133 = pkin(2) * t27;
	t158 = t60 * t51 * (t35 * t31 * (-0.2e1 * pkin(2) * t122 * t6 + t118 * t18 - t40 * t2 - t48) + 0.2e1 * t133 * t131 * t31 * t42) - t76 * t28 * t42 * (t35 * t31 * (pkin(1) * t16 * t1 - t118 * pkin(1) * t62 - 0.2e1 * t116 * t26 * t3 - t41) + 0.2e1 * t133 * t131 * t31 * t50);
	t166 = t97 * t90 * t118 * t17 - 0.2e1 * t108 * t133;
	t169 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t79 * t110 + t79 * t78);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = (t79 * t158 + t79 * t166);
	unknown(2,1) = 0;
	unknown(2,2) = (-t169 * t110 - t169 * t78);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = (-t169 * t158 - t169 * t166);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	Jg_rot = unknown;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:51
	% EndTime: 2020-04-11 22:24:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1859->62), mult. (1872->124), div. (96->12), fcn. (480->6), ass. (0->69)
	unknown=NaN(3,5);
	t1 = sin(qJ(5));
	t2 = pkin(1) * t1;
	t3 = -t2 + pkin(2);
	t5 = 0.2e1 * pkin(2) * t2;
	t6 = pkin(1) ^ 2;
	t10 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t14 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t15 = t14 * t10;
	t16 = sqrt(-t15);
	t17 = 0.1e1 / t16;
	t18 = t17 * t3;
	t23 = -0.2e1 * t14 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t10;
	t26 = cos(qJ(5));
	t27 = t26 * pkin(1);
	t28 = qJ(2) + pkin(6);
	t31 = 0.1e1 / t28;
	t33 = pkin(2) ^ 2;
	t34 = -t5 + t6 + t33;
	t35 = 0.1e1 / t34;
	t38 = t28 ^ 2;
	t39 = pkin(3) ^ 2;
	t40 = -t5 + t6 + t33 + t38 - t39;
	t41 = t40 * t27;
	t42 = t16 * t3 + t41;
	t43 = 0.1e1 / t38;
	t48 = pkin(1) * t16 * t26;
	t50 = t40 * t3 - t48;
	t51 = 0.1e1 / t50;
	t54 = t42 ^ 2;
	t55 = t50 ^ 2;
	t56 = 0.1e1 / t55;
	t59 = 0.1e1 / (t56 * t54 + 0.1e1);
	t60 = t59 * t34 * t28;
	t62 = t17 * t26;
	t76 = t59 * t56 * t34;
	t78 = t60 * t51 * (t35 * t31 * (t23 * t18 / 0.2e1 + 0.2e1 * t28 * t27) - t35 * t43 * t42) - t76 * t28 * t42 * (t35 * t31 * (-t23 * pkin(1) * t62 / 0.2e1 + 0.2e1 * t28 * t3) - t35 * t43 * t50);
	t79 = sin(qJ(1));
	t81 = 0.1e1 / pkin(3);
	t89 = t5 - t6 - t33 + t38 + t39;
	t90 = 0.1e1 / t89;
	t93 = t89 ^ 2;
	t94 = 0.1e1 / t93;
	t97 = 0.1e1 / (-t94 * t15 + 0.1e1);
	t108 = t97 * t94 * t16;
	t110 = t97 * pkin(3) * t28 * t90 * (-t16 * t81 * t43 + t23 * t17 * t81 * t31 / 0.2e1) - t108 * pkin(3) * t28 * (0.2e1 * t81 * t31 * t28 - t81 * t43 * t89);
	t116 = pkin(1) * pkin(2);
	t118 = t14 * pkin(2) * t27 + t116 * t26 * t10;
	t122 = t26 ^ 2;
	t130 = t34 ^ 2;
	t131 = 0.1e1 / t130;
	t133 = pkin(2) * t27;
	t158 = t60 * t51 * (t35 * t31 * (-0.2e1 * pkin(2) * t122 * t6 + t118 * t18 - t40 * t2 - t48) + 0.2e1 * t133 * t131 * t31 * t42) - t76 * t28 * t42 * (t35 * t31 * (pkin(1) * t16 * t1 - t118 * pkin(1) * t62 - 0.2e1 * t116 * t26 * t3 - t41) + 0.2e1 * t133 * t131 * t31 * t50);
	t166 = t97 * t90 * t118 * t17 - 0.2e1 * t108 * t133;
	t169 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t79 * t110 + t79 * t78);
	unknown(1,3) = t79;
	unknown(1,4) = 0;
	unknown(1,5) = (t79 * t158 + t79 * t166);
	unknown(2,1) = 0;
	unknown(2,2) = (-t169 * t110 - t169 * t78);
	unknown(2,3) = -t169;
	unknown(2,4) = 0;
	unknown(2,5) = (-t169 * t158 - t169 * t166);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	Jg_rot = unknown;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:55
	% EndTime: 2020-04-11 22:24:55
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (3272->73), mult. (3428->167), div. (228->15), fcn. (896->12), ass. (0->85)
	unknown=NaN(3,5);
	t1 = sin(qJ(5));
	t2 = pkin(1) * t1;
	t3 = -t2 + pkin(2);
	t5 = 0.2e1 * pkin(2) * t2;
	t6 = pkin(1) ^ 2;
	t10 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t14 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t15 = t14 * t10;
	t16 = sqrt(-t15);
	t17 = 0.1e1 / t16;
	t18 = t17 * t3;
	t23 = -0.2e1 * t14 * (-qJ(2) - pkin(6) - pkin(3)) - 0.2e1 * (-qJ(2) - pkin(6) + pkin(3)) * t10;
	t26 = cos(qJ(5));
	t27 = t26 * pkin(1);
	t28 = qJ(2) + pkin(6);
	t31 = 0.1e1 / t28;
	t33 = pkin(2) ^ 2;
	t34 = -t5 + t6 + t33;
	t35 = 0.1e1 / t34;
	t38 = t28 ^ 2;
	t39 = pkin(3) ^ 2;
	t40 = -t5 + t6 + t33 + t38 - t39;
	t41 = t40 * t27;
	t42 = t16 * t3 + t41;
	t43 = 0.1e1 / t38;
	t45 = t35 * t43 * t42;
	t48 = pkin(1) * t16 * t26;
	t50 = t40 * t3 - t48;
	t51 = 0.1e1 / t50;
	t54 = t42 ^ 2;
	t55 = t50 ^ 2;
	t56 = 0.1e1 / t55;
	t59 = 0.1e1 / (t56 * t54 + 0.1e1);
	t60 = t59 * t34 * t28;
	t62 = t17 * t26;
	t71 = t35 * t43 * t50;
	t76 = t59 * t56 * t34;
	t78 = t60 * t51 * (t35 * t31 * (t23 * t18 / 0.2e1 + 0.2e1 * t28 * t27) - t45) - t76 * t28 * t42 * (t35 * t31 * (-t23 * pkin(1) * t62 / 0.2e1 + 0.2e1 * t28 * t3) - t71);
	t79 = sin(qJ(1));
	t81 = 0.1e1 / pkin(3);
	t89 = t5 - t6 - t33 + t38 + t39;
	t90 = 0.1e1 / t89;
	t93 = t89 ^ 2;
	t94 = 0.1e1 / t93;
	t97 = 0.1e1 / (-t94 * t15 + 0.1e1);
	t108 = t97 * t94 * t16;
	t110 = t97 * pkin(3) * t28 * t90 * (-t16 * t81 * t43 + t23 * t17 * t81 * t31 / 0.2e1) - t108 * pkin(3) * t28 * (0.2e1 * t81 * t31 * t28 - t81 * t43 * t89);
	t113 = cos(qJ(1));
	t115 = t35 * t43;
	t116 = t115 * t42 * t113;
	t118 = t34 ^ 2;
	t119 = 0.1e1 / t118;
	t124 = sqrt(t119 * t43 * t54 + t119 * t43 * t55);
	t125 = 0.1e1 / t124;
	t127 = 0.1e1 / t39;
	t133 = sqrt(-t15 * t127 * t43 + t127 * t43 * t93);
	t134 = 0.1e1 / t133;
	t136 = t134 * t16 * t81 * t125;
	t139 = t115 * t50 * t113;
	t142 = t134 * t81 * t89 * t125;
	t145 = cos(qJ(3));
	t150 = sin(qJ(3));
	t156 = pkin(1) * pkin(2);
	t158 = t14 * pkin(2) * t27 + t156 * t26 * t10;
	t162 = t26 ^ 2;
	t171 = pkin(2) * t27;
	t196 = t60 * t51 * (t35 * t31 * (-0.2e1 * pkin(2) * t162 * t6 + t158 * t18 - t40 * t2 - t48) + 0.2e1 * t171 * t119 * t31 * t42) - t76 * t28 * t42 * (t35 * t31 * (pkin(1) * t16 * t1 - t158 * pkin(1) * t62 - 0.2e1 * t156 * t26 * t3 - t41) + 0.2e1 * t171 * t119 * t31 * t50);
	t204 = t97 * t90 * t158 * t17 - 0.2e1 * t108 * t171;
	t211 = t115 * t42 * t79;
	t214 = t115 * t50 * t79;
	unknown(1,1) = 0;
	unknown(1,2) = (t79 * t110 + t79 * t78);
	unknown(1,3) = t79;
	unknown(1,4) = (t145 * (-t136 * t116 + t142 * t139) + t150 * (-t142 * t116 - t136 * t139));
	unknown(1,5) = (t79 * t196 + t79 * t204);
	unknown(2,1) = 0;
	unknown(2,2) = (-t113 * t110 - t113 * t78);
	unknown(2,3) = -t113;
	unknown(2,4) = (t145 * (-t136 * t211 + t142 * t214) + t150 * (-t136 * t214 - t142 * t211));
	unknown(2,5) = (-t113 * t196 - t113 * t204);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = (t145 * (t136 * t71 + t142 * t45) + t150 * (-t136 * t45 + t142 * t71));
	unknown(3,5) = 0;
	Jg_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:24:49
	% EndTime: 2020-04-11 22:24:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->17)
	unknown=NaN(3,5);
	t1 = sin(qJ(1));
	t2 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = 0;
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = t1;
	unknown(2,1) = 0;
	unknown(2,2) = 0;
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = -t2;
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = 0;
	unknown(3,5) = 0;
	Jg_rot = unknown;
else
	Jg_rot=NaN(3,5);
end