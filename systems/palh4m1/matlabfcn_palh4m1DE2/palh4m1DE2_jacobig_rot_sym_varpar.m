% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh4m1DE2
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
% Datum: 2020-04-11 22:54
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = palh4m1DE2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1DE2_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1DE2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1DE2_jacobig_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
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
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
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
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
	% DurationCPUTime: 0.06s
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
	% StartTime: 2020-04-11 22:53:33
	% EndTime: 2020-04-11 22:53:33
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
	% StartTime: 2020-04-11 22:53:33
	% EndTime: 2020-04-11 22:53:33
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
	% StartTime: 2020-04-11 22:53:33
	% EndTime: 2020-04-11 22:53:33
	% DurationCPUTime: 0.07s
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
	% StartTime: 2020-04-11 22:53:35
	% EndTime: 2020-04-11 22:53:35
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2504->67), mult. (2548->148), div. (156->13), fcn. (682->11), ass. (0->80)
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
	t84 = t81 * t31;
	t89 = t5 - t6 - t33 + t38 + t39;
	t90 = 0.1e1 / t89;
	t93 = t89 ^ 2;
	t94 = 0.1e1 / t93;
	t97 = 0.1e1 / (-t94 * t15 + 0.1e1);
	t108 = t97 * t94 * t16;
	t110 = t97 * pkin(3) * t28 * t90 * (-t16 * t81 * t43 + t23 * t17 * t84 / 0.2e1) - t108 * pkin(3) * t28 * (0.2e1 * t81 * t31 * t28 - t81 * t43 * t89);
	t113 = cos(qJ(1));
	t117 = t34 ^ 2;
	t118 = 0.1e1 / t117;
	t123 = sqrt(t118 * t43 * t54 + t118 * t43 * t55);
	t125 = 0.1e1 / t123 * t35;
	t129 = atan2(t16 * t84, t81 * t31 * t89);
	t130 = t129 + qJ(3);
	t131 = sin(t130);
	t132 = t131 * t125;
	t136 = cos(t130);
	t137 = t136 * t125;
	t143 = pkin(1) * pkin(2);
	t145 = t14 * pkin(2) * t27 + t143 * t26 * t10;
	t149 = t26 ^ 2;
	t156 = t31 * t42;
	t158 = pkin(2) * t27;
	t175 = t31 * t50;
	t183 = t60 * t51 * (t35 * t31 * (-0.2e1 * pkin(2) * t149 * t6 + t145 * t18 - t40 * t2 - t48) + 0.2e1 * t158 * t118 * t156) - t76 * t28 * t42 * (t35 * t31 * (pkin(1) * t16 * t1 - t145 * pkin(1) * t62 - 0.2e1 * t143 * t26 * t3 - t41) + 0.2e1 * t158 * t118 * t175);
	t191 = t97 * t90 * t145 * t17 - 0.2e1 * t108 * t158;
	unknown(1,1) = 0;
	unknown(1,2) = (t79 * t110 + t79 * t78);
	unknown(1,3) = t79;
	unknown(1,4) = (-t132 * t31 * t42 * t113 + t137 * t31 * t50 * t113);
	unknown(1,5) = (t79 * t183 + t79 * t191);
	unknown(2,1) = 0;
	unknown(2,2) = (-t113 * t110 - t113 * t78);
	unknown(2,3) = -t113;
	unknown(2,4) = (-t132 * t31 * t42 * t79 + t137 * t31 * t50 * t79);
	unknown(2,5) = (-t113 * t183 - t113 * t191);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = (t132 * t175 + t137 * t156);
	unknown(3,5) = 0;
	Jg_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 22:53:32
	% EndTime: 2020-04-11 22:53:32
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