% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% palh4m1TE
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
% Datum: 2020-04-11 21:48
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = palh4m1TE_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh4m1TE_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh4m1TE_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'palh4m1TE_jacobig_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
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
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
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
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1362->50), mult. (1416->93), div. (56->8), fcn. (386->6), ass. (0->59)
	unknown=NaN(3,5);
	t1 = sin(qJ(5));
	t2 = t1 * pkin(1);
	t3 = -t2 + pkin(2);
	t5 = 0.2e1 * t2 * pkin(2);
	t6 = pkin(1) ^ 2;
	t10 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t14 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t16 = sqrt(-t10 * t14);
	t17 = 0.1e1 / t16;
	t18 = t3 * t17;
	t23 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t14 - 0.2e1 * t10 * (-qJ(2) - pkin(6) + pkin(3));
	t26 = cos(qJ(5));
	t27 = pkin(1) * t26;
	t28 = qJ(2) + pkin(6);
	t31 = 0.1e1 / t28;
	t33 = pkin(2) ^ 2;
	t34 = -t5 + t6 + t33;
	t35 = 0.1e1 / t34;
	t38 = t28 ^ 2;
	t39 = pkin(3) ^ 2;
	t40 = -t5 + t6 + t33 + t38 - t39;
	t41 = t27 * t40;
	t42 = t3 * t16 + t41;
	t43 = 0.1e1 / t38;
	t48 = t26 * t16 * pkin(1);
	t50 = t3 * t40 - t48;
	t51 = 0.1e1 / t50;
	t54 = t42 ^ 2;
	t55 = t50 ^ 2;
	t56 = 0.1e1 / t55;
	t59 = 0.1e1 / (t54 * t56 + 0.1e1);
	t60 = t28 * t34 * t59;
	t62 = t26 * t17;
	t76 = t34 * t56 * t59;
	t78 = ((t18 * t23 / 0.2e1 + 0.2e1 * t27 * t28) * t31 * t35 - t42 * t43 * t35) * t51 * t60 - ((-t62 * pkin(1) * t23 / 0.2e1 + 0.2e1 * t3 * t28) * t31 * t35 - t50 * t43 * t35) * t42 * t28 * t76;
	t79 = sin(qJ(1));
	t84 = pkin(1) * pkin(2);
	t86 = t27 * pkin(2) * t14 + t10 * t26 * t84;
	t90 = t26 ^ 2;
	t98 = t34 ^ 2;
	t99 = 0.1e1 / t98;
	t101 = t27 * pkin(2);
	t126 = ((-0.2e1 * t6 * t90 * pkin(2) + t18 * t86 - t2 * t40 - t48) * t31 * t35 + 0.2e1 * t42 * t31 * t99 * t101) * t51 * t60 - ((t1 * t16 * pkin(1) - t62 * pkin(1) * t86 - 0.2e1 * t3 * t26 * t84 - t41) * t31 * t35 + 0.2e1 * t50 * t31 * t99 * t101) * t42 * t28 * t76;
	t128 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t78 * t79);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = (t126 * t79);
	unknown(2,1) = 0;
	unknown(2,2) = -(t78 * t128);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = -(t126 * t128);
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
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1362->50), mult. (1416->93), div. (56->8), fcn. (386->6), ass. (0->59)
	unknown=NaN(3,5);
	t1 = sin(qJ(5));
	t2 = t1 * pkin(1);
	t3 = -t2 + pkin(2);
	t5 = 0.2e1 * t2 * pkin(2);
	t6 = pkin(1) ^ 2;
	t10 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t14 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t16 = sqrt(-t10 * t14);
	t17 = 0.1e1 / t16;
	t18 = t3 * t17;
	t23 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t14 - 0.2e1 * t10 * (-qJ(2) - pkin(6) + pkin(3));
	t26 = cos(qJ(5));
	t27 = pkin(1) * t26;
	t28 = qJ(2) + pkin(6);
	t31 = 0.1e1 / t28;
	t33 = pkin(2) ^ 2;
	t34 = -t5 + t6 + t33;
	t35 = 0.1e1 / t34;
	t38 = t28 ^ 2;
	t39 = pkin(3) ^ 2;
	t40 = -t5 + t6 + t33 + t38 - t39;
	t41 = t27 * t40;
	t42 = t3 * t16 + t41;
	t43 = 0.1e1 / t38;
	t48 = t26 * t16 * pkin(1);
	t50 = t3 * t40 - t48;
	t51 = 0.1e1 / t50;
	t54 = t42 ^ 2;
	t55 = t50 ^ 2;
	t56 = 0.1e1 / t55;
	t59 = 0.1e1 / (t54 * t56 + 0.1e1);
	t60 = t28 * t34 * t59;
	t62 = t26 * t17;
	t76 = t34 * t56 * t59;
	t78 = ((t18 * t23 / 0.2e1 + 0.2e1 * t27 * t28) * t31 * t35 - t42 * t43 * t35) * t51 * t60 - ((-t62 * pkin(1) * t23 / 0.2e1 + 0.2e1 * t3 * t28) * t31 * t35 - t50 * t43 * t35) * t42 * t28 * t76;
	t79 = sin(qJ(1));
	t84 = pkin(1) * pkin(2);
	t86 = t27 * pkin(2) * t14 + t10 * t26 * t84;
	t90 = t26 ^ 2;
	t98 = t34 ^ 2;
	t99 = 0.1e1 / t98;
	t101 = t27 * pkin(2);
	t126 = ((-0.2e1 * t6 * t90 * pkin(2) + t18 * t86 - t2 * t40 - t48) * t31 * t35 + 0.2e1 * t42 * t31 * t99 * t101) * t51 * t60 - ((t1 * t16 * pkin(1) - t62 * pkin(1) * t86 - 0.2e1 * t3 * t26 * t84 - t41) * t31 * t35 + 0.2e1 * t50 * t31 * t99 * t101) * t42 * t28 * t76;
	t128 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t78 * t79);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = (t126 * t79);
	unknown(2,1) = 0;
	unknown(2,2) = -(t78 * t128);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = -(t126 * t128);
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
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (1858->61), mult. (1872->124), div. (96->12), fcn. (478->6), ass. (0->69)
	unknown=NaN(3,5);
	t1 = sin(qJ(5));
	t2 = t1 * pkin(1);
	t3 = -t2 + pkin(2);
	t5 = 0.2e1 * t2 * pkin(2);
	t6 = pkin(1) ^ 2;
	t10 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t14 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t15 = t10 * t14;
	t16 = sqrt(-t15);
	t17 = 0.1e1 / t16;
	t18 = t3 * t17;
	t23 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t14 - 0.2e1 * t10 * (-qJ(2) - pkin(6) + pkin(3));
	t26 = cos(qJ(5));
	t27 = pkin(1) * t26;
	t28 = qJ(2) + pkin(6);
	t31 = 0.1e1 / t28;
	t33 = pkin(2) ^ 2;
	t34 = -t5 + t6 + t33;
	t35 = 0.1e1 / t34;
	t38 = t28 ^ 2;
	t39 = pkin(3) ^ 2;
	t40 = -t5 + t6 + t33 + t38 - t39;
	t41 = t27 * t40;
	t42 = t3 * t16 + t41;
	t43 = 0.1e1 / t38;
	t48 = t26 * t16 * pkin(1);
	t50 = t3 * t40 - t48;
	t51 = 0.1e1 / t50;
	t54 = t42 ^ 2;
	t55 = t50 ^ 2;
	t56 = 0.1e1 / t55;
	t59 = 0.1e1 / (t54 * t56 + 0.1e1);
	t60 = t28 * t34 * t59;
	t62 = t26 * t17;
	t76 = t34 * t56 * t59;
	t78 = ((t18 * t23 / 0.2e1 + 0.2e1 * t27 * t28) * t31 * t35 - t42 * t43 * t35) * t51 * t60 - ((-t62 * pkin(1) * t23 / 0.2e1 + 0.2e1 * t3 * t28) * t31 * t35 - t50 * t43 * t35) * t42 * t28 * t76;
	t79 = sin(qJ(1));
	t81 = 0.1e1 / pkin(3);
	t89 = t5 - t6 - t33 + t38 + t39;
	t90 = 0.1e1 / t89;
	t93 = t89 ^ 2;
	t94 = 0.1e1 / t93;
	t97 = 0.1e1 / (-t15 * t94 + 0.1e1);
	t108 = t16 * t94 * t97;
	t110 = (-t43 * t81 * t16 + t31 * t81 * t17 * t23 / 0.2e1) * t90 * t28 * pkin(3) * t97 - (0.2e1 * t28 * t31 * t81 - t89 * t43 * t81) * t28 * pkin(3) * t108;
	t116 = pkin(1) * pkin(2);
	t118 = t27 * pkin(2) * t14 + t10 * t26 * t116;
	t122 = t26 ^ 2;
	t130 = t34 ^ 2;
	t131 = 0.1e1 / t130;
	t133 = t27 * pkin(2);
	t158 = ((-0.2e1 * t6 * t122 * pkin(2) + t18 * t118 - t2 * t40 - t48) * t31 * t35 + 0.2e1 * t42 * t31 * t131 * t133) * t51 * t60 - ((t1 * t16 * pkin(1) - t62 * pkin(1) * t118 - 0.2e1 * t3 * t26 * t116 - t41) * t31 * t35 + 0.2e1 * t50 * t31 * t131 * t133) * t42 * t28 * t76;
	t166 = t17 * t118 * t90 * t97 - 0.2e1 * t133 * t108;
	t169 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t110 * t79 + t78 * t79);
	unknown(1,3) = 0;
	unknown(1,4) = 0;
	unknown(1,5) = (t158 * t79 + t166 * t79);
	unknown(2,1) = 0;
	unknown(2,2) = (-t110 * t169 - t78 * t169);
	unknown(2,3) = 0;
	unknown(2,4) = 0;
	unknown(2,5) = (-t158 * t169 - t166 * t169);
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
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (1859->62), mult. (1872->124), div. (96->12), fcn. (480->6), ass. (0->69)
	unknown=NaN(3,5);
	t1 = sin(qJ(5));
	t2 = t1 * pkin(1);
	t3 = -t2 + pkin(2);
	t5 = 0.2e1 * t2 * pkin(2);
	t6 = pkin(1) ^ 2;
	t10 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t14 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t15 = t10 * t14;
	t16 = sqrt(-t15);
	t17 = 0.1e1 / t16;
	t18 = t3 * t17;
	t23 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t14 - 0.2e1 * t10 * (-qJ(2) - pkin(6) + pkin(3));
	t26 = cos(qJ(5));
	t27 = pkin(1) * t26;
	t28 = qJ(2) + pkin(6);
	t31 = 0.1e1 / t28;
	t33 = pkin(2) ^ 2;
	t34 = -t5 + t6 + t33;
	t35 = 0.1e1 / t34;
	t38 = t28 ^ 2;
	t39 = pkin(3) ^ 2;
	t40 = -t5 + t6 + t33 + t38 - t39;
	t41 = t27 * t40;
	t42 = t3 * t16 + t41;
	t43 = 0.1e1 / t38;
	t48 = t26 * t16 * pkin(1);
	t50 = t3 * t40 - t48;
	t51 = 0.1e1 / t50;
	t54 = t42 ^ 2;
	t55 = t50 ^ 2;
	t56 = 0.1e1 / t55;
	t59 = 0.1e1 / (t54 * t56 + 0.1e1);
	t60 = t28 * t34 * t59;
	t62 = t26 * t17;
	t76 = t34 * t56 * t59;
	t78 = ((t18 * t23 / 0.2e1 + 0.2e1 * t27 * t28) * t31 * t35 - t42 * t43 * t35) * t51 * t60 - ((-t62 * pkin(1) * t23 / 0.2e1 + 0.2e1 * t3 * t28) * t31 * t35 - t50 * t43 * t35) * t42 * t28 * t76;
	t79 = sin(qJ(1));
	t81 = 0.1e1 / pkin(3);
	t89 = t5 - t6 - t33 + t38 + t39;
	t90 = 0.1e1 / t89;
	t93 = t89 ^ 2;
	t94 = 0.1e1 / t93;
	t97 = 0.1e1 / (-t15 * t94 + 0.1e1);
	t108 = t16 * t94 * t97;
	t110 = (-t43 * t81 * t16 + t31 * t81 * t17 * t23 / 0.2e1) * t90 * t28 * pkin(3) * t97 - (0.2e1 * t28 * t31 * t81 - t89 * t43 * t81) * t28 * pkin(3) * t108;
	t116 = pkin(1) * pkin(2);
	t118 = t27 * pkin(2) * t14 + t10 * t26 * t116;
	t122 = t26 ^ 2;
	t130 = t34 ^ 2;
	t131 = 0.1e1 / t130;
	t133 = t27 * pkin(2);
	t158 = ((-0.2e1 * t6 * t122 * pkin(2) + t18 * t118 - t2 * t40 - t48) * t31 * t35 + 0.2e1 * t42 * t31 * t131 * t133) * t51 * t60 - ((t1 * t16 * pkin(1) - t62 * pkin(1) * t118 - 0.2e1 * t3 * t26 * t116 - t41) * t31 * t35 + 0.2e1 * t50 * t31 * t131 * t133) * t42 * t28 * t76;
	t166 = t17 * t118 * t90 * t97 - 0.2e1 * t133 * t108;
	t169 = cos(qJ(1));
	unknown(1,1) = 0;
	unknown(1,2) = (t110 * t79 + t78 * t79);
	unknown(1,3) = t79;
	unknown(1,4) = 0;
	unknown(1,5) = (t158 * t79 + t166 * t79);
	unknown(2,1) = 0;
	unknown(2,2) = (-t110 * t169 - t78 * t169);
	unknown(2,3) = -t169;
	unknown(2,4) = 0;
	unknown(2,5) = (-t158 * t169 - t166 * t169);
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
	% StartTime: 2020-04-11 21:48:21
	% EndTime: 2020-04-11 21:48:22
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (2324->71), mult. (2372->160), div. (132->12), fcn. (620->8), ass. (0->79)
	unknown=NaN(3,5);
	t1 = sin(qJ(5));
	t2 = t1 * pkin(1);
	t3 = -t2 + pkin(2);
	t5 = 0.2e1 * t2 * pkin(2);
	t6 = pkin(1) ^ 2;
	t10 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) + pkin(3)) * (pkin(2) - qJ(2) - pkin(6) - pkin(3));
	t14 = -t5 + t6 + (pkin(2) + qJ(2) + pkin(6) - pkin(3)) * (pkin(2) - qJ(2) - pkin(6) + pkin(3));
	t15 = t10 * t14;
	t16 = sqrt(-t15);
	t17 = 0.1e1 / t16;
	t18 = t3 * t17;
	t23 = -0.2e1 * (-qJ(2) - pkin(6) - pkin(3)) * t14 - 0.2e1 * t10 * (-qJ(2) - pkin(6) + pkin(3));
	t26 = cos(qJ(5));
	t27 = pkin(1) * t26;
	t28 = qJ(2) + pkin(6);
	t31 = 0.1e1 / t28;
	t33 = pkin(2) ^ 2;
	t34 = -t5 + t6 + t33;
	t35 = 0.1e1 / t34;
	t38 = t28 ^ 2;
	t39 = pkin(3) ^ 2;
	t40 = -t5 + t6 + t33 + t38 - t39;
	t41 = t27 * t40;
	t42 = t3 * t16 + t41;
	t43 = 0.1e1 / t38;
	t44 = t42 * t43;
	t48 = t26 * t16 * pkin(1);
	t50 = t3 * t40 - t48;
	t51 = 0.1e1 / t50;
	t54 = t42 ^ 2;
	t55 = t50 ^ 2;
	t56 = 0.1e1 / t55;
	t59 = 0.1e1 / (t54 * t56 + 0.1e1);
	t60 = t28 * t34 * t59;
	t62 = t26 * t17;
	t70 = t50 * t43;
	t76 = t34 * t56 * t59;
	t78 = ((t18 * t23 / 0.2e1 + 0.2e1 * t27 * t28) * t31 * t35 - t44 * t35) * t51 * t60 - ((-t62 * pkin(1) * t23 / 0.2e1 + 0.2e1 * t3 * t28) * t31 * t35 - t70 * t35) * t42 * t28 * t76;
	t79 = sin(qJ(1));
	t81 = 0.1e1 / pkin(3);
	t89 = t5 - t6 - t33 + t38 + t39;
	t90 = 0.1e1 / t89;
	t93 = t89 ^ 2;
	t94 = 0.1e1 / t93;
	t97 = 0.1e1 / (-t15 * t94 + 0.1e1);
	t108 = t16 * t94 * t97;
	t110 = (-t43 * t81 * t16 + t31 * t81 * t17 * t23 / 0.2e1) * t90 * t28 * pkin(3) * t97 - (0.2e1 * t28 * t31 * t81 - t89 * t43 * t81) * t28 * pkin(3) * t108;
	t113 = cos(qJ(1));
	t115 = t113 * t42 * t43;
	t117 = t35 * t81 * t16;
	t120 = t113 * t50 * t43;
	t122 = t35 * t89 * t81;
	t125 = cos(qJ(3));
	t130 = sin(qJ(3));
	t136 = pkin(1) * pkin(2);
	t138 = t27 * pkin(2) * t14 + t10 * t26 * t136;
	t142 = t26 ^ 2;
	t150 = t34 ^ 2;
	t151 = 0.1e1 / t150;
	t153 = t27 * pkin(2);
	t178 = ((-0.2e1 * t6 * t142 * pkin(2) + t18 * t138 - t2 * t40 - t48) * t31 * t35 + 0.2e1 * t42 * t31 * t151 * t153) * t51 * t60 - ((t1 * t16 * pkin(1) - t62 * pkin(1) * t138 - 0.2e1 * t3 * t26 * t136 - t41) * t31 * t35 + 0.2e1 * t50 * t31 * t151 * t153) * t42 * t28 * t76;
	t186 = t17 * t138 * t90 * t97 - 0.2e1 * t153 * t108;
	t193 = t79 * t42 * t43;
	t196 = t79 * t50 * t43;
	unknown(1,1) = 0;
	unknown(1,2) = (t110 * t79 + t78 * t79);
	unknown(1,3) = t79;
	unknown(1,4) = ((-t115 * t117 + t120 * t122) * t125 / 0.4e1 + (-t115 * t122 - t120 * t117) * t130 / 0.4e1);
	unknown(1,5) = (t178 * t79 + t186 * t79);
	unknown(2,1) = 0;
	unknown(2,2) = (-t110 * t113 - t78 * t113);
	unknown(2,3) = -t113;
	unknown(2,4) = ((-t193 * t117 + t196 * t122) * t125 / 0.4e1 + (-t196 * t117 - t193 * t122) * t130 / 0.4e1);
	unknown(2,5) = (-t178 * t113 - t186 * t113);
	unknown(3,1) = 1;
	unknown(3,2) = 0;
	unknown(3,3) = 0;
	unknown(3,4) = ((t70 * t117 + t44 * t122) * t125 / 0.4e1 + (-t44 * t117 + t70 * t122) * t130 / 0.4e1);
	unknown(3,5) = 0;
	Jg_rot = unknown;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobig_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 1
	% StartTime: 2020-04-11 21:48:18
	% EndTime: 2020-04-11 21:48:18
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