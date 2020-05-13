% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% JR_rot [9x10]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = palh3m2OL_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),uint8(0),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_jacobiR_rot_sym_varpar: qJ has to be [10x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2OL_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_jacobiR_rot_sym_varpar: pkin has to be [16x1] (double)');
JR_rot=NaN(9,10);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:50
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t8 = sin(qJ(2));
	t9 = sin(qJ(1));
	t15 = t9 * t8;
	t11 = cos(qJ(1));
	t14 = t11 * t8;
	t10 = cos(qJ(2));
	t13 = t9 * t10;
	t12 = t11 * t10;
	t1 = [-t13, -t14, 0, 0, 0, 0, 0, 0, 0, 0; t12, -t15, 0, 0, 0, 0, 0, 0, 0, 0; 0, t10, 0, 0, 0, 0, 0, 0, 0, 0; t15, -t12, 0, 0, 0, 0, 0, 0, 0, 0; -t14, -t13, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:50
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->5), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t34 = cos(qJ(1));
	t33 = sin(qJ(1));
	t32 = qJ(2) + qJ(3);
	t31 = cos(t32);
	t30 = sin(t32);
	t29 = t34 * t31;
	t28 = t34 * t30;
	t27 = t33 * t31;
	t26 = t33 * t30;
	t1 = [t27, t28, t28, 0, 0, 0, 0, 0, 0, 0; -t29, t26, t26, 0, 0, 0, 0, 0, 0, 0; 0, -t31, -t31, 0, 0, 0, 0, 0, 0, 0; -t26, t29, t29, 0, 0, 0, 0, 0, 0, 0; t28, t27, t27, 0, 0, 0, 0, 0, 0, 0; 0, t30, t30, 0, 0, 0, 0, 0, 0, 0; t34, 0, 0, 0, 0, 0, 0, 0, 0, 0; t33, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (49->6), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t37 = qJ(2) + qJ(3) + qJ(4);
	t36 = cos(t37);
	t35 = sin(t37);
	t34 = t39 * t36;
	t33 = t39 * t35;
	t32 = t38 * t36;
	t31 = t38 * t35;
	t1 = [t32, t33, t33, t33, 0, 0, 0, 0, 0, 0; -t34, t31, t31, t31, 0, 0, 0, 0, 0, 0; 0, -t36, -t36, -t36, 0, 0, 0, 0, 0, 0; -t31, t34, t34, t34, 0, 0, 0, 0, 0, 0; t33, t32, t32, t32, 0, 0, 0, 0, 0, 0; 0, t35, t35, t35, 0, 0, 0, 0, 0, 0; t39, 0, 0, 0, 0, 0, 0, 0, 0, 0; t38, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:52
	% EndTime: 2020-05-07 04:44:52
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (105->26), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
	t100 = sin(qJ(1));
	t98 = qJ(2) + qJ(3) + qJ(4);
	t97 = cos(t98);
	t111 = t100 * t97;
	t99 = sin(qJ(5));
	t110 = t100 * t99;
	t102 = cos(qJ(1));
	t109 = t102 * t97;
	t108 = t102 * t99;
	t101 = cos(qJ(5));
	t107 = t97 * t101;
	t106 = t100 * t101;
	t105 = t102 * t101;
	t96 = sin(t98);
	t104 = t96 * t110;
	t103 = t96 * t108;
	t95 = t97 * t99;
	t94 = t96 * t105;
	t93 = t96 * t106;
	t92 = t97 * t105 - t110;
	t91 = t97 * t108 + t106;
	t90 = t97 * t106 + t108;
	t89 = t97 * t110 - t105;
	t1 = [t90, t94, t94, t94, t91, 0, 0, 0, 0, 0; -t92, t93, t93, t93, t89, 0, 0, 0, 0, 0; 0, -t107, -t107, -t107, t96 * t99, 0, 0, 0, 0, 0; -t89, -t103, -t103, -t103, t92, 0, 0, 0, 0, 0; t91, -t104, -t104, -t104, t90, 0, 0, 0, 0, 0; 0, t95, t95, t95, t96 * t101, 0, 0, 0, 0, 0; t100 * t96, -t109, -t109, -t109, 0, 0, 0, 0, 0, 0; -t102 * t96, -t111, -t111, -t111, 0, 0, 0, 0, 0, 0; 0, -t96, -t96, -t96, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t8 = sin(qJ(6));
	t9 = sin(qJ(1));
	t15 = t9 * t8;
	t11 = cos(qJ(1));
	t14 = t11 * t8;
	t10 = cos(qJ(6));
	t13 = t9 * t10;
	t12 = t11 * t10;
	t1 = [-t13, 0, 0, 0, 0, -t14, 0, 0, 0, 0; t12, 0, 0, 0, 0, -t15, 0, 0, 0, 0; 0, 0, 0, 0, 0, t10, 0, 0, 0, 0; t15, 0, 0, 0, 0, -t12, 0, 0, 0, 0; -t14, 0, 0, 0, 0, -t13, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t20 = qJ(2) + qJ(7);
	t18 = sin(t20);
	t21 = sin(qJ(1));
	t26 = t21 * t18;
	t19 = cos(t20);
	t25 = t21 * t19;
	t22 = cos(qJ(1));
	t24 = t22 * t18;
	t23 = t22 * t19;
	t1 = [-t25, -t24, 0, 0, 0, 0, -t24, 0, 0, 0; t23, -t26, 0, 0, 0, 0, -t26, 0, 0, 0; 0, t19, 0, 0, 0, 0, t19, 0, 0, 0; t26, -t23, 0, 0, 0, 0, -t23, 0, 0, 0; -t24, -t25, 0, 0, 0, 0, -t25, 0, 0, 0; 0, -t18, 0, 0, 0, 0, -t18, 0, 0, 0; t22, 0, 0, 0, 0, 0, 0, 0, 0, 0; t21, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiR_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (80->15), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t29 = -qJ(2) - qJ(7) + pkin(15) - qJ(8);
	t27 = sin(t29);
	t30 = sin(qJ(1));
	t33 = t30 * t27;
	t31 = cos(qJ(1));
	t32 = t31 * t27;
	t28 = cos(t29);
	t26 = t31 * t28;
	t25 = t30 * t28;
	t1 = [t25, -t32, 0, 0, 0, 0, -t32, -t32, 0, 0; -t26, -t33, 0, 0, 0, 0, -t33, -t33, 0, 0; 0, -t28, 0, 0, 0, 0, -t28, -t28, 0, 0; t33, t26, 0, 0, 0, 0, t26, t26, 0, 0; -t32, t25, 0, 0, 0, 0, t25, t25, 0, 0; 0, -t27, 0, 0, 0, 0, -t27, -t27, 0, 0; t31, 0, 0, 0, 0, 0, 0, 0, 0, 0; t30, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
end